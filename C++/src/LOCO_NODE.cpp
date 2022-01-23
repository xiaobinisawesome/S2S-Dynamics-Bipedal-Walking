//
// Xiaobin, Min
//

#include <ros/ros.h>

// Robot descriptions
#include <cassie_description/cassie_model.hpp>
#include <cassie_common_toolbox/RadioButtonMap.hpp>
#include <cassie_common_toolbox/CassieConstants.hpp>

// Controllers
#include <cassie_controllers/standingControlTSC.hpp>
#include <cassie_controllers/inairControlTSC.hpp>
#include <cassie_controllers/QPcontrolSolver.hpp>
#include <cassie_controllers/OUTPUT_HLIP.hpp>

// Load in message types
#include <cassie_common_toolbox/callback_helpers.hpp>
#include <cassie_common_toolbox/cassie_control_msg.h>
#include <cassie_common_toolbox/cassie_proprioception_msg.h>
#include <cassie_common_toolbox/cassie_velocity_estimation_msg.h>

#include <realtime_utilities/timing.hpp>

// Logging
#include <fstream>
#include <cmath>

using namespace Eigen;
using namespace cassie_model;

// Robotic states
static bool is_initialized = false;
static Cassie robot;
static VectorXd radio(16);
static int mode_command = -1; // Coming from the radio
static double simTime = 0;


// get robot states 
void proprioception_callback(const cassie_common_toolbox::cassie_proprioception_msg::ConstPtr& propmsg)
{
    unpack_proprioception(propmsg, robot.q, robot.dq, radio, robot.gyroscope, robot.accelerometer, robot.torque,
        robot.leftContact, robot.rightContact);
    get_proprioception_orientation(*propmsg, robot.q, robot.dq, robot.quat_pelvis);
    mode_command = (int)radio(RadioButtonMap::SB);
    is_initialized = propmsg->isSimRunning && propmsg->isCalibrated;
    simTime = propmsg->simTime;
}

// get the velocity estimation 
void velocity_estimation_callback(const cassie_common_toolbox::cassie_velocity_estimation_msg::ConstPtr& velmsg)
{
    robot.dq(BasePosX) = velmsg->linear_velocity.x;
    robot.dq(BasePosY) = velmsg->linear_velocity.y;
    robot.dq(BasePosZ) = velmsg->linear_velocity.z;
}

// Main node
int main(int argc, char* argv[])
{
    // Establish the current ROS node and associated timing
    ros::init(argc, argv, "loco_node");
    ros::NodeHandle nh("/cassie/locomotion");

    // get some parameters
    double dt_des;  // control frequency 
    ros::param::get("/cassie/locomotion/dt", dt_des);

    bool is_one_domain = true;
    ros::param::get("/cassie/locomotion/use_one_domain", is_one_domain);
    ros::param::get("/cassie/locomotion/if_rigid", robot.kinematics.ifRigid);

    // Setup ROS publisher/subscriber networks
    ros::Rate looprate(1 / dt_des);
    ros::Publisher controller_pub = nh.advertise<cassie_common_toolbox::cassie_control_msg>("/cassie_control", 1);
    ros::Subscriber vel_est_sub = nh.subscribe("/cassie_velocity_estimation", 1, velocity_estimation_callback, ros::TransportHints().tcpNoDelay(true));
    ros::Subscriber proprioception_sub = nh.subscribe("/cassie_proprioception", 1, proprioception_callback, ros::TransportHints().tcpNoDelay(true));

    cassie_common_toolbox::cassie_control_msg control_message;

    // controller objects
    StandingControlTSC stand_control(nh, robot);

    OUTPUT_HLIP outputWalking(nh, robot);
    QPcontrolSolver qpControllerL(nh, robot), qpControllerR(nh, robot);

    // log files
    std::fstream logfileStand, logfileOutput, logfileQP;

    bool log_controller = false;
    ros::param::get("/cassie/log_controller", log_controller);
    if (log_controller)
    {
        std::string home = getenv("HOME");

        std::string path = home + "/ROBOTLOG/CASSIE/OutputWalk_log.bin";
        logfileOutput.open(path, std::ios::out | std::ios::binary);

        path = home + "/ROBOTLOG/CASSIE/qpSol_walk_log.bin";
        logfileQP.open(path, std::ios::out | std::ios::binary);

        path = home + "/ROBOTLOG/CASSIE/stand_log.bin";
        logfileStand.open(path, std::ios::out | std::ios::binary);
    }

    int cur_mode = locomotion_null;
    VectorXd u = VectorXd::Zero(nTorques);

    while (ros::ok()) {

        ros::spinOnce();

        // Zero the torque and set the message timestamp
        control_message.motor_torque.fill(0);

        if (is_initialized) {        // Main state machine

            // update kinematics and dynamics calculation
            robot.kinematics.update(robot.model, robot.q, robot.dq);
            if (robot.kinematics.ifRigid) {
                robot.q = robot.kinematics.q_rigid;
                robot.dq = robot.kinematics.dq_rigid;
            }
            robot.dynamics.calcHandC(robot.model, robot.q, robot.dq);

            // based on the current mode do the controllers
            switch (cur_mode) {
            case locomotion_null: {
                // do nothing
                break;
            }
            case locomotion_stand: {
                stand_control.updateTime(simTime);
                stand_control.updateControl(radio, u);
                if (log_controller) {
                    VectorXf standLog;
                    standLog = stand_control.getLog();
                    logfileStand.write(reinterpret_cast<char*>(standLog.data()),
                        (standLog.size()) * sizeof(float));
                }
                break;
            }
            case locomotion_walk: {
                outputWalking.updateOutputs(radio, simTime);
                bool qpsolved = false;
                if (outputWalking.updated.stanceLeg == leftStance)
                    qpControllerL.getTorqueQP(u, outputWalking, simTime);
                else
                    qpControllerR.getTorqueQP(u, outputWalking, simTime);
                // if (!qpsolved)  
                    // do a id controller
                if (log_controller) {
                    VectorXf qpSolLog, outputLog;
                    outputLog = outputWalking.getLog();
                    if (outputWalking.updated.stanceLeg == leftStance)
                        qpSolLog = qpControllerL.getLog();
                    else
                        qpSolLog = qpControllerR.getLog();
                    logfileOutput.write(reinterpret_cast<char*>(outputLog.data()), (outputLog.size()) * sizeof(float));
                    logfileQP.write(reinterpret_cast<char*>(qpSolLog.data()), (qpSolLog.size()) * sizeof(float));
                }
                break;
            }
            default:
                ROS_ERROR("STATE MACHINE IS WRONG");
                return 1;
            }

            if (cur_mode != mode_command) {
                // check switching 
                switch (cur_mode) {
                case locomotion_null: {
                    if (mode_command == locomotion_stand) {
                        ROS_INFO("Transitioning to standing  !");
                        stand_control.reset();
                    }
                    break;
                }
                case locomotion_stand: {
                    if (mode_command == locomotion_null) {
                        ROS_INFO("Transitioning to null!");
                    }
                    if ((mode_command == locomotion_walk) && (stand_control.isReadyToTransition())) {
                        ROS_INFO("Transitioning to walking  !");
                        outputWalking.reset();
                        qpControllerL.reset();
                        qpControllerR.reset();
                    }
                    break;
                }
                case locomotion_walk: {
                    if ((mode_command == locomotion_stand) && (outputWalking.isReadyToTransition())) {
                        ROS_INFO("Transitioning to standing control!");
                        stand_control.reset();
                    }
                    if (mode_command == locomotion_null) {
                        ROS_INFO("Transitioning to null control!");
                    }
                    break;
                }
                default:
                    ROS_ERROR("STATE MACHINE IS WRONG");
                    return 1;
                }
                cur_mode = mode_command;
            }
        }

        // process control torques
        for (int i = 0; i < 10; ++i) {
            bool val = std::isnan(control_message.motor_torque[i]);
            if (val) {        // check if NAN torques
                control_message.motor_torque.fill(0.0);
                control_message.header.stamp = ros::Time::now();
                controller_pub.publish(control_message);

                ROS_ERROR("The controller generated a nan torque!!!!!!!");
                return 1;
            }
            control_message.motor_torque[i] = u(i);
        }

        // Publish the control message
        control_message.header.stamp = ros::Time::now();
        controller_pub.publish(control_message);

        looprate.sleep();
    }

    if (log_controller)
    {
        logfileOutput.close();
        logfileQP.close();
        logfileStand.close();
    }
    ros::shutdown();
    return 0;
}
