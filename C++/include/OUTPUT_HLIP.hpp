// Xiaobin, Min
#ifndef OUTPUT_HLIP_HPP
#define OUTPUT_HLIP_HPP

#include <cassie_description/cassie_model.hpp>
#include <ros_utilities/ros_utilities.hpp>
#include <cassie_common_toolbox/PhaseVariable.hpp>
#include <cassie_common_toolbox/smoothing.hpp>
#include <std_srvs/Empty.h>
#include <ros_utilities/timing.hpp>
#include <cassie_controllers/feedforward_controller.hpp>
#include <cassie_controllers/HLIP.hpp>
#include <control_utilities/filters.hpp>
#include <control_utilities/limits.hpp>
#include <cassie_common_toolbox/CassieConstants.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>


using namespace control_utilities;
using namespace Eigen;
// this class should be providing the desired outputs and contact idx to the controller based on time and current state of the robot 

class OUTPUT_HLIP {
public:

    OUTPUT_HLIP(ros::NodeHandle& nh, cassie_model::Cassie& robot);

    void updateOutputs(VectorXd& radio, double t); //update output
    void reset();
    void configure(ros::NodeHandle& nh);
    void init();
    bool isReadyToTransition() { return updated.readyToStop; }
    void getLog(VectorXf& dbg);

    // H-LIP class for planning step size
    HLIP HLIP_sagittal, HLIP_laterral;
    // Gait phasing variable
    PhaseVariable phase;

    double logTime = 0;
    struct Params {
        double grav = 9.81;
        bool isSim = true;
        // basic walking params, can be updated by joystick
        double Ts = 0.35, Td = 0.05, zNorm = 0.75; // nominal walking
        double zsw_max = 0.2, zsw_neg = -0.01;

        double Tnorm = Ts + Td;
        double T = Tnorm, z0 = zNorm;
        double stepWidthNorminal = 0.25;
        VectorXd bezierSwingHorizontal;
        VectorXd bezierComVertical;
        VectorXd bezierSwingVertical;

        // offset for hardware mainly
        double vx_offset = 0;
        double vy_offset = 0;

        // walking safety params 
        double maxStepSize = 1, velXmax = 2, velYmax = 0.5;

        // Parameter checker
        ros_utilities::ParamChecker paramChecker;
    } params;


    struct Updated {

        double tNstep0 = -1;
        //from SymFunction
        VectorXd ya = VectorXd::Zero(9);
        VectorXd dya = VectorXd::Zero(9);
        MatrixXd Jya = MatrixXd::Zero(9, nConfigSpace);
        MatrixXd dJya = MatrixXd::Zero(9, nConfigSpace);;

        // post-impact condition
        VectorXd y0 = VectorXd::Zero(9);
        VectorXd dy0 = VectorXd::Zero(9);

        //compute from bezier polynomial
        VectorXd ydes = VectorXd::Zero(9);
        VectorXd dydes = VectorXd::Zero(9);
        VectorXd ddydes = VectorXd::Zero(9);

        //stepping stone related
        double udes;
        double Ly;
        double Lx;

        // current robot COM state w.r.t. stance ankle
        Vector3d pCOM = Vector3d::Zero();
        Vector3d vCOM = Vector3d::Zero();
        Vector2d Lpivot = Vector2d::Zero();
        double z0LIP, dz0LIP, dz0LIP_f; // reserve for non-flat ground walking 

        double xcLIP_minus;
        double Ly_minus;
        double z0_minus;
        double dzc_minus;

        //contact domain
        int stanceLeg;
        int ithStep = 0;
        bool requestTransition = false;
        bool isDSP = true;
        bool readyToStop = false;
    } updated;


    struct LOGData
    { /// tempurarly for debugging quanties
        int doAddStuff;
    } logData;

private:

    // Pointer to the controlling nodehandle and related ROS things
    ros::NodeHandle* nh;

    // filters for desired walking beahviors from the joysticks 
    LowPassFilter lowpassVelXdes = LowPassFilter(0.001, 0.75);
    LowPassFilter lowpassVelYdes = LowPassFilter(0.001, 0.25);
    LowPassFilter lowpassVdesTurn = LowPassFilter(0.001, 0.25);

    // possibly being used for the pelvis velocities that come from the simulator; should put elsewhere
    LowPassFilter lowpassVelx = LowPassFilter(0.0005, 0.01);
    LowPassFilter lowpassVely = LowPassFilter(0.0005, 0.01);

    // Pointer to robot model
    cassie_model::Cassie* robot;

    // Private methods
    void timebasedSwitch(VectorXd& radio, double t);
    void computeActual();
    void computeDesired();

    void getSwingXdes();
    void getSwingYdes();
    void getSwingZdes();
    void getComZdes();

    void updateTargetWalking(VectorXd& radio);
    void update_turning_yaw();

    void getCOMstate();
    void updateAngMomentum();

    Vector2d predictPreImpact(double time2impact, double p, double v);

};


#endif // OUTPUT_HLIP_HPP
