#ifndef QPcontrolSOLVER_HPP
#define QPcontrolSOLVER_HPP

#include <cassie_description/cassie_model.hpp>
#include <ros_utilities/ros_utilities.hpp>
#include <cassie_common_toolbox/PhaseVariable.hpp>
#include <cassie_common_toolbox/smoothing.hpp>
#include <std_srvs/Empty.h>
#include <ros_utilities/timing.hpp>
#include <cassie_controllers/feedforward_controller.hpp>

#include <control_utilities/filters.hpp>
#include <control_utilities/limits.hpp>

#include <Eigen/Sparse>

#include <cassie_controllers/OUTPUT_HLIP.hpp>

#include <osqp/osqp.h>
#include <OsqpEigen/OsqpEigen.h> 

using namespace Eigen;

/////////////////////// Jan 13, keep opt variables and nCon constant across differen domains
class QPcontrolSolver {
public:

    QPcontrolSolver(ros::NodeHandle& nh, cassie_model::Cassie& robot);

    void getLog(VectorXf& dbg);
    void reset();
    bool getTorqueQP(VectorXd& u, OUTPUT_HLIP& output, double time);

private:

    ros_utilities::ParamChecker paramChecker;
    ros::NodeHandle* nh;
    cassie_model::Cassie* robot;    // Pointer to robot model
    OsqpEigen::Solver solver;   // osqp solver 


    void init(cassie_model::Cassie& Robot);

    bool isSim, isRigid;
    int stanceLeg = 0;  // -1 is left stance; 1 is right stance
    double logTime = 0;
    /////////////// single support walking QP controller. //////////////////////////

    int nOutputs = 9; // number of outputs
    int nJc = 2 + 4 + 6; // number of holonomic constraints 
    // 2 push rods, 4 spring joints, 6 single foot contact

    int idx = 0; // indexing the QP problems
    int maxQPIter = 50;
    double rho, alpha;
    bool warmstart, verbose, adaptRho;

    bool qp_initialized = false;

    ///////////////////// optimization variables /////////////////////////////
    struct optVar {
        int idx;
        int length;
    };

    struct optVar opt_ddq = { 0, 22 };
    struct optVar optU = { opt_ddq.idx + opt_ddq.length, 10 }; //22, 10
    // F = [F_rod;  F_leftFoot; F_rightFoot; F_springs;]
    struct optVar optF_rod = { optU.idx + optU.length, 2 }; // 32, 2
    struct optVar optLeftF = { optF_rod.idx + optF_rod.length, 6 }; // 34, 6
    struct optVar optRightF = { optLeftF.idx + optLeftF.length, 6 };  // 40, 6
    struct optVar optSpringF = { optRightF.idx + optRightF.length, 4 };
    int nVar = optSpringF.idx + optSpringF.length; //50

    ///////////////////// constraint matrix /////////////////////////////////////
    MatrixXd A_friction = MatrixXd::Zero(20, nVar); // two point contact will be 10 on each foot
    VectorXd b_lb_friction = VectorXd::Zero(20);
    VectorXd b_ub_friction = VectorXd::Zero(20);

    MatrixXd A_constContactForce = MatrixXd::Zero(6, nVar); // 6 contact force vectors to be zero 
    VectorXd b_constContactForce = VectorXd::Zero(6);

    int nCon = opt_ddq.length + A_friction.rows() + nJc + A_constContactForce.rows() + 1;
    //  EOM, 22  +  20  +  12  +  6
    //   60;  

    /////////////////////////////////////////////////////////////////
    //////////// solver related parameters //////////////////////////

    VectorXd u_prev = VectorXd::Zero(10);
    VectorXd u_sol = VectorXd::Zero(10);

    // QP Tuning
    double w_uChatter;
    double w_outputs;

    VectorXd OutputKP = VectorXd::Zero(nOutputs);
    VectorXd OutputKD = VectorXd::Zero(nOutputs); // pd on outputs;

    VectorXd JointKP = VectorXd::Zero(nTorques);
    VectorXd JointKD = VectorXd::Zero(nTorques);;  // pd on joints;

    // QP //  cost = 1x' H x + g' x; (Ax - b); A does not change, b updates
    MatrixXd H = MatrixXd::Zero(nVar, nVar);
    VectorXd g = VectorXd::Zero(nVar);

    MatrixXd A_holonomics = MatrixXd::Zero(nJc, nVar);
    VectorXd b_holonomics = VectorXd::Zero(nJc);
    MatrixXd Aeq_dynamics = MatrixXd::Zero(opt_ddq.length, nVar);
    VectorXd beq_dynamics = VectorXd::Zero(opt_ddq.length);

    MatrixXd A_chatter = MatrixXd::Zero(optU.length, nVar);
    VectorXd b_chatter = VectorXd::Zero(optU.length);
    MatrixXd A_y = MatrixXd::Zero(nOutputs, nVar);
    VectorXd b_y = VectorXd::Zero(nOutputs);
    MatrixXd A_springs = MatrixXd::Zero(4, nVar);
    VectorXd b_springs = VectorXd::Zero(4);
    MatrixXd A_reg = MatrixXd::Zero(nVar, nVar);
    VectorXd b_reg = VectorXd::Zero(nVar);

    MatrixXd A_toe = MatrixXd::Zero(1, nVar);

    MatrixXd AconstrLarge = MatrixXd::Zero(nCon + nVar, nVar);

    VectorXd sol = VectorXd::Zero(nVar);

    VectorXd VarLowerBounds = -OSQP_INFTY * VectorXd::Ones(nVar);
    VectorXd VarUpperBounds = OSQP_INFTY * VectorXd::Ones(nVar);
};















#endif //IDCLFQPSOLVER_HPP
