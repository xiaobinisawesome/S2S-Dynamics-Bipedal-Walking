//
// Created by dm on 8/12/21.
//
#include <cassie_controllers/QPcontrolSolver.hpp>
#include <cassie_common_toolbox/motion_transition.hpp>
#include <cassie_common_toolbox/RadioButtonMap.hpp>
#include <cassie_controllers/feedforward_controller.hpp>
#include <cassie_common_toolbox/bezier_tools.hpp>
#include <roslib_utilities/ros_package_utilities.hpp>
#include <control_utilities/limits.hpp>
#include <cassie_common_toolbox/linear_interpolation.hpp>

using namespace Eigen;
using namespace cassie_model;
using namespace control_utilities;
using namespace std;

QPcontrolSolver::QPcontrolSolver(ros::NodeHandle& nh, cassie_model::Cassie& robot) {
    this->robot = &robot;

    isRigid = false;
    ros::param::get("/cassie/if_rigid", isRigid);
    init(robot);
}


void QPcontrolSolver::init(cassie_model::Cassie& Robot) {

    // constamt matrix  || u - u_prev ||
    A_chatter.block(0, optU.idx, optU.length, optU.length) <<
        MatrixXd::Identity(optU.length, optU.length);

    // Friction cone // MUST MAKE SURE THE SIZE 
    MatrixXd Acone = Robot.PointContactConst;  // 10x6

    A_friction.block(0, optLeftF.idx, Acone.rows(), Acone.cols()) << Acone;  // left foot 
    A_friction.block(Acone.rows(), optRightF.idx, Acone.rows(), Acone.cols()) << Acone;            // right foot
    b_lb_friction.setConstant(-OSQP_INFTY); // -inf
    b_ub_friction.setZero();

    // Build all constraint matrices
    VarLowerBounds << -OSQP_INFTY * VectorXd::Ones(opt_ddq.length),           // ddq
        -Robot.MotorJointTorqueBounds,        // torque   hundreds 
        -OSQP_INFTY * VectorXd::Ones(optF_rod.length + optLeftF.length + optRightF.length + optSpringF.length);

    VarUpperBounds << OSQP_INFTY * VectorXd::Ones(opt_ddq.length),        // ddq
        Robot.MotorJointTorqueBounds,
        OSQP_INFTY* VectorXd::Ones(optF_rod.length + optLeftF.length + optRightF.length + optSpringF.length);

    paramChecker.init("/cassie/locomotion/stepping");
    // Get all QP tunable parameters
    paramChecker.checkAndUpdate("qp/rho", rho);
    paramChecker.checkAndUpdate("qp/alpha", alpha);
    paramChecker.checkAndUpdate("qp/warmstart", warmstart);
    paramChecker.checkAndUpdate("qp/adaptRho", adaptRho);
    paramChecker.checkAndUpdate("qp/verbose", verbose);
    paramChecker.checkAndUpdate("qp/maxQPIter", maxQPIter);

    paramChecker.checkAndUpdate("qp/w_u_chatter", w_uChatter);
    paramChecker.checkAndUpdate("qp/w_outputs", w_outputs);

    paramChecker.checkAndUpdateYaml("qp/OutputKP", OutputKP);
    paramChecker.checkAndUpdateYaml("qp/OutputKD", OutputKD);

    VectorXd oneLegKP(5), oneLegKD(5);
    paramChecker.checkAndUpdateYaml("qp/JointKP", oneLegKP);
    paramChecker.checkAndUpdateYaml("qp/JointKD", oneLegKD);

    JointKP << oneLegKP, oneLegKP;
    JointKD << oneLegKD, oneLegKD;

    // setting the QP solver
    solver.data()->setNumberOfVariables(nVar);
    int nSolverCon = nCon + nVar; // due to osqp's convention
    solver.data()->setNumberOfConstraints(nSolverCon);

    Eigen::SparseMatrix<double> Hcost_sparse(nVar, nVar);
    Hcost_sparse.setIdentity();
    VectorXd gCost = VectorXd::Zero(nVar);
    solver.data()->setHessianMatrix(Hcost_sparse);
    solver.data()->setGradient(gCost);

    MatrixXd Aconstr(nSolverCon, nVar);
    Aconstr << MatrixXd::Identity(nSolverCon, nVar);

    VectorXd lbAconstr(nSolverCon);
    lbAconstr.setConstant(-OSQP_INFTY);
    VectorXd ubAconstr(nSolverCon);
    ubAconstr.setConstant(OSQP_INFTY);

    Eigen::SparseMatrix<double> AconstrSparse = Aconstr.sparseView();
    solver.data()->setLinearConstraintsMatrix(AconstrSparse);
    solver.data()->setLowerBound(lbAconstr);
    solver.data()->setUpperBound(ubAconstr);

    solver.settings()->setVerbosity(verbose);
    solver.settings()->setPrimalInfeasibilityTollerance(1e-5);
    solver.settings()->setDualInfeasibilityTollerance(1e-5);
    solver.settings()->setMaxIteration(maxQPIter);
    solver.settings()->setAlpha(alpha);
    solver.settings()->setRho(rho);
    solver.settings()->setAdaptiveRho(adaptRho);
    solver.settings()->setWarmStart(warmstart);

    bool initSolve = solver.initSolver();
    solver.solve();
    // TODO: still fails in the first QP, 
}


void QPcontrolSolver::reset() {
    idx = 0;
    u_prev.setZero();
    u_sol.setZero();
    g.setZero();
    qp_initialized = false;
}

bool QPcontrolSolver::getTorqueQP(VectorXd& u, OUTPUT_HLIP& output, double time) {
    // Formulate the QP: States are X = [ddq; u; F] 
    logTime = time;
    // Task space terms
    // || A_y var - b_y ||
    //  Jy ddq + dJy dq = ddy_des - Kp y - Kd dy (y = y^a -y^d);

       // Output error
    VectorXd eta(nOutputs * 2);
    eta << output.updated.ya - output.updated.ydes,
        output.updated.dya - output.updated.dydes;

    stanceLeg = output.updated.stanceLeg;

    cout << "----------------------------------------" << endl;
    cout << "stanceLeg: " << stanceLeg << endl;
    cout << "output time: " << output.updated.tNstep0 << endl;
    cout << "tau:       " << output.phase.tau << endl;
    cout << "actual outputs: " << output.updated.ya.transpose() << endl;
    cout << "desird outputs: " << output.updated.ydes.transpose() << endl;
    cout << "----------------------------------------" << endl;

    //////////////////////// Task space ///////////////////////
    A_y.block(0, 0, nOutputs, opt_ddq.length) << output.updated.Jya;  // Jy
    b_y << output.updated.ddydes
        - output.updated.dJya * robot->dq  //
        - OutputKP.cwiseProduct(eta.head(nOutputs))
        - OutputKD.cwiseProduct(eta.tail(nOutputs));

    // Compute the robot constraints
    MatrixXd J_hol(nJc, 22), dJ_hol(nJc, 22);
    MatrixXd J_leftC(6, 22), J_rightC(6, 22);
    J_leftC << robot->kinematics.cache.J_leftToe,
        robot->kinematics.cache.J_leftHeel;
    J_rightC << robot->kinematics.cache.J_rightToe,
        robot->kinematics.cache.J_rightHeel;

    switch (stanceLeg)
    {
    case leftStance:
        J_hol << robot->kinematics.cache.J_achilles,
            robot->kinematics.J_springs,
            J_leftC;
        dJ_hol << robot->kinematics.cache.Jdot_achilles,
            robot->kinematics.Jdot_springs,
            robot->kinematics.cache.dJ_leftToe,
            robot->kinematics.cache.dJ_leftHeel;
        A_constContactForce.setZero();
        A_constContactForce.block(0, optRightF.idx, optRightF.length, optRightF.length) <<
            MatrixXd::Identity(optRightF.length, optRightF.length);
        A_toe.setZero();
        A_toe(0, opt_ddq.length + motorLeftFootPitch) = 1;
        break;
    case rightStance:
        J_hol << robot->kinematics.cache.J_achilles,
            robot->kinematics.J_springs,
            J_rightC;
        dJ_hol << robot->kinematics.cache.Jdot_achilles,
            robot->kinematics.Jdot_springs,
            robot->kinematics.cache.dJ_rightToe,
            robot->kinematics.cache.dJ_rightHeel;
        A_constContactForce.setZero();
        A_constContactForce.block(0, optLeftF.idx, optLeftF.length, optLeftF.length) <<
            MatrixXd::Identity(optLeftF.length, optLeftF.length);
        A_toe.setZero();
        A_toe(0, opt_ddq.length + motorRightFootPitch) = 1;

        break;
    default:
        ROS_ERROR("stance Leg is wrong");
        break;
    }

    A_holonomics.block(0, 0, nJc, opt_ddq.length) << J_hol;
    b_holonomics << -dJ_hol * robot->dq;

    // H*ddq - B*u - J^T F (Fc and Fs ) = - C
    Aeq_dynamics.block(0, 0, opt_ddq.length, opt_ddq.length) << robot->dynamics.H;
    Aeq_dynamics.block(0, optU.idx, opt_ddq.length, optU.length) << -robot->Be1;
    Aeq_dynamics.block(0, optF_rod.idx, opt_ddq.length, optF_rod.length) << -robot->kinematics.cache.J_achilles.transpose();
    Aeq_dynamics.block(0, optSpringF.idx, opt_ddq.length, optSpringF.length) << -robot->kinematics.J_springs.transpose();

    Aeq_dynamics.block(0, optLeftF.idx, opt_ddq.length, optLeftF.length) << -J_leftC.transpose();
    Aeq_dynamics.block(0, optRightF.idx, opt_ddq.length, optRightF.length) << -J_rightC.transpose(); // it's okay, for single stance, the forces on the other foot will be zero
    beq_dynamics << -robot->dynamics.C;  // shin and heel springs 

    // Torque chatter cost
    b_chatter << u_prev;

    // play with some combinations:e.g.: putting constriants in the cost for relaxation
    /*
        Aeq_dynamic = b, A_holonomics  = b , A_chatter = b, A_y = b, A_springs = b, A_friction < b
    */
    ///////////////////////////////////////// constraints ///////////////////////////////////////////////////////////////////////
    int Nconst = Aeq_dynamics.rows() + // 22
        A_friction.rows() +  // 20
        A_holonomics.rows() + //12
        A_constContactForce.rows() + 1; // 6

    if (Nconst != nCon) {
        cout << "number of constraints are not correct";
        return false;
    }

    MatrixXd Aconstr(Nconst, nVar);
    VectorXd lbAconstr(Nconst);
    VectorXd ubAconstr(Nconst);
    Aconstr << Aeq_dynamics,
        A_constContactForce,
        A_friction,
        A_holonomics,
        A_toe;

    lbAconstr << beq_dynamics,
        b_constContactForce,
        b_lb_friction,
        b_holonomics, 0;

    ubAconstr << beq_dynamics,
        b_constContactForce,
        b_ub_friction,
        b_holonomics, 0;

    ////////////////////////////////////////// cost setup ///////////////////////////////////////////////////////////////////
    // Construct the cost
    // 0.5*x'*H*x + g'*x

    VectorXd W_OUTPUTS(nOutputs), W_uCHATTER(10);
    W_uCHATTER << w_uChatter * VectorXd::Ones(10);
    W_OUTPUTS << w_outputs * VectorXd::Ones(nOutputs);

    // Build the QP cost formulation
    long n_cost_terms = W_uCHATTER.size() +
        W_OUTPUTS.size();
    /// || A x - b||_w
    MatrixXd A(n_cost_terms, nVar); // Gets passed directly to
    VectorXd b(n_cost_terms);
    VectorXd w(n_cost_terms);

    A << A_chatter, A_y;

    b << b_chatter, b_y;

    w << W_uCHATTER, W_OUTPUTS;

    for (int i = 0; i < w.size(); i++) {
        b(i) *= w(i);
        A.row(i) *= w(i);
    }

    // Need to declare row major for all matrices getting mapped to normal real_t arrays in qpoases.
    MatrixXd G(nVar, nVar);
    VectorXd g(nVar);
    G << A.transpose() * A;
    g << -A.transpose() * b;

    // config for OSQP    
    MatrixXd AconstrLarge = MatrixXd::Zero(nCon + nVar, nVar);
    VectorXd lbALarge = VectorXd::Zero(nCon + nVar);
    VectorXd ubALarge = VectorXd::Zero(nCon + nVar);

    AconstrLarge << Aconstr, MatrixXd::Identity(nVar, nVar);
    // for the lower and upper bounds of the opt variables
    lbALarge << lbAconstr, VarLowerBounds;
    ubALarge << ubAconstr, VarUpperBounds;

    // Convert Matrix to SparseMatrix
    Eigen::SparseMatrix<double> Gsparse = G.sparseView();
    Eigen::SparseMatrix<double> Aconstrsparse = AconstrLarge.sparseView();
    // update the QP
    solver.updateHessianMatrix(Gsparse);
    solver.updateGradient(g);

    solver.updateLinearConstraintsMatrix(Aconstrsparse);
    // lbALarge = lbALarge - 1e-5 * VectorXd::Ones(lbALarge.size());
    solver.updateBounds(lbALarge, ubALarge);

    if (!solver.isInitialized())
        solver.initSolver();

    int ifSolved = solver.solve();    // Solve the QP
    // Get solution
    if (ifSolved == 1) {
        sol.setZero();
        sol << solver.getSolution();
        u_sol << sol.segment(optU.idx, optU.length);
        u_prev << u_sol;
        // cout << "----------------------------------------" << endl;
        cout << "OSQP solution: " << endl;
        cout << "ddq: " << sol.segment(opt_ddq.idx, opt_ddq.length).transpose() << endl;
        cout << "u      : " << sol.segment(optU.idx, optU.length).transpose() << "    ";
        cout << "F rod  : " << sol.segment(optF_rod.idx, optF_rod.length).transpose() << endl;
        cout << "F left : " << sol.segment(optLeftF.idx, optLeftF.length).transpose() << "     ";
        cout << "F right: " << sol.segment(optRightF.idx, optRightF.length).transpose() << "     ";
        cout << "Fspring :" << sol.segment(optSpringF.idx, optSpringF.length).transpose() << endl;
        // cout << "c al Fs  : " << robot->kinematics.cache.Fsprings.transpose() << endl;
        cout << "----------------------------------------" << endl;
        cout << "OSQP cost terms: " << endl;
        VectorXd errorV(nOutputs), weightedErrorV(nOutputs);
        errorV << A_y * sol - b_y;
        weightedErrorV << W_OUTPUTS.cwiseProduct(errorV);
        cout << "cost:          " << 0.5 * sol.transpose() * Gsparse * sol + g.transpose() * sol << "    ";
        cout << "error vector Norm:  " << errorV.transpose() * errorV << endl;
        cout << "eta vec:  " << eta.transpose() << endl;
        cout << "error vector:  " << errorV.transpose() << endl;
        // cout << "error W Output:" << weightedErrorV.transpose() * weightedErrorV << endl;
        cout << "----------------------------------------" << endl;
        u << u_sol;
    }
    else {//(ifSolved != 1) {
        ROS_WARN("THE QP in locomotion FAILED!");
        // TODO: safe action 
        // u_sol = u_prev;
        u_sol.setZero();
        // u_prev.setZero();
    }

    return (ifSolved == 1);
}

void QPcontrolSolver::getLog(VectorXf& log) {
    // Use floats for logging size and speed

    log << static_cast<float>(logTime),
        sol.cast <float>();
}


