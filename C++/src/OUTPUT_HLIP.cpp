/**
 Xiaobin, Min
*/

#include <cassie_controllers/OUTPUT_HLIP.hpp>
#include <cassie_common_toolbox/motion_transition.hpp>
#include <cassie_common_toolbox/RadioButtonMap.hpp>
#include <cassie_controllers/feedforward_controller.hpp>
#include <cassie_common_toolbox/bezier_tools.hpp>
#include <roslib_utilities/ros_package_utilities.hpp>
#include <control_utilities/limits.hpp>
#include <cassie_common_toolbox/linear_interpolation.hpp>
#include <cassie_common_toolbox/CassieConstants.hpp>

using namespace control_utilities;
using namespace bezier_tools;
using namespace std;


// TODO: put another angular momentum version too
void OUTPUT_HLIP::init() {
    // init is the parames that has to be reset after reuse or intially use 

    updated.isDSP = true; // walking typically start from DSP 

    updated.tNstep0 = -1; // negative value. 
    updated.stanceLeg = leftStance;
    updated.ithStep = 0;
    updated.readyToStop = false;

    HLIP_sagittal.init(params.zNorm, params.Ts, params.Td, P1orbit, 0, 0);
    HLIP_laterral.init(params.zNorm, params.Ts, params.Td, P2orbit, 0, params.stepWidthNorminal);

    Vector2d PhaseRange;
    PhaseRange << 0, params.T;
    phase.reconfigure(PhaseRange, 1);
}

void OUTPUT_HLIP::reset() {
    // if the same class object will be used. 
    // do something 
    init();
}

void OUTPUT_HLIP::configure(ros::NodeHandle& nh) {
    // configure is something lives with this class that does not change at all
    params.paramChecker.init(nh.getNamespace() + "/stepping");

    ros::param::get("/cassie/is_simulation", params.isSim);

    params.paramChecker.checkAndUpdate("vx_offset", params.vx_offset);
    params.paramChecker.checkAndUpdate("vy_offset", params.vy_offset);

    params.paramChecker.checkAndUpdateYaml("bezierSwingHorizontal", params.bezierSwingHorizontal);
    params.paramChecker.checkAndUpdateYaml("bezierComVertical", params.bezierComVertical);

    params.paramChecker.checkAndUpdate("zsw_max", params.zsw_max);
    params.paramChecker.checkAndUpdate("zsw_neg", params.zsw_neg);

    params.paramChecker.checkAndUpdate("zNorm", params.zNorm);
    params.paramChecker.checkAndUpdate("stepWidthNorminal", params.stepWidthNorminal);
}

OUTPUT_HLIP::OUTPUT_HLIP(ros::NodeHandle& nh, cassie_model::Cassie& Robot) {
    robot = &Robot;
    configure(nh);
    init();
}

void OUTPUT_HLIP::updateOutputs(VectorXd& radio, double t) {
    // TODOOOOO:: check the SLRT implementation about where did I update waling params

    if (updated.tNstep0 < 0)
        updated.tNstep0 = t;

    cout << "sim Time: " << t << endl;
    phase.update(t - updated.tNstep0);

    logTime = t;
    timebasedSwitch(radio, t);

    computeActual();
    if (phase.tau == 0) {
        updated.y0 = updated.ya;
        updated.dy0 = updated.dya;
    }
    getCOMstate();
    updateAngMomentum();

    // Get the outputs  
    computeDesired();

    // Pull the radio outputs for movement
    updated.requestTransition = (radio(SB) < 1.0);
}

//time-based domain switching
void OUTPUT_HLIP::timebasedSwitch(VectorXd& radio, double t) {
    if (phase.tau == 1) {
        updated.readyToStop = true;
        updated.stanceLeg = (updated.stanceLeg == leftStance) ? rightStance : leftStance;
        updated.ithStep++;
        updateTargetWalking(radio);
        updated.tNstep0 = t;
        phase.tau = 0;
    }
}

void OUTPUT_HLIP::updateTargetWalking(VectorXd& radio) {

    // TODO: use a modifiable way to interpret the radio
    lowpassVelXdes.update(params.velXmax * radio(LV));
    lowpassVelYdes.update(params.velYmax * radio(LH));

    double desiredVx = lowpassVelXdes.getValue() + params.vx_offset;
    double desiredVy = lowpassVelYdes.getValue() + params.vy_offset;

    double stepWidth = params.stepWidthNorminal + 0.2 * radio(S2);
    stepWidth = clamp(stepWidth, 0.05, 0.5); // range of the realizable step Width

    params.T = params.Tnorm + 0.15 * radio(RS); // walking period 
    params.z0 = params.zNorm + 0.2 * radio(LS); // walking height

    //reconfigure phase
    Vector2d PhaseRange;
    PhaseRange << 0, this->params.T;
    phase.reconfigure(PhaseRange, 1);

    // update the H-LIP, currently simply scale the potential Tssp and Tdsp. 
    HLIP_sagittal.updateHLIP(params.z0, params.T / params.Tnorm * params.Ts, params.T / params.Tnorm * params.Td);
    HLIP_laterral.updateHLIP(params.z0, params.T / params.Tnorm * params.Ts, params.T / params.Tnorm * params.Td);

    // update the H-LIP target walking behaviors 
    HLIP_sagittal.updateDesiredWalking(desiredVx, 0);
    HLIP_laterral.updateDesiredWalking(desiredVy, stepWidth);

}

// get actual and desired output function
void OUTPUT_HLIP::computeActual() {
    // TODO: fix these bad conventions 

    // MatrixXd Dya = MatrixXd::Zero(9, nConfigSpace * 2);
    // MatrixXd DLfya = MatrixXd::Zero(9, nConfigSpace * 2);
    // if (updated.stanceLeg == rightStance) {
    //     VectorWrap ya_(updated.ya), dya_(updated.dya);
    //     SymFunction::yaRightStance_new(ya_, robot->q);
    //     SymFunction::dyaRightStance_new(dya_, robot->q, robot->dq);
    //     SymFunction::Dya_RightStanceActual_new(Dya, robot->q);
    //     SymFunction::DLfya_RightStanceActual_new(DLfya, robot->q, robot->dq);
    // }
    // else {
    //     VectorWrap ya_(updated.ya), dya_(updated.dya);
    //     SymFunction::yaLeftStance_new(ya_, robot->q);
    //     SymFunction::dyaLeftStance_new(dya_, robot->q, robot->dq);
    //     SymFunction::Dya_LeftStanceActual_new(Dya, robot->q);
    //     SymFunction::DLfya_LeftStanceActual_new(DLfya, robot->q, robot->dq);
    // }
    // updated.Jya << Dya.block(0, 0, 9, nConfigSpace);
    // updated.dJya << DLfya.block(0, 0, 9, nConfigSpace);


    if (updated.stanceLeg == rightStance) {
        VectorWrap ya_(updated.ya), dya_(updated.dya);
        SymFunction::yaRightStanceFull(ya_, robot->q);
        SymFunction::dyaRightStanceFull(dya_, robot->q, robot->dq);
        SymFunction::J_yaRightStanceFull(updated.Jya, robot->q);
        SymFunction::Jdot_yRightStanceFull(updated.dJya, robot->q, robot->dq);
    }
    else {
        VectorWrap ya_(updated.ya), dya_(updated.dya);
        SymFunction::yaLeftStanceFull(ya_, robot->q);
        SymFunction::dyaLeftStanceFull(dya_, robot->q, robot->dq);
        SymFunction::J_yaLeftStanceFull(updated.Jya, robot->q);
        SymFunction::Jdot_yLeftStanceFull(updated.dJya, robot->q, robot->dq);
    }

}

void OUTPUT_HLIP::computeDesired() {

    getSwingXdes();
    getSwingYdes();
    getSwingZdes();
    getComZdes();

    // pelvis RPY set to 0
    updated.ydes(y_deltaRoll_idx) = 0;
    updated.ydes(y_deltaPitch_idx) = 0;
    updated.ydes(y_stanceHipYaw_idx) = 0;
    updated.ydes(y_swingHipYaw_idx) = 0;
    updated.ydes(y_deltaSwingFoot_idx) = 0;

    updated.dydes(y_deltaRoll_idx) = 0;
    updated.dydes(y_deltaPitch_idx) = 0;
    updated.dydes(y_stanceHipYaw_idx) = 0;
    updated.dydes(y_swingHipYaw_idx) = 0;
    updated.dydes(y_deltaSwingFoot_idx) = 0;

    updated.ddydes(y_deltaRoll_idx) = 0;
    updated.ddydes(y_deltaPitch_idx) = 0;
    updated.ddydes(y_stanceHipYaw_idx) = 0;
    updated.ddydes(y_swingHipYaw_idx) = 0;
    updated.ddydes(y_deltaSwingFoot_idx) = 0;

    // update_turning_yaw();
}

void OUTPUT_HLIP::getSwingYdes() {

    // TODO: add a functionality to be able to use the predicted preimpact states 
    // predictPreImpact(this->phase.pActual, this->param.Ts);

    Vector2d step = HLIP_laterral.getDesiredStepSizeDeadbeat(updated.pCOM(y_in_xyz3),
        updated.vCOM(y_in_xyz3), updated.stanceLeg);

    cout << "comY : " << updated.pCOM(y_in_xyz3) << "   " << updated.vCOM(y_in_xyz3) << "    "
        "u_y_des : " << step(0) << endl;
    step(0) = (updated.stanceLeg == rightStance) ? clamp(step(0), 0.06, 0.6) : clamp(step(0), -0.6, -0.06);
    cout << "HLIP left des : " << HLIP_laterral.p2.XleftDes.transpose() << "    " <<
        "HLIP right des : " << HLIP_laterral.p2.XrightDes.transpose() << endl;

    // step(0) = (updated.stanceLeg == rightStance) ? HLIP_laterral.p2.UrightDes : HLIP_laterral.p2.UleftDes;
    // step(1) = 0;

    double bht, dbht, d2bht;
    bht = bezier(params.bezierSwingHorizontal, phase.tau);
    dbht = dtimeBezier(params.bezierSwingHorizontal, phase.tau, phase.dtau);
    d2bht = dtime2Bezier(params.bezierSwingHorizontal, phase.tau, phase.dtau);

    updated.ydes(y_swingStepy_idx) = (1. - bht) * updated.y0(y_swingStepy_idx) + bht * step(0);
    updated.dydes(y_swingStepy_idx) = -dbht * updated.y0(y_swingStepy_idx) + dbht * step(0) + bht * step(1);
    updated.ddydes(y_swingStepy_idx) = -d2bht * updated.y0(y_swingStepy_idx) + d2bht * step(0) + dbht * step(1) +
        dbht * step(1);// + bht*this->config.g/this->cache.z0LIP*this->cache.ycLIP;
}

void OUTPUT_HLIP::getSwingXdes() {

    Vector2d step = HLIP_sagittal.getDesiredStepSizeDeadbeat(updated.pCOM(x_in_xyz3),
        updated.vCOM(x_in_xyz3), updated.stanceLeg);

    cout << "comX : " << updated.pCOM(x_in_xyz3) << "   " << updated.vCOM(x_in_xyz3) << "    "
        "u_x_des : " << step(0) << endl;

    step(0) = clamp(step(0), -params.maxStepSize, params.maxStepSize);

    // step.setZero();

    double bht, dbht, d2bht;
    bht = bezier(params.bezierSwingHorizontal, phase.tau);
    dbht = dtimeBezier(params.bezierSwingHorizontal, phase.tau, phase.dtau);
    d2bht = dtime2Bezier(params.bezierSwingHorizontal, phase.tau, phase.dtau);

    updated.ydes(y_swingStepx_idx) = (1 - bht) * updated.y0(y_swingStepx_idx) + bht * step(0);
    updated.dydes(y_swingStepx_idx) = -dbht * updated.y0(y_swingStepx_idx) + dbht * step(0) + bht * step(1);
    updated.ddydes(y_swingStepx_idx) = -d2bht * updated.y0(y_swingStepx_idx) + d2bht * step(0) + dbht * step(1) + dbht * step(1);
}

void OUTPUT_HLIP::getSwingZdes() {
    //TODO: find a better way to do this
    VectorXd bv(8);
    bv << updated.y0(y_swingStepz_idx), updated.y0(y_swingStepz_idx), params.zsw_max / 3,
        params.zsw_max, params.zsw_max, params.zsw_max / 2, params.zsw_neg,
        params.zsw_neg;
    updated.ydes(y_swingStepz_idx) = bezier(bv, phase.tau);
    updated.dydes(y_swingStepz_idx) = dtimeBezier(bv, phase.tau, phase.dtau);
    updated.ddydes(y_swingStepz_idx) = dtime2Bezier(bv, phase.tau, phase.dtau);
}

void OUTPUT_HLIP::getComZdes() {
    double bht, dbht, d2bht;
    bht = bezier(params.bezierComVertical, this->phase.tau);
    dbht = dtimeBezier(params.bezierComVertical, phase.tau, phase.dtau);
    d2bht = dtime2Bezier(params.bezierComVertical, phase.tau, phase.dtau);

    updated.ydes(y_zCOM_idx) = (1 - bht) * updated.y0(y_zCOM_idx) + bht * params.z0;
    updated.dydes(y_zCOM_idx) = -dbht * updated.y0(y_zCOM_idx) + dbht * params.z0;
    updated.ddydes(y_zCOM_idx) = -d2bht * updated.y0(y_zCOM_idx) + d2bht * params.z0;
}

//turning controller
void OUTPUT_HLIP::update_turning_yaw() {
    // TODO: translate the SLRT implementaiotn here 
    double YawDes = 0;
    double K = 0.5;
    double yaw_error = robot->q(BaseRotZ) - YawDes;
    //    if (abs(yaw_error) > 0.05) {
    double ToYaw = K * yaw_error;
    updated.ydes(y_stanceHipYaw_idx) += ToYaw * bezier(params.bezierSwingHorizontal, phase.tau);
    updated.dydes(y_stanceHipYaw_idx) += ToYaw * dtimeBezier(params.bezierSwingHorizontal, phase.tau, phase.dtau);
    updated.ddydes(y_stanceHipYaw_idx) += ToYaw * dtime2Bezier(params.bezierSwingHorizontal, phase.tau, phase.dtau);

    //    }
    updated.ydes(y_stanceHipYaw_idx) += updated.y0(y_stanceHipYaw_idx) * bezier(params.bezierSwingHorizontal, 1.0 - phase.tau);
    updated.dydes(y_stanceHipYaw_idx) += updated.y0(y_stanceHipYaw_idx) * dtimeBezier(params.bezierSwingHorizontal, 1.0 - phase.tau, phase.dtau);
    updated.ddydes(y_stanceHipYaw_idx) += updated.y0(y_stanceHipYaw_idx) * dtime2Bezier(params.bezierSwingHorizontal, 1.0 - phase.tau, phase.dtau);
}

void OUTPUT_HLIP::getCOMstate() {

    MatrixXd p_COM(3, 1); //assumed pelvis yaw angle =0
    VectorXd qNoYaw(22);
    qNoYaw << robot->q;
    qNoYaw(5) = 0;
    if (updated.stanceLeg == rightStance) {
        SymFunction::p_com_RightStance(p_COM, qNoYaw);
    }
    else {
        SymFunction::p_com_LeftStance(p_COM, qNoYaw);
    }

    updated.pCOM << p_COM(0, 0), p_COM(1, 0), p_COM(2, 0);

    MatrixXd dp_COM(3, 1);
    SymFunction::dp_com_absolute(dp_COM, qNoYaw, robot->dq);

    updated.vCOM << dp_COM(0, 0), dp_COM(1, 0), dp_COM(2, 0);
}

void OUTPUT_HLIP::updateAngMomentum() {
    // update angular momentum about contact pivot
    Vector3d L3d;
    MatrixXd Lcom(6, 1); // centrodial
    VectorXd x(44, 1);
    x << robot->q, robot->dq;

    MatrixXd p_com_abs(3, 1);
    SymFunction::p_com_absolute(p_com_abs, robot->q);
    SymFunction::CentroidalMomentum(Lcom, x, p_com_abs);

    L3d = updated.pCOM.cross(updated.vCOM);
    updated.Lpivot << Lcom(3) / 33 + L3d(0),
        Lcom(4) / 33 + L3d(1);
}

Vector2d OUTPUT_HLIP::predictPreImpact(double time2impact, double p, double v) {
    // which HLIP not matter
    return HLIP_sagittal.get_LIPsol(time2impact, params.z0, p, v);
}

void OUTPUT_HLIP::getLog(VectorXf& log) {

    double tsec = static_cast<double>(ros::Time::now().sec);
    double tnsec = 1e-9 * static_cast<double>(ros::Time::now().nsec);
    // Move zero to closer time rather than 1970 so it fits in a float
    // Still bad practice to put time in floats, but this is just for logging
    if (!params.isSim)
        tsec -= 1631496319;//1.6074545e9;
    else
        tsec = logTime;

    // Use floats for logging size and speed
    log << static_cast<float>(tsec),
        static_cast<float>(tnsec), // 2
        static_cast<float>(updated.stanceLeg),      //1
        static_cast<float>(this->phase.tau),   // 1
        static_cast<float>(this->phase.dtau),  // 1
        updated.ya.cast<float>(),        // 9
        updated.dya.cast<float>(),       // 9
        updated.ydes.cast<float>(),        // 9
        updated.dydes.cast<float>(),       // 9
        updated.ddydes.cast<float>(),       //9 
        updated.pCOM.cast<float>(),//6
        updated.vCOM.cast<float>();  //
    //for walking 
}

