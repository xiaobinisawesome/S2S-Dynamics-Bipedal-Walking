

#include <cassie_controllers/HLIP.hpp>
#include <ros_utilities/ros_utilities.hpp>
#include <cassie_common_toolbox/hyperbolic.hpp>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

HLIP::HLIP(){};

void HLIP::init(double z0, double Ts, double Td, int orbitPeriod, double vel, double stepwidth = 0)
{

    this->orbitPeriod = orbitPeriod;
    updateHLIP(z0, Ts, Td);
    updateDesiredWalking(vel, stepwidth);
}

void HLIP::updateHLIP(double z0, double Ts, double Td)
{
    params.z0 = z0;
    params.Ts = Ts;
    params.Td = Td;

    params.T = Ts + Td;
    params.lambda = sqrt(params.grav / params.z0);

    // S2S dynamics
    MatrixXd eATs(2, 2);
    eATs << cosh(Ts * params.lambda), sinh(Ts * params.lambda) / params.lambda,
        params.lambda * sinh(Ts * params.lambda), cosh(Ts * params.lambda);

    VectorXd temp2(2);
    MatrixXd temp22(2, 2);
    temp2 << -1, 0;
    temp22 << 1, Td, 0, 1;
    params.B_S2S = eATs * temp2;
    params.A_S2S = eATs * temp22;

    Kdeadbeat << 1, params.Td + coth(params.Ts * params.lambda) / params.lambda;

    // Orbit slopes
    p1.sigma1 = params.lambda * coth(params.Ts / 2 * params.lambda);
    p2.sigma2 = params.lambda * tanh(params.Ts / 2 * params.lambda);

    p1.Kdeadbeat = Kdeadbeat;
    p2.Kdeadbeat = Kdeadbeat;
}

void HLIP::updateDesiredWalking(double vel, double stepWidth)
{
    // only p2 orbit need to update step width
    if (orbitPeriod == P2orbit)
        p2.setP2orbitUleft(-stepWidth); // left stance u is negative typically
    setDesiredVelocity(vel);
}

void HLIP::setDesiredVelocity(double vel)
{
    params.velDes = vel;

    double DistSum = vel * params.T; // total traveled distance in SSPand DSP.

    switch (orbitPeriod)
    {
    case P1orbit:
        double vxdes;
        vxdes = DistSum / (2 / p1.sigma1 + params.T);
        p1.Xdes << vxdes / p1.sigma1, vxdes;
        p1.Udes = p1.Xdes(0) * 2 + p1.Xdes(1) * params.Td;
        break;
    case P2orbit:
        p2.d2 = pow(params.lambda, 2) * pow((sech(params.lambda * params.Ts / 2)), 2) * params.velDes *
                (params.Ts + params.Td) / (pow(params.lambda, 2) * params.Td + 2 * p2.sigma2);

        double pleftdes, vleftdes;
        pleftdes = (p2.UleftDes - params.Td * p2.d2) / (2 + params.Td * p2.sigma2);
        vleftdes = p2.sigma2 * pleftdes + p2.d2;

        p2.XleftDes << pleftdes, vleftdes;

        // calculate the right stance base boundary states
        p2.UrightDes = DistSum * 2 - p2.UleftDes;
        double prightdes, vrightdes;
        prightdes = (p2.UrightDes - params.Td * p2.d2) / (2 + params.Td * p2.sigma2);
        vrightdes = p2.sigma2 * prightdes + p2.d2;
        p2.XrightDes << prightdes, vrightdes;

        break;
    default:
        ROS_WARN("orbit type is wrong!");
        break;
    }
}

void HLIP::P2::setP2orbitUleft(double uLeftDes)
{
    this->UleftDes = uLeftDes;
}

Vector2d HLIP::getDesiredStepSizeDeadbeat(double p, double v, int stanceLeg)
{
    Vector2d out = (orbitPeriod == P1orbit) ? p1.getDeadbeatStepSize(p, v, params.lambda) : p2.getDeadbeatStepSize(p, v, params.lambda, stanceLeg);
    return out;
};

Vector2d HLIP::P1::getDeadbeatStepSize(double p, double v, double lambda)
{
    Vector2d Xnow, dXnow, Ustepsize;
    Xnow << p, v;
    dXnow << v, pow(lambda, 2) * p;
    double dstepLength = 0;

    double stepLength = Kdeadbeat.transpose() * (Xnow - this->Xdes) + this->Udes;

    dstepLength = Kdeadbeat.transpose() * dXnow;

    Ustepsize << stepLength, dstepLength;
    return Ustepsize;
}

Vector2d HLIP::P2::getDeadbeatStepSize(double p, double v, double lambda, int stanceLegIdx)
{

    Vector2d Xnow, dXnow, Ustepsize;
    Xnow << p, v;
    dXnow << v, pow(lambda, 2) * p;
    double dstepLength = 0;
    double stepLength = 0;

    switch (stanceLegIdx)
    {
    case leftStance: // left stance
        stepLength = Kdeadbeat.transpose() * (Xnow - XleftDes) + UleftDes;
        break;
    case rightStance:
        stepLength = Kdeadbeat.transpose() * (Xnow - XrightDes) + UrightDes;
        break;
    default:
        ROS_WARN("leg idx is wrong!");
        break;
    }
    dstepLength = Kdeadbeat.transpose() * dXnow;
    Ustepsize << stepLength, dstepLength;
    return Ustepsize;
}

///////////////////////////////// other LIP related helper functions ////////////////////////
double HLIP::solve_Ts()
{
    double Ts = 0;

    ///////// TODO: Time to impact.
    // double xcdes = this->param.xratio * this->param.ldes;
    // double lam = sqrt(this->config.g / this->cache.z0LIP);
    // double a = this->cache.xcLIP;
    // double b = this->cache.Ly / lam / this->cache.z0LIP;
    // double c = xcdes;
    // Ts = 1 / lam * log((c + sqrt(-pow(a, 2) + pow(b, 2) + pow(c, 2))) / (a + b));

    return Ts;
}

Vector2d HLIP::get_LIPsol(double t, double z0, double x0, double v0)
{

    Vector2d sol;

    MatrixXd A(2, 2);
    A << 0, 1,
        pow(params.lambda, 2), 0;

    MatrixXd At = A * t;
    Vector2d IC;
    IC << x0, v0;
    sol = At.exp() * IC;
    return sol;
}

double HLIP::getOrbitalEnergy(double p, double Ly)
{
    // get the orbital energy in SSP
    // p is the postion w.r.t. contact pivot, Ly is the momentum about the contact pivot
    return pow(Ly / params.z0, 2) - params.grav /
                                        params.z0 * pow(p, 2);
}