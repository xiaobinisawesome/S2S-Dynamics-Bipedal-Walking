#include <cassie_common_toolbox/CassieConstants.hpp>
#include <Eigen/Dense>


using namespace Eigen;

class HLIP {

public:

    HLIP();
    void init(double z0, double Ts, double Td, int orbitPeriod, double vel, double stepwidth);

    int orbitPeriod = 0; // 1 is P1, 2 is P2 

    struct Params
    {
        double z0, Ts, Td, T, lambda;
        double grav = 9.81;

        MatrixXd A_S2S = MatrixXd::Zero(2, 2);
        VectorXd B_S2S = VectorXd::Zero(2);

        double LIPacc = 0;
        double velDes = 0;
    } params;

    Vector2d Kdeadbeat;
    Vector2d X_HLIP;   // TODO: add two layers using a virtual HLIP

    void updateHLIP(double z0, double Ts, double Td);

    void updateDesiredWalking(double vel, double uLeftDes);

    void setDesiredVelocity(double vdes);

    Vector2d getDesiredStepSizeDeadbeat(double p, double v, int stanceLeg);

    struct P1 {
        // desired walking state for P1 orbits
        VectorXd Xdes = VectorXd::Zero(2);
        double Udes = 0;
        double sigma1 = 0;
        Vector2d Kdeadbeat;

        Vector2d getDeadbeatStepSize(double p, double v, double lambda);
    } p1;

    struct P2 {
        // desired walking state for P2 orbits
        VectorXd XleftDes = VectorXd::Zero(2);
        VectorXd XrightDes = VectorXd::Zero(2);
        double UleftDes = 0;
        double UrightDes = 0;
        double sigma2 = 0;
        double d2 = 0;
        Vector2d Kdeadbeat;

        Vector2d getDeadbeatStepSize(double p, double v, double lambda, int stanceLegIdx);
        void setP2orbitUleft(double uLeftDes);
    } p2;

    double solve_Ts();
    // RODO: fix the solution inconsistenc y (velocity or momentum?) . .
    Vector2d get_LIPsol(double t, double z0, double p0, double v0);

    double getOrbitalEnergy(double p, double Ly);

};