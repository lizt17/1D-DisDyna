#include <vector>
#include <cmath>

struct mobilityLaw_W {
    // Constants
    const double kB_eV = 8.6173303e-5; // eV/K
    const double Tm = 3695.0;
    const double a_SI = 3.16e-10;
    const double b_SI = a_SI * std::sqrt(3.0) / 2.0;
    const double mu_SI = 161e9;
    const double h_SI = a_SI * std::sqrt(2.0 / 3.0);
    const double w_SI = a_SI * 25.0;
    const double rho_SI = 19250.0;
    const double cs_SI = std::sqrt(mu_SI / rho_SI);

    // Dislocation parameters
    const double dH0 = 1.63;
    const double p = 0.86;
    const double q = 1.69;
    const double T0 = 0.8 * Tm;
    const double tauP_SI = 2.03e9;
    const double a0 = 1.5;
    const double Bk_SI = 8.3e-5;
    const double B0_SI = 9.8e-4;
    const double L_SI = 1e-7;

    // Normalized constants
    const double a = a_SI / b_SI;
    const double h = h_SI / b_SI;
    const double L = L_SI / b_SI;
    const double w = w_SI / b_SI;
    const double tauP = tauP_SI / mu_SI;
    const double Bk = Bk_SI * cs_SI / (mu_SI * b_SI);
    const double B0 = B0_SI * cs_SI / (mu_SI * b_SI);
    const double factor1 = 2.0 * a;
    const double factor2 = 1.0 / (2.0 * h * L);
    const double inv2kB = 1.0 / (2.0 * kB_eV);

    static double sigmoid(const double & x)
    {
        return 2.0/(1.0+exp(2.0*x));
    }
    std::pair<bool, double> velocity(double rss, double Temp) const
    {
        double tau_abs = std::abs(rss);
        int sign = (rss>0.0)? 1.0 : -1.0;
        double Theta = tau_abs / tauP / a0;
        const double dGkp = (Theta<1.0)? (std::pow(1.0-std::pow(Theta,p),q)-Temp/T0) : 0.0;
        const double dGkp1 = (dGkp>0.0)? dGkp : 0.0;
        const double expCoeff = exp(-dH0*dGkp1/(2.0*kB_eV*Temp));
        // Compute screw drag coeff
        const double sgm=0.5*sigmoid(-0.5*(0.05-dGkp1)/0.05);
        const double Bs=Bk*w/(2.0*h)*(1.0-sgm)+B0*sgm; //kink-dominated to drag-dominated interpolation
        // Compute screw velocity
        double vs=std::fabs(tau_abs)/Bs*expCoeff;
        bool athermal = dGkp <= 0.0;

        return std::make_pair(athermal, vs*sign);
    }

    std::pair<bool, double> velocity_Q(double rss, double Temp)
    {
        double v0_SI = 5.8e5; // [m/s]
        double Q_eV = 1.75;
        double VA = 20.0;
        double v0 = v0_SI / cs_SI;

        double dH = Q_eV - rss * VA;
        dH = dH<=0.0 ? 0.0 : dH;
        double vel = v0 * std::exp(-dH / (kB_eV * Temp));
        return std::make_pair(false, vel);

    }

    std::pair<bool, double> velocity_Inf(double rss, double Temp)
    {
        double B_inf_SI = 1e8;
        double B_inf = B_inf_SI * cs_SI / (mu_SI * b_SI);
        double v_inf = std::abs(rss) / B_inf;
        return std::make_pair(false, v_inf);

    }

    const double dG_kinkpair(double _rss, double _Temp)
    {
        double tau_abs = std::abs(_rss);
        double Theta = tau_abs / tauP / a0;
        const double dGkp = (Theta<1.0)? (std::pow(1.0-std::pow(Theta,p),q)-_Temp/T0) : 0.0;
        const double dGkp1 = (dGkp>0.0)? dGkp : 0.0;

        return dGkp1;
    }

};
