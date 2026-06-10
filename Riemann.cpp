#include "Riemann.H"

#include <algorithm>
#include <cmath>

namespace
{
    constexpr double gam = 1.4;
    constexpr double TOL = 1.e-6;
    constexpr int NRITER = 20;
    constexpr double xdisc = 0.5;

    constexpr double g[9] = {
        gam,
        (gam - 1.0) / (2.0 * gam),
        (gam + 1.0) / (2.0 * gam),
        2.0 * gam / (gam - 1.0),
        2.0 / (gam - 1.0),
        2.0 / (gam + 1.0),
        (gam - 1.0) / (gam + 1.0),
        0.5 * (gam - 1.0),
        gam - 1.0
    };

    double guessp (double rhol, double ul, double pl, double cl,
                   double rhor, double ur, double pr, double cr)
    {
        constexpr double quser = 2.0;

        double const cup = 0.25 * (rhol + rhor) * (cl + cr);
        double ppv = 0.5 * (pl + pr) + 0.5 * (ul - ur) * cup;
        ppv = std::max(0.0, ppv);
        double const pmin = std::min(pl, pr);
        double const pmax = std::max(pl, pr);
        double const qmax = pmax / pmin;

        if (qmax <= quser && pmin <= ppv && ppv <= pmax) {
            return ppv;
        } else if (ppv < pmin) {
            double const pq = std::pow(pl / pr, g[1]);
            double const um = (pq * ul / cl + ur / cr + g[4] * (pq - 1.0))
                              / (pq / cl + 1.0 / cr);
            double const ptL = 1.0 + g[7] * (ul - um) / cl;
            double const ptR = 1.0 + g[7] * (um - ur) / cr;
            return 0.5 * (pl * std::pow(ptL, g[3]) + pr * std::pow(ptR, g[3]));
        } else {
            double const geL = std::sqrt((g[5] / rhol) / (g[6] * pl + ppv));
            double const geR = std::sqrt((g[5] / rhor) / (g[6] * pr + ppv));
            return (geL * pl + geR * pr - (ur - ul)) / (geL + geR);
        }
    }

    void prefun (double p, double dk, double pk, double ck,
                 double& f, double& fd)
    {
        if (p <= pk) {
            double const prat = p / pk;
            f = g[4] * ck * (std::pow(prat, g[1]) - 1.0);
            fd = std::pow(prat, -g[2]) / (dk * ck);
        } else {
            double const ak = g[5] / dk;
            double const bk = g[6] * pk;
            double const qrt = std::sqrt(ak / (bk + p));
            f = (p - pk) * qrt;
            fd = (1.0 - 0.5 * (p - pk) / (bk + p)) * qrt;
        }
    }

    void star_state (double rhol, double ul, double pl,
                     double rhor, double ur, double pr,
                     double& p, double& u)
    {
        double const cl = std::sqrt(gam * pl / rhol);
        double const cr = std::sqrt(gam * pr / rhor);
        double const uDiff = ur - ul;

        double pOld = guessp(rhol, ul, pl, cl, rhor, ur, pr, cr);
        double fL = 0.0, fR = 0.0, fLd = 0.0, fRd = 0.0;

        for (int k = 0; k < NRITER; ++k) {
            prefun(pOld, rhol, pl, cl, fL, fLd);
            prefun(pOld, rhor, pr, cr, fR, fRd);
            p = pOld - (fL + fR + uDiff) / (fLd + fRd);
            double const change = 2.0 * std::abs((p - pOld) / (p + pOld));
            if (change <= TOL) {
                break;
            }
            if (p < 0.0) {
                p = TOL;
            }
            pOld = p;
        }

        u = 0.5 * (ul + ur + fR - fL);
    }

    void sample (double S,
                 double rhol, double ul, double pl, double cl,
                 double rhor, double ur, double pr, double cr,
                 double p, double u,
                 double& rho, double& vel, double& press)
    {
        if (S <= u) {
            if (p <= pl) {
                double const shL = ul - cl;
                if (S <= shL) {
                    rho = rhol;
                    vel = ul;
                    press = pl;
                } else {
                    double const cmL = cl * std::pow(p / pl, g[1]);
                    double const stL = u - cmL;
                    if (S > stL) {
                        rho = rhol * std::pow(p / pl, 1.0 / gam);
                        vel = u;
                        press = p;
                    } else {
                        double const c = g[5] * (cl + g[7] * (ul - S));
                        rho = rhol * std::pow(c / cl, g[4]);
                        vel = g[5] * (cl + g[7] * ul + S);
                        press = pl * std::pow(c / cl, g[3]);
                    }
                }
            } else {
                double const pmL = p / pl;
                double const SL = ul - cl * std::sqrt(g[2] * pmL + g[1]);
                if (S <= SL) {
                    rho = rhol;
                    vel = ul;
                    press = pl;
                } else {
                    rho = rhol * (pmL + g[6]) / (pmL * g[6] + 1.0);
                    vel = u;
                    press = p;
                }
            }
        } else {
            if (p > pr) {
                double const pmR = p / pr;
                double const SR = ur + cr * std::sqrt(g[2] * pmR + g[1]);
                if (S >= SR) {
                    rho = rhor;
                    vel = ur;
                    press = pr;
                } else {
                    rho = rhor * (pmR + g[6]) / (pmR * g[6] + 1.0);
                    vel = u;
                    press = p;
                }
            } else {
                double const shR = ur + cr;
                if (S >= shR) {
                    rho = rhor;
                    vel = ur;
                    press = pr;
                } else {
                    double const cmR = cr * std::pow(p / pr, g[1]);
                    double const stR = u + cmR;
                    if (S <= stR) {
                        rho = rhor * std::pow(p / pr, 1.0 / gam);
                        vel = u;
                        press = p;
                    } else {
                        double const c = g[5] * (cr - g[7] * (ur - S));
                        rho = rhor * std::pow(c / cr, g[4]);
                        vel = g[5] * (-cr + g[7] * ur + S);
                        press = pr * std::pow(c / cr, g[3]);
                    }
                }
            }
        }
    }
}

void exact_riemann (std::vector<double>& x, std::vector<double>& p,
                    std::vector<double>& rho, std::vector<double>& v, double t,
                    double pl, double rhol, double vl,
                    double pr, double rhor, double vr)
{
    std::size_t const n = x.size();
    if (n == 0) {
        return;
    }

    if (t <= 0.0) {
        for (std::size_t i = 0; i < n; ++i) {
            if (x[i] < xdisc) {
                rho[i] = rhol;
                v[i] = vl;
                p[i] = pl;
            } else {
                rho[i] = rhor;
                v[i] = vr;
                p[i] = pr;
            }
        }
        return;
    }

    double const cl = std::sqrt(gam * pl / rhol);
    double const cr = std::sqrt(gam * pr / rhor);
    double pstar = 0.0;
    double ustar = 0.0;
    star_state(rhol, vl, pl, rhor, vr, pr, pstar, ustar);

    double const invt = 1.0 / t;
    for (std::size_t i = 0; i < n; ++i) {
        double const S = (x[i] - xdisc) * invt;
        sample(S, rhol, vl, pl, cl, rhor, vr, pr, cr, pstar, ustar,
               rho[i], v[i], p[i]);
    }
}
