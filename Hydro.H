#ifndef HYDRO_H_
#define HYDRO_H_

#include "FAB.H"

#include <vector>

#ifdef USE_FLOAT
using Real = float;
#else
using Real = double;
#endif

namespace {
    constexpr int URHO  = 0;
    constexpr int UMX   = 1;
    constexpr int UMY   = 2;
    constexpr int UMZ   = 3;
    constexpr int UEDEN = 4;
    constexpr int UEINT = 5;
    constexpr int NCONS = 6;

    constexpr int QRHO   = 0;
    constexpr int QU     = 1;
    constexpr int QV     = 2;
    constexpr int QW     = 3;
    constexpr int QPRES  = 4;
    constexpr int QREINT = 5;
    constexpr int NPRIM  = 6;

    constexpr int GDRHO  = 0;
    constexpr int GDU    = 1;
    constexpr int GDV    = 2;
    constexpr int GDW    = 3;
    constexpr int GDPRES = 4;
    constexpr int NGDNV  = 5;

    constexpr auto gam = Real(1.4);
    constexpr auto smallr = Real(1.e-10);
    constexpr auto smallp = Real(1.e-10);
    constexpr auto smallu = Real(1.e-10);
    constexpr auto small  = Real(1.e-10);

    constexpr int NUM_GROW = 4;
}

namespace Sod {
    constexpr auto pl   = Real(1.0);
    constexpr auto rhol = Real(1.0);
    constexpr auto ul   = Real(0.0);

    constexpr auto pr   = Real(0.1);
    constexpr auto rhor = Real(0.125);
    constexpr auto ur   = Real(0.0);
}

void init_sod (FAB<Real>& state);
Real compute_dt (FAB<Real> const& state);
void fill_boundary (FAB<Real>& state);
void cons_to_prim (FAB<Real> const& state, FAB<Real>& prim);
void godunov (FAB<Real>& state, FAB<Real> const& prim, Real dt);
Real compute_l1_error (Real t, FAB<Real> const& prim,
                       std::vector<double>& x, std::vector<double>& p,
                       std::vector<double>& rho, std::vector<double>& v);

#endif
