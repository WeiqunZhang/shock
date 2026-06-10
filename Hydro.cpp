#include "Hydro.H"
#include "Riemann.H"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
    int nx = 0;
    int ny = 0;
    int nz = 0;
    Real dx = 0;
    Real dy = 0;
    Real dz = 0;
}

void init_sod (FAB<Real>& state)
{
    auto const glo = state.lo();
    auto const ghi = state.hi();

    nx = ghi[0]-glo[0]+1-2*NUM_GROW;
    ny = ghi[1]-glo[1]+1-2*NUM_GROW;
    nz = ghi[2]-glo[2]+1-2*NUM_GROW;
    dx = Real(1) / Real(nx);
    dy = Real(1) / Real(ny);
    dz = Real(1) / Real(nz);

    int const imid = nx/2;

    for (int k = glo[2]; k <= ghi[2]; ++k) {
    for (int j = glo[1]; j <= ghi[1]; ++j) {
    for (int i = glo[0]; i <= ghi[0]; ++i) {
        Real Pt, rhot, ut;
        if (i < imid) {
            Pt = Sod::pl;
            rhot = Sod::rhol;
            ut = Sod::ul;
        } else {
            Pt = Sod::pr;
            rhot = Sod::rhor;
            ut = Sod::ur;
        }
        auto vt = Real(0.1);
        auto wt = Real(0.2);
        state(i,j,k,URHO) = rhot;
        state(i,j,k,UMX ) = rhot*ut;
        state(i,j,k,UMY ) = rhot*vt;
        state(i,j,k,UMZ ) = rhot*wt;
        Real eit = Pt/(gam-Real(1));
        Real ekt = Real(0.5) * rhot * (ut*ut + vt*vt + wt*wt);
        state(i,j,k,UEINT) = eit;
        state(i,j,k,UEDEN) = eit + ekt;
    }}}
}

Real compute_dt (FAB<Real> const& state)
{
    Real dt = std::numeric_limits<Real>::max();
    for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
        Real rho = state(i,j,k,URHO);
        Real mx  = state(i,j,k,UMX);
        Real my  = state(i,j,k,UMY);
        Real mz  = state(i,j,k,UMZ);
        Real ei  = state(i,j,k,UEINT);
        Real rhoinv = Real(1.0)/std::max(rho, smallr);
        Real vx = mx*rhoinv;
        Real vy = my*rhoinv;
        Real vz = mz*rhoinv;
        Real p = std::max((gam-Real(1.0))*ei, smallp);
        Real cs = std::sqrt(gam*p*rhoinv);
        Real dtx = dx/(std::abs(vx)+cs);
        Real dty = dy/(std::abs(vy)+cs);
        Real dtz = dz/(std::abs(vz)+cs);
        dt = std::min({dt,dtx,dty,dtz});
    }}}
    return dt;
}

// Fill ghost cells. Periodic boundaries in the y- and z-directions, and
// first-order extrapolation in the x-direction.
void fill_boundary (FAB<Real>& state)
{
    auto const glo = state.lo();
    auto const ghi = state.hi();

    for (int n = 0; n < NCONS; ++n) {
        for (int k = glo[2]; k <= ghi[2]; ++k) {
        for (int j = glo[1]; j <  0     ; ++j) {
        for (int i = glo[0]; i <= ghi[0]; ++i) {
            state(i,j,k,n) = state(i,j+ny,k,n);
        }}}

        for (int k = glo[2]; k <= ghi[2]; ++k) {
        for (int j = ny    ; j <= ghi[1]; ++j) {
        for (int i = glo[0]; i <= ghi[0]; ++i) {
            state(i,j,k,n) = state(i,j-ny,k,n);
        }}}

        for (int k = glo[2]; k <  0     ; ++k) {
        for (int j = glo[1]; j <= ghi[1]; ++j) {
        for (int i = glo[0]; i <= ghi[0]; ++i) {
            state(i,j,k,n) = state(i,j,k+nz,n);
        }}}

        for (int k = nz    ; k <= ghi[2]; ++k) {
        for (int j = glo[1]; j <= ghi[1]; ++j) {
        for (int i = glo[0]; i <= ghi[0]; ++i) {
            state(i,j,k,n) = state(i,j,k-nz,n);
        }}}

        for (int k = glo[2]; k <= ghi[2]; ++k) {
        for (int j = glo[1]; j <= ghi[1]; ++j) {
            for (int i = glo[0]; i < 0; ++i) {
                state(i,j,k,n) = state(0,j,k,n);
            }
            for (int i = nx; i <= ghi[0]; ++i) {
                state(i,j,k,n) = state(nx-1,j,k,n);
            }
        }}
    }
}

// conservative to primitive
void cons_to_prim (FAB<Real> const& state, FAB<Real>& prim)
{
    auto const glo = state.lo();
    auto const ghi = state.hi();

    for (int k = glo[2]; k <= ghi[2]; ++k) {
    for (int j = glo[1]; j <= ghi[1]; ++j) {
    for (int i = glo[0]; i <= ghi[0]; ++i) {
        Real rho = std::max(state(i,j,k,URHO), smallr);
        Real rhoinv = Real(1.0)/rho;
        Real ux = state(i,j,k,UMX)*rhoinv;
        Real uy = state(i,j,k,UMY)*rhoinv;
        Real uz = state(i,j,k,UMZ)*rhoinv;
        Real kineng = Real(0.5)*rho*(ux*ux + uy*uy + uz*uz);
        Real ei = state(i,j,k,UEDEN) - kineng;
        if (ei <= Real(0.0)) { ei = state(i,j,k,UEINT); }
        ei = std::max(ei, smallp/(gam-Real(1.0)));
        Real p = std::max((gam-Real(1.0))*ei, smallp);

        prim(i,j,k,QRHO  ) = rho;
        prim(i,j,k,QU    ) = ux;
        prim(i,j,k,QV    ) = uy;
        prim(i,j,k,QW    ) = uz;
        prim(i,j,k,QPRES ) = p;
        prim(i,j,k,QREINT) = ei;
    }}}
}

namespace {
// Compute the coefficients of a parabolic reconstruction of the data in a
// zone. This uses the standard PPM limiters described in Colella & Woodward
// (1984)
FORCE_INLINE void ppm_reconstruct (const Real s[5], Real flatn, Real& sm, Real& sp)
{
  // Compute van Leer slopes

    Real dsl = Real(2.0) * (s[1] - s[0]);
    Real dsr = Real(2.0) * (s[2] - s[1]);

    Real dsvl_l = Real(0.0);
    if (dsl * dsr > Real(0.0)) {
        Real dsc = Real(0.5) * (s[2] - s[0]);
        dsvl_l = std::copysign(Real(1.0), dsc) *
            std::min({std::abs(dsc), std::abs(dsl), std::abs(dsr)});
    }

    dsl = Real(2.0) * (s[2] - s[1]);
    dsr = Real(2.0) * (s[3] - s[2]);

    Real dsvl_r = Real(0.0);
    if (dsl * dsr > Real(0.0)) {
        Real dsc = Real(0.5) * (s[3] - s[1]);
        dsvl_r = std::copysign(Real(1.0), dsc) *
            std::min({std::abs(dsc), std::abs(dsl), std::abs(dsr)});
    }

    // Interpolate s to edges
    sm = Real(0.5) * (s[2] + s[1]) - Real(1.0 / 6.0) * (dsvl_r - dsvl_l);

    // Make sure sedge lies in between adjacent cell-centered values
    sm = std::max(sm, std::min(s[2], s[1]));
    sm = std::min(sm, std::max(s[2], s[1]));

    // Compute van Leer slopes
    dsl = Real(2.0) * (s[2] - s[1]);
    dsr = Real(2.0) * (s[3] - s[2]);

    dsvl_l = Real(0.0);
    if (dsl * dsr > Real(0.0)) {
        Real dsc = Real(0.5) * (s[3] - s[1]);
        dsvl_l = std::copysign(Real(1.0), dsc) *
            std::min({std::abs(dsc), std::abs(dsl), std::abs(dsr)});
    }

    dsl = Real(2.0) * (s[3] - s[2]);
    dsr = Real(2.0) * (s[4] - s[3]);

    dsvl_r = Real(0.0);
    if (dsl * dsr > Real(0.0)) {
        Real dsc = Real(0.5) * (s[4] - s[2]);
        dsvl_r = std::copysign(Real(1.0), dsc) *
            std::min({std::abs(dsc), std::abs(dsl), std::abs(dsr)});
    }

    // Interpolate s to edges
    sp = Real(0.5) * (s[3] + s[2]) - Real(1.0 / 6.0) * (dsvl_r - dsvl_l);

    // Make sure sedge lies in between adjacent cell-centered values
    sp = std::max(sp, std::min(s[3], s[2]));
    sp = std::min(sp, std::max(s[3], s[2]));

    // Flatten the parabola
    sm = flatn * sm + (Real(1.0) - flatn) * s[2];
    sp = flatn * sp + (Real(1.0) - flatn) * s[2];

    // Modify using quadratic limiters -- note this version of the limiting comes
    // from Colella and Sekora (2008), not the original PPM paper.
    if ((sp - s[2]) * (s[2] - sm) <= Real(0.0)) {
        sp = s[2];
        sm = s[2];
    } else if (
        std::abs(sp - s[2]) >= Real(2.0) * std::abs(sm - s[2])) {
        sp = Real(3.0) * s[2] - Real(2.0) * sm;
    } else if (
        std::abs(sm - s[2]) >= Real(2.0) * std::abs(sp - s[2])) {
        sm = Real(3.0) * s[2] - Real(2.0) * sp;
    }
}

// Integrate under the parabola using from the left and right edges
// with the wave speeds u-c, u, u+c
FORCE_INLINE void ppm_int_profile (Real sm, Real sp, Real sc, Real u, Real c, Real dtdx,
                                   Real* RESTRICT Ip, Real* RESTRICT Im)
{
    // Integrate the parabolic profile to the edge of the cell.

    // compute x-component of Ip and Im
    Real s6 = Real(6.0) * sc - Real(3.0) * (sm + sp);

    // Ip/m is the integral under the parabola for the extent
    // that a wave can travel over a timestep

    // Ip integrates to the right edge of a cell
    // Im integrates to the left edge of a cell

    // u-c wave
    Real speed = u - c;
    Real sigma = std::abs(speed) * dtdx;

    // if speed == ZERO, then either branch is the same
    if (speed <= Real(0.0)) {
        Ip[0] = sp;
        Im[0] = sm + Real(0.5) * sigma * (sp - sm + (Real(1.0) - Real(2.0 / 3.0) * sigma) * s6);
    } else {
        Ip[0] = sp - Real(0.5) * sigma * (sp - sm - (Real(1.0) - Real(2.0 / 3.0) * sigma) * s6);
        Im[0] = sm;
    }

    // u wave
    speed = u;
    sigma = std::abs(speed) * dtdx;

    if (speed <= Real(0.0)) {
        Ip[1] = sp;
        Im[1] = sm + Real(0.5) * sigma * (sp - sm + (Real(1.0) - Real(2.0 / 3.0) * sigma) * s6);
    } else {
        Ip[1] = sp - Real(0.5) * sigma * (sp - sm - (Real(1.0) - Real(2.0 / 3.0) * sigma) * s6);
        Im[1] = sm;
    }

    // u+c wave
    speed = u + c;
    sigma = std::abs(speed) * dtdx;

    if (speed <= Real(0.0)) {
        Ip[2] = sp;
        Im[2] = sm + Real(0.5) * sigma * (sp - sm + (Real(1.0) - Real(2.0 / 3.0) * sigma) * s6);
    } else {
        Ip[2] = sp - Real(0.5) * sigma * (sp - sm - (Real(1.0) - Real(2.0 / 3.0) * sigma) * s6);
        Im[2] = sm;
    }
}

void trace_ppm (int dir, FAB<Real> const& q, FAB<Real>& qm, FAB<Real>& qp, Real dtdx)
{
    int QUN = 0;
    int QUT = 0;
    int QUTT = 0;
    if (dir == 0) {
        QUN = QU;
        QUT = QV;
        QUTT = QW;
    } else if (dir == 1) {
        QUN = QV;
        QUT = QW;
        QUTT = QU;
    } else if (dir == 2) {
        QUN = QW;
        QUT = QU;
        QUTT = QV;
    }

    // Trace to left and right edges using upwind PPM
    for (int k = -2; k < nz+2; ++k) {
    for (int j = -2; j < ny+2; ++j) {
    for (int i = -2; i < nx+2; ++i) {
        int im1 = (dir == 0) ? i-1 : i;
        int im2 = (dir == 0) ? i-2 : i;
        int jm1 = (dir == 1) ? j-1 : j;
        int jm2 = (dir == 1) ? j-2 : j;
        int km1 = (dir == 2) ? k-1 : k;
        int km2 = (dir == 2) ? k-2 : k;
        int ip1 = (dir == 0) ? i+1 : i;
        int ip2 = (dir == 0) ? i+2 : i;
        int jp1 = (dir == 1) ? j+1 : j;
        int jp2 = (dir == 1) ? j+2 : j;
        int kp1 = (dir == 2) ? k+1 : k;
        int kp2 = (dir == 2) ? k+2 : k;

        Real cc = std::sqrt(gam * q(i,j,k,QPRES)/q(i,j,k,QRHO));
        Real un = q(i,j,k,QUN);

        // do the parabolic reconstruction and compute the
        // integrals under the characteristic waves

        Real flat = Real(1.0); // TODO: add flattening

        Real Ip[NPRIM][3];
        Real Im[NPRIM][3];

        for (int n = 0; n < NPRIM; n++)
        {
            Real s[5];
            s[0] = q(im2,jm2,km2,n);
            s[1] = q(im1,jm1,km1,n);
            s[2] = q(i  ,j  , k, n);
            s[3] = q(ip1,jp1,kp1,n);
            s[4] = q(ip2,jp2,kp2,n);
            Real sm, sp;
            ppm_reconstruct(s, flat, sm, sp);
            ppm_int_profile(sm, sp, s[2], un, cc, dtdx, Ip[n], Im[n]);
        }

        // plus state on face i
        if ((dir == 0 && i >= -1) ||
            (dir == 1 && j >= -1) ||
            (dir == 2 && k >= -1)) {
            // Set the reference state
            // This will be the fastest moving state to the left --
            // this is the method that Miller & Colella and Colella &
            // Woodward use
            Real rho_ref = Im[QRHO][0];
            Real un_ref = Im[QUN][0];

            Real p_ref = Im[QPRES][0];
            Real rhoe_g_ref = Im[QREINT][0];

            rho_ref = std::max( rho_ref, smallr);
            Real rho_ref_inv = Real(1.0) / rho_ref;
            p_ref = std::max(p_ref, smallp);

            // For tracing
            Real cc_ref = std::sqrt(gam*p_ref/rho_ref);
            Real csq_ref = cc_ref * cc_ref;
            Real cc_ref_inv = Real(1.0) / cc_ref;
            Real h_g_ref = (p_ref + rhoe_g_ref) * rho_ref_inv / csq_ref;

            // *m are the jumps carried by un-c
            // *p are the jumps carried by un+c

            // Note: for the transverse velocities, the jump is carried
            //       only by the u wave (the contact)

            // we also add the sources here so they participate in the tracing
            Real dum = un_ref - Im[QUN][0];
            Real dptotm = p_ref - Im[QPRES][0];

            Real drho = rho_ref - Im[QRHO][1];
            Real dptot = p_ref - Im[QPRES][1];
            Real drhoe_g = rhoe_g_ref - Im[QREINT][1];

            Real dup = un_ref - Im[QUN][2];
            Real dptotp = p_ref - Im[QPRES][2];

            // {rho, u, p, (rho e)} eigensystem

            // These are analogous to the beta's from the original PPM
            // paper (except we work with rho instead of tau).  This is
            // simply (l . dq), where dq = qref - I(q)

            Real alpham = Real(0.5) * (dptotm * rho_ref_inv * cc_ref_inv - dum)
                * rho_ref * cc_ref_inv;
            Real alphap = Real(0.5) * (dptotp * rho_ref_inv * cc_ref_inv + dup)
                * rho_ref * cc_ref_inv;
            Real alpha0r = drho - dptot / csq_ref;
            Real alpha0e_g = drhoe_g - dptot * h_g_ref;

            if (un-cc > 0) {
                alpham = Real(0);
            } else if (un-cc < Real(0)) {
                alpham *= Real(-1.0);
            } else {
                alpham *= Real(-0.5);
            }

            if (un+cc > 0) {
                alphap = Real(0);
            } else if (un+cc < Real(0)) {
                alphap *= Real(-1.0);
            } else {
                alphap *= Real(-0.5);
            }

            if (un > 0) {
                alpha0r   = Real(0);
                alpha0e_g = Real(0);
            } else if (un < Real(0)) {
                alpha0r   *= Real(-1.0);
                alpha0e_g *= Real(-1.0);
            } else {
                alpha0r   *= Real(-0.5);
                alpha0e_g *= Real(-0.5);
            }

            // The final interface states are just
            // q_s = q_ref - sum(l . dq) r
            // note that the a{mpz}right as defined above have the minus already
            qp(i,j,k,QRHO) = rho_ref + alphap + alpham + alpha0r;
            qp(i,j,k,QUN) = un_ref + (alphap - alpham) * cc_ref * rho_ref_inv;
            qp(i,j,k,QPRES) = p_ref + (alphap + alpham) * csq_ref;

            qp(i,j,k,QRHO)  = std::max( qp(i,j,k,QRHO ), smallr);
            qp(i,j,k,QPRES) = std::max( qp(i,j,k,QPRES), smallp);
            // Transverse velocities -- there's no projection here, so we
            // don't need a reference state.  We only care about the state
            // traced under the middle wave

            // Recall that I already takes the limit of the parabola
            // in the event that the wave is not moving toward the
            // interface
            qp(i,j,k,QUT) = Im[QUT][1];
            qp(i,j,k,QUTT) = Im[QUTT][1];

            // This allows the (rho e) to take advantage of (pressure > smallp)
            qp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref*csq_ref + alpha0e_g;
        }

        // minus state on face i + 1
        if ( (dir == 0 && i <= nx) ||
             (dir == 1 && j <= ny) ||
             (dir == 2 && k <= nz)) {
            // Set the reference state
            // This will be the fastest moving state to the right
            Real rho_ref = Ip[QRHO][2];
            Real un_ref = Ip[QUN][2];

            Real p_ref = Ip[QPRES][2];
            Real rhoe_g_ref = Ip[QREINT][2];

            rho_ref = std::max( rho_ref, smallr );
            Real rho_ref_inv = Real(1.0) / rho_ref;
            p_ref = std::max( p_ref, smallp );

            // For tracing
            Real cc_ref = std::sqrt(gam*p_ref/rho_ref);

            Real csq_ref = cc_ref * cc_ref;
            Real cc_ref_inv = Real(1.0) / cc_ref;
            Real h_g_ref = (p_ref + rhoe_g_ref) * rho_ref_inv / csq_ref;

            // *m are the jumps carried by u-c
            // *p are the jumps carried by u+c

            Real dum = un_ref - Ip[QUN][0];
            Real dptotm = p_ref - Ip[QPRES][0];

            Real drho = rho_ref - Ip[QRHO][1];
            Real dptot = p_ref - Ip[QPRES][1];
            Real drhoe_g = rhoe_g_ref - Ip[QREINT][1];

            Real dup = un_ref - Ip[QUN][2];
            Real dptotp = p_ref - Ip[QPRES][2];

            // {rho, u, p, (rho e)} eigensystem

            // These are analogous to the beta's from the original PPM
            // paper (except we work with rho instead of tau).  This is
            // simply (l . dq), where dq = qref - I(q)

            Real alpham = Real(0.5) * (dptotm * rho_ref_inv * cc_ref_inv - dum)
                * rho_ref * cc_ref_inv;
            Real alphap = Real(0.5) * (dptotp * rho_ref_inv * cc_ref_inv + dup)
                * rho_ref * cc_ref_inv;
            Real alpha0r = drho - dptot / csq_ref;
            Real alpha0e_g = drhoe_g - dptot * h_g_ref;

            if (un-cc > 0) {
                alpham *= Real(-1.0);
            } else if (un-cc < Real(0)) {
                alpham = Real(0);
            } else {
                alpham *= Real(-0.5);
            }

            if (un+cc > 0) {
                alphap *= Real(-1.0);
            } else if (un+cc < Real(0)) {
                alphap = 0.;
            } else {
                alphap *= Real(-0.5);
            }

            if (un > 0) {
                alpha0r   *= Real(-1.0);
                alpha0e_g *= Real(-1.0);
            } else if (un < Real(0)) {
                alpha0r   = Real(0);
                alpha0e_g = Real(0);
            } else {
                alpha0r   *= Real(-0.5);
                alpha0e_g *= Real(-0.5);
            }

            // The final interface states are just
            // q_s = q_ref - sum (l . dq) r
            // note that the a{mpz}left as defined above have the minus already
            qm(ip1,jp1,kp1,QRHO) = rho_ref + alphap + alpham + alpha0r;
            qm(ip1,jp1,kp1,QUN) = un_ref + (alphap - alpham) * cc_ref * rho_ref_inv;
            qm(ip1,jp1,kp1,QPRES) = p_ref + (alphap + alpham) * csq_ref;

            qm(ip1,jp1,kp1,QRHO)  = std::max( qm(ip1,jp1,kp1,QRHO), smallr);
            qm(ip1,jp1,kp1,QPRES) = std::max( qm(ip1,jp1,kp1,QPRES), smallp);

            // transverse velocities
            qm(ip1,jp1,kp1,QUT) = Ip[QUT][1];
            qm(ip1,jp1,kp1,QUTT) = Ip[QUTT][1];

            // This allows the (rho e) to take advantage of (pressure >
            // smallp)
            qm(ip1,jp1,kp1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref*csq_ref + alpha0e_g;
        }
  }}}
}

FORCE_INLINE
void riemannus (Real rl_in, Real ul, Real vl, Real v2l, Real pl_in, Real rel,
                Real rr_in, Real ur, Real vr, Real v2r, Real pr_in, Real rer,
                Real& uflx_rho,
                Real& uflx_u, Real& uflx_v, Real& uflx_w,
                Real& uflx_eden, Real& uflx_eint,
                Real& qint_iu, Real& qint_iv1, Real& qint_iv2, Real& qint_gdpres)
{
    Real wsmall = smallu * smallr;
    Real rl = std::max(rl_in,smallr);
    Real pl = std::max(pl_in,smallp);
    Real rr = std::max(rr_in,smallr);
    Real pr = std::max(pr_in,smallp);

    Real wl = std::max(wsmall, std::sqrt(std::abs(gam * pl * rl)));
    Real wr = std::max(wsmall, std::sqrt(std::abs(gam * pr * rr)));

    Real pstar = ((wr * pl + wl * pr) + wl * wr * (ul - ur)) / (wl + wr);
    pstar = std::max(pstar,smallp);
    Real ustar = ((wl * ul + wr * ur) + (pl - pr)) / (wl + wr);

    Real ro,uo,po,reo;

    bool mask = std::abs(ustar) < smallu * Real(0.5) * (std::abs(ul) + std::abs(ur));
    ustar = mask ? Real(0.0) : ustar;

    if (ustar > Real(0.)) {
        ro = rl;
        uo = ul;
        po = pl;
        reo = rel;
    } else if (ustar < Real(0.)) {
        ro = rr;
        uo = ur;
        po = pr;
        reo = rer;
    } else {
        ro = Real(0.5)*(rl+rr);
        uo = Real(0.5)*(ul+ur);
        po = Real(0.5)*(pl+pr);
        reo = Real(0.5)*(rel+rer);
    }

    ro = std::max(ro,smallr);
    Real roinv = Real(1.0)/ro;

    Real co,co2inv,rstar,cstar,entho,estar,sgnm,spout,spin,ushock,frac,drho;

    co = std::sqrt(std::abs(gam*po*roinv));
    co = std::max(smallu,co);
    co2inv = Real(1.0)/(co*co);

    drho = (pstar - po)*co2inv;
    rstar = ro + drho;
    rstar = std::max(smallr,rstar);

    entho = (reo + po)*roinv*co2inv;
    estar = reo + (pstar - po)*entho;
    cstar = std::sqrt(std::abs(gam*pstar/rstar));
    cstar = std::max(cstar,smallu);

    sgnm = (ustar < Real(0.)) ? Real(-1.0) : Real(1.0);
    spout = co - sgnm*uo;
    spin = cstar - sgnm*ustar;
    ushock = Real(0.5)*(spin + spout);

    if (pstar-po > Real(0.0)) {
        spin = ushock;
        spout = ushock;
    }

    Real scr;
    if (spout-spin == Real(0.0)){
        Real cav = Real(0.5) * (std::sqrt(gam*pl/rl) + std::sqrt(gam*pr/rr));
        scr = smallu*cav;
    } else {
        scr = spout-spin;
    }

    frac = (Real(1.0) + (spout + spin)/scr)*Real(0.5);
    frac = std::max(Real(0.0),std::min(Real(1.0),frac));

    mask = ustar > Real(0.0);
    qint_iv1 = mask ? vl : vr;
    qint_iv2 = mask ? v2l : v2r;

    mask = (ustar == Real(0.0));
    qint_iv1 = mask ? Real(0.5) * (vl + vr) : qint_iv1;
    qint_iv2 = mask ? Real(0.5) * (v2l + v2r) : qint_iv2;

    Real rho_gd = frac*rstar+(Real(1.0)-frac)*ro;
    qint_iu = frac * ustar + (Real(1.0) - frac) * uo;
    qint_gdpres = frac * pstar + (Real(1.0) - frac) * po;
    Real re_gd = frac*estar + (Real(1.0)-frac)*reo;

    mask = (spout < Real(0.0));
    rho_gd = mask ? ro : rho_gd;
    qint_iu = mask ? uo : qint_iu;
    qint_gdpres = mask ? po : qint_gdpres;
    re_gd = mask ? reo : re_gd;

    mask = (spin >= Real(0.0));
    rho_gd = mask ? rstar : rho_gd;
    qint_iu = mask ? ustar : qint_iu;
    qint_gdpres = mask ? pstar : qint_gdpres;
    re_gd = mask ? estar : re_gd;

    qint_gdpres = std::max(qint_gdpres,smallp);

    uflx_rho = rho_gd * qint_iu;

    uflx_u = uflx_rho * qint_iu + qint_gdpres;
    uflx_v = uflx_rho * qint_iv1;
    uflx_w = uflx_rho * qint_iv2;
    Real rhoetot = re_gd +
        Real(0.5) * rho_gd * (qint_iu * qint_iu + qint_iv1 * qint_iv1 + qint_iv2 * qint_iv2);
    uflx_eden = qint_iu * (rhoetot + qint_gdpres);
    uflx_eint = qint_iu * re_gd;
}

void cmpflx (int dir, FAB<Real>& flx, FAB<Real> const& ql, FAB<Real> const& qr,
             FAB<Real> & q)
{
    int IU, IV, IV2, GU, GV, GV2;
    int f_idx[3];
    if (dir == 0) {
        IU = QU;
        IV = QV;
        IV2 = QW;
        GU = GDU;
        GV = GDV;
        GV2 = GDW;
        f_idx[0] = UMX;
        f_idx[1] = UMY;
        f_idx[2] = UMZ;
    } else if (dir == 1) {
        IU = QV;
        IV = QU;
        IV2 = QW;
        GU = GDV;
        GV = GDU;
        GV2 = GDW;
        f_idx[0] = UMY;
        f_idx[1] = UMX;
        f_idx[2] = UMZ;
    } else {
        IU = QW;
        IV = QU;
        IV2 = QV;
        GU = GDW;
        GV = GDU;
        GV2 = GDV;
        f_idx[0] = UMZ;
        f_idx[1] = UMX;
        f_idx[2] = UMY;
    }

    auto const flo = flx.lo();
    auto const fhi = flx.hi();

    for (int k = flo[2]; k <= fhi[2]; ++k) {
    for (int j = flo[1]; j <= fhi[1]; ++j) {
    for (int i = flo[0]; i <= fhi[0]; ++i) {
        Real ul = ql(i,j,k,IU);
        Real vl = ql(i,j,k,IV);
        Real v2l = ql(i,j,k,IV2);

        Real ur = qr(i,j,k,IU);
        Real vr = qr(i,j,k,IV);
        Real v2r = qr(i,j,k,IV2);

        riemannus(ql(i,j,k,QRHO), ul, vl, v2l, ql(i,j,k,QPRES), ql(i,j,k,QREINT),
                  qr(i,j,k,QRHO), ur, vr, v2r, qr(i,j,k,QPRES), qr(i,j,k,QREINT),
                  flx(i,j,k,URHO),
                  flx(i,j,k,f_idx[0]), flx(i,j,k,f_idx[1]), flx(i,j,k,f_idx[2]),
                  flx(i,j,k,UEDEN), flx(i,j,k,UEINT),
                  q(i,j,k,GU), q(i,j,k,GV), q(i,j,k,GV2), q(i,j,k,GDPRES));
    }}}
}

// Transverse Correction for Predicted dir-states, using other_dir-Flux
void transdo (int dir, int other_dir, FAB<Real>& qm, FAB<Real>& qp,
              FAB<Real> const& qnormm, FAB<Real> const& qnormp,
              FAB<Real> const& flxx, FAB<Real> const& qint, Real cdtdx)
{
    int ilo = (dir == 0) ? -2 : -1;
    int jlo = (dir == 1) ? -2 : -1;
    int klo = (dir == 2) ? -2 : -1;
    int ihi = (dir == 0) ? nx+1 : nx;
    int jhi = (dir == 1) ? ny+1 : ny;
    int khi = (dir == 2) ? nz+1 : nz;

    for (int k = klo; k <= khi; ++k) {
    for (int j = jlo; j <= jhi; ++j) {
    for (int i = ilo; i <= ihi; ++i) {
        int in = (      dir == 0) ? i+1 : i;
        int jn = (      dir == 1) ? j+1 : j;
        int kn = (      dir == 2) ? k+1 : k;
        int it = (other_dir == 0) ? i+1 : i;
        int jt = (other_dir == 1) ? j+1 : j;
        int kt = (other_dir == 2) ? k+1 : k;
        int qvidx = (other_dir == 0) ? GDU : ((other_dir == 1) ? GDV : GDW);

        Real flxrho = cdtdx * (flxx(it,jt,kt,URHO ) - flxx(i,j,k,URHO ));
        Real flxu   = cdtdx * (flxx(it,jt,kt,UMX  ) - flxx(i,j,k,UMX  ));
        Real flxv   = cdtdx * (flxx(it,jt,kt,UMY  ) - flxx(i,j,k,UMY  ));
        Real flxw   = cdtdx * (flxx(it,jt,kt,UMZ  ) - flxx(i,j,k,UMZ  ));
        Real flxe   = cdtdx * (flxx(it,jt,kt,UEDEN) - flxx(i,j,k,UEDEN));

        // Update hydro vars

        Real pggp = qint(it,jt,kt,GDPRES);
        Real pggm = qint(i,j,k,GDPRES);
        Real ugp = qint(it,jt,kt,qvidx);
        Real ugm = qint(i,j,k,qvidx);
        Real dAup = pggp * ugp - pggm * ugm;
        Real pav = Real(0.5) * (pggp + pggm);
        Real dAu = ugp - ugm;

        // QP
        if (((dir == 0) && (i > -2)) ||
            ((dir == 1) && (j > -2)) ||
            ((dir == 2) && (k > -2)))
        {
            // Convert to conservative
            Real rrr = qnormp(i,j,k,QRHO);
            Real rur = qnormp(i,j,k,QU);
            Real rvr = qnormp(i,j,k,QV);
            Real rwr = qnormp(i,j,k,QW);
            Real ekinr = Real(0.5) * rrr * (rur * rur + rvr * rvr + rwr * rwr);
            rur *= rrr;
            rvr *= rrr;
            rwr *= rrr;

            Real rer = qnormp(i,j,k,QREINT) + ekinr;
            // Add transverse predictor
            Real rrnewr = rrr - flxrho;
            Real runewr = rur - flxu;
            Real rvnewr = rvr - flxv;
            Real rwnewr = rwr - flxw;
            Real renewr = rer - flxe;

            // Convert back to primitive
            qp(i,j,k,QRHO) = rrnewr;
            qp(i,j,k,QU) = runewr / rrnewr;
            qp(i,j,k,QV) = rvnewr / rrnewr;
            qp(i,j,k,QW) = rwnewr / rrnewr;

            Real rhoekinr = Real(0.5) * (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;

            qp(i,j,k,QREINT) = renewr - rhoekinr;

            Real pnewr = qnormp(i,j,k,QPRES) - cdtdx * (dAup + pav * dAu * (gam - Real(1.0)));
            qp(i,j,k,QPRES) = std::max(pnewr, smallp);
        }

        // QM
        if (((dir == 0) && (i <= nx)) ||
            ((dir == 1) && (j <= ny)) ||
            ((dir == 2) && (k <= nz)))
        {
            // Conversion to Conservative
            Real rrl = qnormm(in,jn,kn,QRHO);
            Real rul = qnormm(in,jn,kn,QU);
            Real rvl = qnormm(in,jn,kn,QV);
            Real rwl = qnormm(in,jn,kn,QW);
            Real ekinl = Real(0.5) * rrl * (rul * rul + rvl * rvl + rwl * rwl);
            rul *= rrl;
            rvl *= rrl;
            rwl *= rrl;
            Real rel = qnormm(in,jn,kn,QREINT) + ekinl;

            // Transverse fluxes
            Real rrnewl = rrl - flxrho;
            Real runewl = rul - flxu;
            Real rvnewl = rvl - flxv;
            Real rwnewl = rwl - flxw;
            Real renewl = rel - flxe;

            qm(in,jn,kn,QRHO) = rrnewl;
            qm(in,jn,kn,QU) = runewl / rrnewl;
            qm(in,jn,kn,QV) = rvnewl / rrnewl;
            qm(in,jn,kn,QW) = rwnewl / rrnewl;

            Real rhoekinl = Real(0.5) * (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;

            qm(in,jn,kn,QREINT) = renewl - rhoekinl;

            Real pnewl = qnormm(in,jn,kn,QPRES) - cdtdx * (dAup + pav * dAu * (gam - Real(1.0)));
            qm(in,jn,kn,QPRES) = std::max(pnewl, smallp);
        }
    }}}
}

// dir corrected from other two dirs
void transdd (int dir, FAB<Real>& qm, FAB<Real>& qp,
              FAB<Real> const& qnormm, FAB<Real> const& qnormp,
              FAB<Real> const& flxx, FAB<Real> const& flxy,
              FAB<Real> const& qx, FAB<Real> const& qy,
              Real cdtdx0, Real cdtdx1)
{
    int ilo = (dir == 0) ? -1 : 0;
    int jlo = (dir == 1) ? -1 : 0;
    int klo = (dir == 2) ? -1 : 0;
    int ihi = nx;
    int jhi = ny;
    int khi = nz;

    for (int k = klo; k <= khi; ++k) {
    for (int j = jlo; j <= jhi; ++j) {
    for (int i = ilo; i <= ihi; ++i) {
        int in, jn, kn, it0, jt0, kt0, it1, jt1, kt1, qvidx0, qvidx1;
        if (dir == 0) {
            in = i+1;
            jn = j;
            kn = k;
            it0 = i;
            jt0 = j+1;
            kt0 = k;
            it1 = i;
            jt1 = j;
            kt1 = k+1;
            qvidx0 = GDV;
            qvidx1 = GDW;
        } else if (dir == 1) {
            in = i;
            jn = j+1;
            kn = k;
            it0 = i+1;
            jt0 = j;
            kt0 = k;
            it1 = i;
            jt1 = j;
            kt1 = k+1;
            qvidx0 = GDU;
            qvidx1 = GDW;
        } else {
            in = i;
            jn = j;
            kn = k+1;
            it0 = i+1;
            jt0 = j;
            kt0 = k;
            it1 = i;
            jt1 = j+1;
            kt1 = k;
            qvidx0 = GDU;
            qvidx1 = GDV;
        }

        Real flxrho = cdtdx0 * (flxx(it0,jt0,kt0,URHO ) - flxx(i,j,k,URHO ))
            +         cdtdx1 * (flxy(it1,jt1,kt1,URHO ) - flxy(i,j,k,URHO ));
        Real flxu   = cdtdx0 * (flxx(it0,jt0,kt0,UMX  ) - flxx(i,j,k,UMX  ))
            +         cdtdx1 * (flxy(it1,jt1,kt1,UMX  ) - flxy(i,j,k,UMX  ));
        Real flxv   = cdtdx0 * (flxx(it0,jt0,kt0,UMY  ) - flxx(i,j,k,UMY  ))
            +         cdtdx1 * (flxy(it1,jt1,kt1,UMY  ) - flxy(i,j,k,UMY  ));
        Real flxw   = cdtdx0 * (flxx(it0,jt0,kt0,UMZ  ) - flxx(i,j,k,UMZ  ))
            +         cdtdx1 * (flxy(it1,jt1,kt1,UMZ  ) - flxy(i,j,k,UMZ  ));
        Real flxe   = cdtdx0 * (flxx(it0,jt0,kt0,UEDEN) - flxx(i,j,k,UEDEN))
            +         cdtdx1 * (flxy(it1,jt1,kt1,UEDEN) - flxy(i,j,k,UEDEN));

        // Update hydro vars

        Real pggpx = qx(it0,jt0,kt0,GDPRES);
        Real pggmx = qx(i  ,j  ,k  ,GDPRES);
        Real  ugpx = qx(it0,jt0,kt0,qvidx0);
        Real  ugmx = qx(i  ,j  ,k  ,qvidx0);

        Real dAupx = pggpx * ugpx - pggmx * ugmx;
        Real pavx = Real(0.5) * (pggpx + pggmx);
        Real dAux = ugpx - ugmx;

        Real pggpy = qy(it1,jt1,kt1,GDPRES);
        Real pggmy = qy(i  ,j  ,k  ,GDPRES);
        Real ugpy  = qy(it1,jt1,kt1,qvidx1);
        Real ugmy  = qy(i  ,j  ,k  ,qvidx1);

        Real dAupy = pggpy * ugpy - pggmy * ugmy;
        Real pavy = Real(0.5) * (pggpy + pggmy);
        Real dAuy = ugpy - ugmy;

        Real pxnew = cdtdx0 * (dAupx + pavx * dAux * (gam - Real(1.0)));
        Real pynew = cdtdx1 * (dAupy + pavy * dAuy * (gam - Real(1.0)));

        // qp state
        if (((dir == 0) && (i >= 0)) ||
            ((dir == 1) && (j >= 0)) ||
            ((dir == 2) && (k >= 0)))
        {
            Real rrr = qnormp(i,j,k,QRHO);
            Real rrnewr = rrr - flxrho;
            Real rur = rrr * qnormp(i,j,k,QU);
            Real rvr = rrr * qnormp(i,j,k,QV);
            Real rwr = rrr * qnormp(i,j,k,QW);
            Real ekinr = Real(0.5) * (rur * rur + rvr * rvr + rwr * rwr) / rrr;
            Real rer = qnormp(i,j,k,QREINT) + ekinr;

            Real runewr = rur - flxu;
            Real rvnewr = rvr - flxv;
            Real rwnewr = rwr - flxw;
            Real renewr = rer - flxe;

            qp(i,j,k,QRHO) = rrnewr;
            qp(i,j,k,QU) = runewr / rrnewr;
            qp(i,j,k,QV) = rvnewr / rrnewr;
            qp(i,j,k,QW) = rwnewr / rrnewr;

            Real rhoekinr = Real(0.5) *
                (runewr * runewr + rvnewr * rvnewr + rwnewr * rwnewr) / rrnewr;
            qp(i,j,k,QREINT) = renewr - rhoekinr;

            qp(i,j,k,QPRES) =  qnormp(i,j,k,QPRES) - pxnew - pynew;

            qp(i,j,k,QPRES) = std::max( qp(i,j,k,QPRES), smallp);
        }

        // qm state
        if (((dir == 0) && (i < nx)) ||
            ((dir == 1) && (j < ny)) ||
            ((dir == 2) && (k < nz)))
        {
            Real rrl = qnormm(in,jn,kn,QRHO);
            Real rrnewl = rrl - flxrho;
            Real rul = rrl * qnormm(in,jn,kn,QU);
            Real rvl = rrl * qnormm(in,jn,kn,QV);
            Real rwl = rrl * qnormm(in,jn,kn,QW);
            Real ekinl = Real(0.5) * (rul * rul + rvl * rvl + rwl * rwl) / rrl;
            Real rel = qnormm(in,jn,kn,QREINT) + ekinl;

            Real runewl = rul - flxu;
            Real rvnewl = rvl - flxv;
            Real rwnewl = rwl - flxw;
            Real renewl = rel - flxe;

            qm(in,jn,kn,QRHO) = rrnewl;
            qm(in,jn,kn,QU) = runewl / rrnewl;
            qm(in,jn,kn,QV) = rvnewl / rrnewl;
            qm(in,jn,kn,QW) = rwnewl / rrnewl;

            Real rhoekinl = Real(0.5) *
                (runewl * runewl + rvnewl * rvnewl + rwnewl * rwnewl) / rrnewl;
            qm(in,jn,kn,QREINT) = renewl - rhoekinl;

            qm(in,jn,kn,QPRES) = qnormm(in,jn,kn,QPRES) - pxnew - pynew;

            qm(in,jn,kn,QPRES) = std::max(qm(in,jn,kn,QPRES), smallp);
        }
    }}}
}

}

// Godunov method using unsplit PPM
void godunov (FAB<Real>& state, FAB<Real> const& prim, Real dt)
{
    Real const dtdx = dt/dx;
    Real const dtdy = dt/dy;
    Real const dtdz = dt/dz;

    Real const hdtdx = Real(0.5)*dtdx;
    Real const hdtdy = Real(0.5)*dtdy;
    Real const hdtdz = Real(0.5)*dtdz;

    Real const cdtdx = Real(1./3.)*dtdx;
    Real const cdtdy = Real(1./3.)*dtdy;
    Real const cdtdz = Real(1./3.)*dtdz;

    FAB<Real> qxm({-1,-2,-2}, {nx+1,ny+1,nz+1}, NPRIM);
    FAB<Real> qxp({-1,-2,-2}, {nx+1,ny+1,nz+1}, NPRIM);

    FAB<Real> qym({-2,-1,-2}, {nx+1,ny+1,nz+1}, NPRIM);
    FAB<Real> qyp({-2,-1,-2}, {nx+1,ny+1,nz+1}, NPRIM);

    FAB<Real> qzm({-2,-2,-1}, {nx+1,ny+1,nz+1}, NPRIM);
    FAB<Real> qzp({-2,-2,-1}, {nx+1,ny+1,nz+1}, NPRIM);

    trace_ppm(0, prim, qxm, qxp, dtdx);
    trace_ppm(1, prim, qym, qyp, dtdy);
    trace_ppm(2, prim, qzm, qzp, dtdz);

    FAB<Real> fx({-1,-2,-2}, {nx+1,ny+1,nz+1}, NCONS);
    FAB<Real> fy({-2,-1,-2}, {nx+1,ny+1,nz+1}, NCONS);
    FAB<Real> fz({-2,-2,-1}, {nx+1,ny+1,nz+1}, NCONS);

    FAB<Real> qgdx({-1,-2,-2}, {nx+1,ny+1,nz+1}, NGDNV);
    FAB<Real> qgdy({-2,-1,-2}, {nx+1,ny+1,nz+1}, NGDNV);
    FAB<Real> qgdz({-2,-2,-1}, {nx+1,ny+1,nz+1}, NGDNV);

    // The first flux estimates as per the corner-transport-upwind method

    // X initial fluxes
    cmpflx(0, fx, qxm, qxp, qgdx);

    // Y initial fluxes
    cmpflx(1, fy, qym, qyp, qgdy);

    // Z initial fluxes
    cmpflx(2, fz, qzm, qzp, qgdz);

    // X interface corrections
    FAB<Real> qmxy({-1,-1,-1}, {nx+1,ny,nz}, NPRIM);
    FAB<Real> qpxy({-1,-1,-1}, {nx+1,ny,nz}, NPRIM);
    FAB<Real> qmxz({-1,-1,-1}, {nx+1,ny,nz}, NPRIM);
    FAB<Real> qpxz({-1,-1,-1}, {nx+1,ny,nz}, NPRIM);

    transdo(0, 1, qmxy, qpxy, qxm, qxp, fy, qgdy, cdtdy);
    transdo(0, 2, qmxz, qpxz, qxm, qxp, fz, qgdz, cdtdz);

    FAB<Real> flxy({-1,-1,-1}, {nx+1,ny,nz}, NCONS);
    FAB<Real> flxz({-1,-1,-1}, {nx+1,ny,nz}, NCONS);
    FAB<Real>  qxy({-1,-1,-1}, {nx+1,ny,nz}, NGDNV);
    FAB<Real>  qxz({-1,-1,-1}, {nx+1,ny,nz}, NGDNV);

    // Riemann problem X|Y, X|Z
    cmpflx(0, flxy, qmxy, qpxy, qxy);
    cmpflx(0, flxz, qmxz, qpxz, qxz);

    qmxy.clear();
    qpxy.clear();
    qmxz.clear();
    qpxz.clear();

    // Y interface corrections
    FAB<Real> qmyx({-1,-1,-1}, {nx,ny+1,nz}, NPRIM);
    FAB<Real> qpyx({-1,-1,-1}, {nx,ny+1,nz}, NPRIM);
    FAB<Real> qmyz({-1,-1,-1}, {nx,ny+1,nz}, NPRIM);
    FAB<Real> qpyz({-1,-1,-1}, {nx,ny+1,nz}, NPRIM);

    transdo(1, 0, qmyx, qpyx, qym, qyp, fx, qgdx, cdtdx);
    transdo(1, 2, qmyz, qpyz, qym, qyp, fz, qgdz, cdtdz);

    fz.clear();
    qgdz.clear();

    FAB<Real> flyx({-1,-1,-1}, {nx,ny+1,nz}, NCONS);
    FAB<Real> flyz({-1,-1,-1}, {nx,ny+1,nz}, NCONS);
    FAB<Real>  qyx({-1,-1,-1}, {nx,ny+1,nz}, NGDNV);
    FAB<Real>  qyz({-1,-1,-1}, {nx,ny+1,nz}, NGDNV);

    // Riemann problem Y|X, Y|Z
    cmpflx(1, flyx, qmyx, qpyx, qyx);
    cmpflx(1, flyz, qmyz, qpyz, qyz);

    qmyx.clear();
    qpyx.clear();
    qmyz.clear();
    qpyz.clear();

    // Z interface corrections
    FAB<Real> qmzx({-1,-1,-1}, {nx,ny,nz+1}, NPRIM);
    FAB<Real> qpzx({-1,-1,-1}, {nx,ny,nz+1}, NPRIM);
    FAB<Real> qmzy({-1,-1,-1}, {nx,ny,nz+1}, NPRIM);
    FAB<Real> qpzy({-1,-1,-1}, {nx,ny,nz+1}, NPRIM);

    transdo(2, 0, qmzx, qpzx, qzm, qzp, fx, qgdx, cdtdx);
    transdo(2, 1, qmzy, qpzy, qzm, qzp, fy, qgdy, cdtdy);

    fx.clear();
    fy.clear();
    qgdx.clear();
    qgdy.clear();

    FAB<Real> flzx({-1,-1,-1}, {nx,ny,nz+1}, NCONS);
    FAB<Real> flzy({-1,-1,-1}, {nx,ny,nz+1}, NCONS);
    FAB<Real>  qzx({-1,-1,-1}, {nx,ny,nz+1}, NGDNV);
    FAB<Real>  qzy({-1,-1,-1}, {nx,ny,nz+1}, NGDNV);

    // Riemann problem Z|X, Z|Y
    cmpflx(2, flzx, qmzx, qpzx, qzx);
    cmpflx(2, flzy, qmzy, qpzy, qzy);

    qmzx.clear();
    qpzx.clear();
    qmzy.clear();
    qpzy.clear();

    FAB<Real> qm({0,0,0}, {nx,ny,nz}, NPRIM);
    FAB<Real> qp({0,0,0}, {nx,ny,nz}, NPRIM);

    // X | Y&Z
    transdd(0, qm, qp, qxm, qxp, flyz, flzy, qyz, qzy, hdtdy, hdtdz);

    qxm.clear();
    qxp.clear();
    flyz.clear();
    flzy.clear();
    qyz.clear();
    qzy.clear();

    // Final X flux
    FAB<Real> flx1({0,0,0}, {nx,ny,nz}, NCONS);
    FAB<Real>   q1({0,0,0}, {nx,ny,nz}, NGDNV);
    cmpflx(0, flx1, qm, qp, q1);

    // Y | X&Z
    transdd(1, qm, qp, qym, qyp, flxz, flzx, qxz, qzx, hdtdx, hdtdz);

    qym.clear();
    qyp.clear();
    flxz.clear();
    flzx.clear();
    qxz.clear();
    qzx.clear();

    // Final Y flux
    FAB<Real> flx2({0,0,0}, {nx,ny,nz}, NCONS);
    FAB<Real>   q2({0,0,0}, {nx,ny,nz}, NGDNV);
    cmpflx(1, flx2, qm, qp, q2);

    // Z | X&Y
    transdd(2, qm, qp, qzm, qzp, flxy, flyx, qxy, qyx, hdtdx, hdtdy);

    qzm.clear();
    qzp.clear();
    flxy.clear();
    flyx.clear();
    qxy.clear();
    qyx.clear();

    // Final Z flux
    FAB<Real> flx3({0,0,0}, {nx,ny,nz}, NCONS);
    FAB<Real>   q3({0,0,0}, {nx,ny,nz}, NGDNV);
    cmpflx(2, flx3, qm, qp, q3);

    // TODO: add artificial viscosity

    for (int n = 0; n < NCONS; ++n) {
    for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
        if (n != UEINT) {
            for (int i = 0; i < nx; ++i) {
                state(i,j,k,n) += dtdx * (flx1(i,j,k,n) - flx1(i+1,j  ,k  ,n))
                    +             dtdy * (flx2(i,j,k,n) - flx2(i  ,j+1,k  ,n))
                    +             dtdz * (flx3(i,j,k,n) - flx3(i  ,j  ,k+1,n));
            }
        } else {
            for (int i = 0; i < nx; ++i) {
                state(i,j,k,n) += dtdx * (flx1(i,j,k,n) - flx1(i+1,j  ,k  ,n))
                    +             dtdy * (flx2(i,j,k,n) - flx2(i  ,j+1,k  ,n))
                    +             dtdz * (flx3(i,j,k,n) - flx3(i  ,j  ,k+1,n))
                    + Real(0.5)*(dtdx*((q1(i,j,k,GDPRES) + q1(i+1,j  ,k  ,GDPRES))*
                                       (q1(i,j,k,GDU   ) - q1(i+1,j  ,k  ,GDU   )))+
                                 dtdy*((q2(i,j,k,GDPRES) + q2(i  ,j+1,k  ,GDPRES))*
                                       (q2(i,j,k,GDV   ) - q2(i  ,j+1,k  ,GDV   )))+
                                 dtdz*((q3(i,j,k,GDPRES) + q3(i  ,j  ,k+1,GDPRES))*
                                       (q3(i,j,k,GDW   ) - q3(i  ,j  ,k+1,GDW   ))));
            }
        }
    }}}
}

Real compute_l1_error (Real t, FAB<Real> const& prim,
                       std::vector<double>& x, std::vector<double>& p,
                       std::vector<double>& rho, std::vector<double>& v)
{
    x.resize(nx);
    p.resize(nx);
    rho.resize(nx);
    v.resize(nx);
    for (int i = 0; i < nx; ++i) {
        x[i] = (double(i) + 0.5) * double(dx);
    }

    exact_riemann(x, p, rho, v, double(t),
                  double(Sod::pl), double(Sod::rhol), double(Sod::ul),
                  double(Sod::pr), double(Sod::rhor), double(Sod::ur));

    double err = 0;
    for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
        err += std::abs(double(prim(i,j,k,QRHO)) - rho[i]);
    }}}
    return Real(err) * dx * dy * dz;
}
