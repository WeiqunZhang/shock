
#include "FAB.H"
#include "Hydro.H"

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

int main (int argc, char* argv[])
{
    constexpr int N_CELL_YZ = 8;

    if (argc != 2) {
        std::cerr << "usage: " << argv[0] << " n_cell_x\n";
        return 1;
    }

    int n_cell_x = 0;
    try {
        n_cell_x = std::stoi(argv[1]);
    } catch (std::exception const&) {
        std::cerr << "usage: " << argv[0] << " n_cell_x\n";
        return 1;
    }

    if (n_cell_x <= 0) {
        std::cerr << "usage: " << argv[0] << " n_cell_x\n";
        return 1;
    }

    std::array<int,3> const lo{0,0,0};
    std::array<int,3> const hi{n_cell_x-1, N_CELL_YZ-1, N_CELL_YZ-1};

    std::array<int,3> glo, ghi;
    for (int i = 0; i < 3; ++i) {
        glo[i] = lo[i] - NUM_GROW;
        ghi[i] = hi[i] + NUM_GROW;
    }

    FAB<Real> state(glo, ghi, NCONS);
    FAB<Real> prim(glo, ghi, NPRIM);

    // Set up Sod shocktube problem
    init_sod(state);

    Real t = 0;
    Real const stop_time = Real(0.2);
    Real const cfl = Real(0.8);

    for (int istep = 0; istep < 1000000; ++istep) {
        // compute dt
        Real dt = compute_dt(state);
        if (istep == 0) {
            dt *= 0.2;
        } else {
            dt *= cfl;
        }

        Real tnew = t + dt;
        bool to_break = false;
        if ((tnew + Real(1.0e-4)*dt) >= stop_time) {
            tnew = stop_time;
            dt = tnew - t;
            to_break = true;
        }

        std::cout << " Step " << istep << ": t = " << t << ", dt = " << dt << '\n';

        fill_boundary(state);

        cons_to_prim(state, prim);

        godunov(state, prim, dt);

        t = tnew;
        if (to_break) { break; }
    }

    std::cout << " Stop at t = " << t << '\n';

    cons_to_prim(state, prim);

    std::vector<double> xexact, pexact, rhoexact, vexact;
    Real l1_error =  compute_l1_error(t, prim, xexact, pexact, rhoexact, vexact);

    std::cout << " Number of cells = " << n_cell_x
              << ", L1 error of rho = " << l1_error << '\n';

    std::ofstream ofs("output.txt");
    int const w_i = 5;
    int const w_val = 14;
    ofs << std::setw(w_i) << "i"
        << std::setw(w_val) << "rho"
        << std::setw(w_val) << "v"
        << std::setw(w_val) << "p"
        << std::setw(w_val) << "rho_exact"
        << std::setw(w_val) << "v_exact"
        << std::setw(w_val) << "p_exact"
        << '\n';
    ofs << std::fixed << std::setprecision(6);
    for (int i = 0; i <= hi[0]; ++i) {
        ofs << std::setw(w_i) << i
            << std::setw(w_val) << prim(i,0,0,QRHO)
            << std::setw(w_val) << prim(i,0,0,QU)
            << std::setw(w_val) << prim(i,0,0,QPRES)
            << std::setw(w_val) << rhoexact[i]
            << std::setw(w_val) << vexact[i]
            << std::setw(w_val) << pexact[i]
            << '\n';
    }
}
