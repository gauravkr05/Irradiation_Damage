// ============================================================
//  1D Cahn-Hilliard equation  (CONSERVED order parameter)
//
//      dc/dt = M * d2/dx2 [ mu ]
//      mu    = df/dc - 2*kappa*d2c/dx2
//
//  Bulk free energy :  f(c) = A c^2 (1 - c)^2
//      df/dc = A ( 2*c - 6*c^2 + 4*c^3 )
//
//  Parameters from the sheet :  kappa = 1 , M = 1 , A = 1
//
//  I.C. :  background c = 0.2 everywhere,
//          with a block of c = 1.0 of width 50, centred in the domain.
//          Domain length = 512.
//
//  NOTE : Cahn-Hilliard IS conserved -- sum(c) can never change.
//         The 0.2 background is NOT an equilibrium phase (the double well
//         f = A c^2 (1-c)^2 has its common tangent along the c-axis, so the
//         only coexisting compositions are c = 0 and c = 1). So the
//         background must drain away to c = 0, but the solute it carries
//         cannot disappear -- it is pumped into the block, which widens.
//
//         Average composition = (50*1.0 + 462*0.2)/512 = 0.278
//         => final state: one block of c=1 about 0.278*512 = 142 cells wide,
//            sitting in a background of c=0.
//
//         This is the whole contrast with Allen-Cahn, where the same 0.2
//         background just decays locally to 0 and the block stays 50 wide.
//
//         Implemented as two passes (mu, then lap(mu)) so only the 3-point
//         laplacian is ever needed -- no explicit 5-point biharmonic.
// ============================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/stat.h>
using namespace std;

int main()
{
    // ---------------- Input ----------------
    double M     = 1.0;        // mobility
    double kappa = 1.0;        // gradient energy coefficient
    double A     = 1.0;        // free energy barrier height

    double dx = 1.0;
    double dt = 0.01;
    int    nx = 512;                               // domain 0 .. 512
    double total_time = 400.0;
    long   total_timesteps = (long)(total_time / dt + 0.5);
    int    n_saves = 5;
    int    save_times[5] = {0, 1, 5, 200, 400};

    double c_background = 0.2;                     // starting value outside the block
    double c_block      = 1.0;                     // value inside the block
    int    block_width  = 50;                      // the "gap" of 50
    int    block_lo = (nx - block_width) / 2;      // centred: 231
    int    block_hi = block_lo + block_width;      //          281

    const string tag = "ch";

    // ---------------- Stability check ----------------
    // Linearised growth rate of mode q :  M * q2 * ( f'' + 2*kappa*q2 )
    // Discrete laplacian symbol q2 = (2 - 2 cos(k dx))/dx^2  in [0, 4/dx^2].
    // Worst case uses |f''|max = 2A (attained at c = 0 and c = 1).
    double q2max  = 4.0 / (dx * dx);
    double lam    = M * q2max * (2.0 * A + 2.0 * kappa * q2max);
    cout << "Cahn-Hilliard : M = " << M << ", kappa = " << kappa << ", A = " << A
         << ", dt = " << dt << ", nx = " << nx << "\n";
    cout << "dt * lambda_max = " << dt * lam << "  (needs <= 2)\n";
    cout << "block : i = " << block_lo << " .. " << block_hi
         << "  (width " << block_width << ")\n";
    if (dt * lam > 2.0) { cout << "!! UNSTABLE : reduce dt\n"; return 1; }

    // ---------------- Variable Decl ----------------
    vector<double> c    (nx, 0.0);
    vector<double> c_new(nx, 0.0);
    vector<double> mu   (nx, 0.0);

    // ---------------- I.C. ----------------
    for (int i = 0; i < nx; i++)
    {
        if (i >= block_lo && i < block_hi)  c[i] = c_block;
        else                                c[i] = c_background;
    }

    // ---------------- Time loop ----------------
    for (long t = 0; t <= total_timesteps; t++)
    {
        if (t > 0)
        {
            // ---- pass 1 : chemical potential mu = df/dc - 2*kappa*lap(c)
            for (int i = 0; i < nx; i++)
            {
                // B.C. : periodic
                int im1 = (i - 1 + nx) % nx;
                int ip1 = (i + 1) % nx;

                double ci = c[i];

                // df/dc
                double d1f = A * (2.0*ci - 6.0*ci*ci + 4.0*ci*ci*ci);

                // central difference laplacian of c
                double lap_c = (c[ip1] - 2.0*ci + c[im1]) / (dx*dx);

                mu[i] = d1f - 2.0*kappa*lap_c;
            }

            // ---- pass 2 : dc/dt = M * lap(mu)
            for (int i = 0; i < nx; i++)
            {
                int im1 = (i - 1 + nx) % nx;
                int ip1 = (i + 1) % nx;

                double lap_mu = (mu[ip1] - 2.0*mu[i] + mu[im1]) / (dx*dx);

                // D.E.
                c_new[i] = c[i] + dt * M * lap_mu;
            }
            c.swap(c_new);
        }

        // ---------------- Output ----------------
        double time = t * dt;
        for (int k = 0; k < n_saves; k++)
        {
            if (t == (long)(save_times[k] / dt + 0.5))
            {
                ostringstream folder;
                folder << tag << "_time_" << save_times[k];
                mkdir(folder.str().c_str(), 0777);

                ostringstream fname;
                fname << folder.str() << "/c_time_" << save_times[k] << ".dat";
                ofstream out(fname.str().c_str());
                out << "x c\n";
                for (int i = 0; i < nx; i++)
                    out << i*dx << " " << c[i] << "\n";
                out.close();

                // mass check : this number must not move (that is the point)
                double total = 0.0;
                for (int i = 0; i < nx; i++) total += c[i]*dx;

                cout << "saved time = " << time
                     << "   sum(c)dx = " << total
                     << "   <c> = " << total/(nx*dx) << "\n";
            }
        }
    }
    return 0;
}
