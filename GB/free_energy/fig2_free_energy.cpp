// ============================================================
//  FIG 2 : Local free energy density f(eta_1, eta_2) across the boundary
//
//  Model: continuum field (multi-order-parameter) description of a
//  polycrystal. p field variables eta_1..eta_p, each one an orientation.
//
//      f = sum_i ( -alpha/2 eta_i^2 + beta/4 eta_i^4 )
//          + gamma * sum_i sum_{j>i} eta_i^2 eta_j^2          (2p wells if gamma > beta/2)
//
//      d(eta_i)/dt = -L dF/d(eta_i)
//                  = -L [ -alpha*eta_i + beta*eta_i^3
//                         + 2*gamma*eta_i*sum_{j!=i} eta_j^2
//                         - kappa * lap(eta_i) ]
//
//  Laplacian: 9-point isotropic stencil,
//      lap(eta) = 1/dx^2 [ 1/2 sum_NN (eta_j - eta_i) + 1/4 sum_diag (eta_j' - eta_i) ]
//
//  Setup: 200x200 cell, periodic BCs, left half eta_1 = 1, right half
//         eta_2 = 1, relaxed for 5000 steps. Two straight boundaries result
//         (x = 200 in the middle, x = 0 at the cell border).
//
//  OUTPUT OF THIS PROGRAM: free_energy.dat   ->  x, f(x), f(x)-f_min along a line scan
//  Plot it with: gnuplot plot_fig2.gp
// ============================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/stat.h>
using namespace std;

const int NP = 2;                 // number of field variables (orientations)

// flattened 2D index with periodic wrap
inline int id(int i, int j, int nx, int ny)
{
    int ii = (i + nx) % nx;
    int jj = (j + ny) % ny;
    return jj * nx + ii;
}

int main()
{
    // ---------------- Input ----------------
    double alpha = 1.0;           // phenomenological parameters of f
    double beta  = 1.0;
    double gam   = 1.0;           // needs gam > beta/2
    double kappa = 2.0;           // gradient energy coefficient
    double Lmob  = 1.0;           // kinetic coefficient (mobility)

    double dx = 2.0;
    double dt = 0.1;
    int    nx = 200, ny = 200;
    long   total_timesteps = 5000;

    const string outdir = "fig2_out";
    mkdir(outdir.c_str(), 0777);

    // ---------------- Stability check ----------------
    double r = Lmob * kappa * dt / (dx * dx);
    cout << "alpha = " << alpha << ", beta = " << beta << ", gamma = " << gam
         << ", kappa = " << kappa << ", L = " << Lmob << "\n";
    cout << "gamma > beta/2 ? " << (gam > 0.5 * beta ? "yes (2p wells)" : "NO") << "\n";
    cout << "L*kappa*dt/dx^2 = " << r << "  (keep <= 0.25)\n";
    if (r > 0.25) { cout << "!! UNSTABLE : reduce dt\n"; return 1; }

    // ---------------- Variable Decl ----------------
    int N = nx * ny;
    vector< vector<double> > eta    (NP, vector<double>(N, 0.0));
    vector< vector<double> > eta_new(NP, vector<double>(N, 0.0));

    // ---------------- I.C. ----------------
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
        {
            int k = j * nx + i;
            if (i < nx / 2) { eta[0][k] = 1.0; eta[1][k] = 0.0; }
            else            { eta[0][k] = 0.0; eta[1][k] = 1.0; }
        }

    // ---------------- Time loop ----------------
    for (long t = 1; t <= total_timesteps; t++)
    {
        for (int p = 0; p < NP; p++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                {
                    int k = j * nx + i;
                    double e = eta[p][k];

                    double sNN = eta[p][id(i+1,j  ,nx,ny)] + eta[p][id(i-1,j  ,nx,ny)]
                               + eta[p][id(i  ,j+1,nx,ny)] + eta[p][id(i  ,j-1,nx,ny)]
                               - 4.0 * e;
                    double sD  = eta[p][id(i+1,j+1,nx,ny)] + eta[p][id(i-1,j+1,nx,ny)]
                               + eta[p][id(i+1,j-1,nx,ny)] + eta[p][id(i-1,j-1,nx,ny)]
                               - 4.0 * e;
                    double lap = (0.5 * sNN + 0.25 * sD) / (dx * dx);

                    double sum_other = 0.0;
                    for (int q = 0; q < NP; q++)
                        if (q != p) sum_other += eta[q][k] * eta[q][k];

                    double d1f = -alpha * e + beta * e * e * e
                                 + 2.0 * gam * e * sum_other;

                    eta_new[p][k] = e + dt * ( -Lmob * (d1f - kappa * lap) );
                }
        for (int p = 0; p < NP; p++) eta[p].swap(eta_new[p]);
    }

    // ---------------- Output : local free energy density along a line scan ----------------
    //   f = sum_i ( -alpha/2 eta_i^2 + beta/4 eta_i^4 ) + gamma sum_{i<j} eta_i^2 eta_j^2
    //   f_min = -alpha/2 + beta/4  inside a grain (eta = 1, others 0)
    int jmid = ny / 2;
    double fmin = -0.5 * alpha + 0.25 * beta;

    ofstream out((outdir + "/free_energy.dat").c_str());
    out << "# x  f  f-fmin\n";
    for (int i = 0; i < nx; i++)
    {
        int k = jmid * nx + i;
        double f = 0.0;
        for (int p = 0; p < NP; p++)
        {
            double e = eta[p][k];
            f += -0.5 * alpha * e * e + 0.25 * beta * e * e * e * e;
        }
        for (int p = 0; p < NP; p++)
            for (int q = p + 1; q < NP; q++)
                f += gam * eta[p][k] * eta[p][k] * eta[q][k] * eta[q][k];

        out << i * dx << " " << f << " " << f - fmin << "\n";
    }
    out.close();

    // grain boundary energy, sigma = Int [ (f - fmin) + sum_i kappa/2 (d eta_i/dx)^2 ] dx
    // integrated over the left half of the cell, which holds only the middle boundary
    double sigma = 0.0;
    for (int i = nx / 4; i < 3 * nx / 4; i++)
    {
        int k = jmid * nx + i;
        double f = 0.0, grad2 = 0.0;
        for (int p = 0; p < NP; p++)
        {
            double e = eta[p][k];
            f += -0.5 * alpha * e * e + 0.25 * beta * e * e * e * e;
            double de = (eta[p][id(i+1,jmid,nx,ny)] - eta[p][id(i-1,jmid,nx,ny)]) / (2.0*dx);
            grad2 += 0.5 * kappa * de * de;
        }
        for (int p = 0; p < NP; p++)
            for (int q = p + 1; q < NP; q++)
                f += gam * eta[p][k] * eta[p][k] * eta[q][k] * eta[q][k];
        sigma += ((f - fmin) + grad2) * dx;
    }
    cout << "f_min = " << fmin << "\n";
    cout << "grain boundary energy sigma_12 = " << sigma << "\n";

    cout << "files written to " << outdir << "/\n";
    return 0;
}
