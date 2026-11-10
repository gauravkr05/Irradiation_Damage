// ============================================================
//  FIG 4 : Microstructure with two straight grain boundaries (grey-level map)
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
//  OUTPUT OF THIS PROGRAM: phi_map.dat       ->  x, y, phi, dominant orientation over the whole cell
//  Plot it with: gnuplot plot_fig4.gp
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

    const string outdir = "fig4_out";
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

    // ---------------- Output : 2D map of phi(r) for the grey-level microstructure ----------------
    //   Column 3 is phi (bright grains, dark boundaries).
    //   Column 4 is the dominant orientation, useful when the interface is still
    //   sharp (phi = 1 everywhere) and the phi map would look blank.
    //   The map is shifted by nx/4 so that BOTH boundaries sit inside the frame
    //   instead of one being split across the periodic edge.
    int shift = nx / 4;
    ofstream out((outdir + "/phi_map.dat").c_str());
    out << "# x y phi orientation\n";
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            int is = (i + shift) % nx;
            int k = j * nx + is;
            double phi = 0.0; int best = 0; double bmax = -1.0;
            for (int p = 0; p < NP; p++)
            {
                phi += eta[p][k] * eta[p][k];
                if (fabs(eta[p][k]) > bmax) { bmax = fabs(eta[p][k]); best = p; }
            }
            out << i * dx << " " << j * dx << " " << phi << " " << best << "\n";
        }
        out << "\n";                          // blank line -> gnuplot pm3d row break
    }
    out.close();
    cout << "map is shifted by " << shift * dx << " in x so both boundaries are visible\n";

    cout << "files written to " << outdir << "/\n";
    return 0;
}
