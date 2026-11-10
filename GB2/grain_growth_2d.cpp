// ============================================================
//  FIG 2 : 2D polycrystal grain growth (multi-order-parameter Allen-Cahn)
//
//  Same continuum field model as Fig 1, extended to a real 2D microstructure.
//  p field variables eta_0..eta_{NP-1}, each one an orientation.
//
//      f = sum_i ( -alpha/2 eta_i^2 + beta/4 eta_i^4 )
//          + gamma * sum_i sum_{j>i} eta_i^2 eta_j^2          (2p wells if gamma > beta/2)
//
//      d(eta_i)/dt = -L dF/d(eta_i)
//                  = -L [ -alpha*eta_i + beta*eta_i^3
//                         + 2*gamma*eta_i*sum_{j!=i} eta_j^2
//                         - kappa * lap(eta_i) ]
//
//  Laplacian: 9-point isotropic stencil (same as Fig 1),
//      lap(eta) = 1/dx^2 [ 1/2 sum_NN (eta_j - eta_i) + 1/4 sum_diag (eta_j' - eta_i) ]
//
//  INITIAL CONDITION (this is the new part, matching the sketch):
//    1. Scatter NGRAIN seed points at random over the grid.
//    2. Voronoi tessellation: every cell belongs to its nearest seed
//       (periodic / minimum-image distance, so grains wrap across the border).
//    3. Each grain is given one order parameter; that field = 1 inside the grain.
//    4. "(1) + 0.01 RN" : add a small random noise 0.01*rand to EVERY field at
//       EVERY cell, to break the symmetry so boundaries can move.
//
//  OUTPUT:
//    fig2_out/phi_XXXXX.csv    -> phi(x,y) = sum_i eta_i^2   (grain interiors ~1,
//                                 boundaries dip below 1 -> boundaries show up dark)
//    fig2_out/grain_XXXXX.csv  -> argmax_i eta_i             (integer grain/orientation map)
//    fig2_out/eta_profiles.dat -> final line scan across the mid row (the "X" crossing)
//
//  Visualise with:  python plot_grain_growth.py fig2_out
// ============================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <cstdio>
#include <sys/stat.h>
using namespace std;

const int NP     = 25;    // number of field variables (orientations)
const int NGRAIN = 25;    // number of Voronoi seeds. Keep <= NP so each grain gets
                          // its own order parameter (no artificial coalescence).

// flattened 2D index with periodic wrap
inline int id(int i, int j, int nx, int ny)
{
    int ii = (i + nx) % nx;
    int jj = (j + ny) % ny;
    return jj * nx + ii;
}

// periodic minimum-image squared distance between (i,j) and a seed (sx,sy)
inline double pdist2(int i, int j, double sx, double sy, int nx, int ny)
{
    double dxp = fabs(i - sx); if (dxp > nx / 2.0) dxp = nx - dxp;
    double dyp = fabs(j - sy); if (dyp > ny / 2.0) dyp = ny - dyp;
    return dxp * dxp + dyp * dyp;
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
    int    nx = 256, ny = 256;
    long   total_timesteps = 5000;
    long   nprint = 100;          // write a snapshot every nprint steps

    double noise_amp = 0.01;      // the "0.01 RN" perturbation amplitude
    unsigned seed = 12345;        // RNG seed (fixed -> reproducible microstructure)

    const string outdir = "fig2_out";
    mkdir(outdir.c_str(), 0777);

    // ---------------- Stability check ----------------
    double r = Lmob * kappa * dt / (dx * dx);
    cout << "alpha=" << alpha << " beta=" << beta << " gamma=" << gam
         << " kappa=" << kappa << " L=" << Lmob << "\n";
    cout << "grid " << nx << "x" << ny << ", NP=" << NP << ", grains=" << NGRAIN << "\n";
    cout << "gamma > beta/2 ? " << (gam > 0.5 * beta ? "yes (2p wells)" : "NO") << "\n";
    cout << "L*kappa*dt/dx^2 = " << r << "  (keep <= 0.25)\n";
    if (r > 0.25)   { cout << "!! UNSTABLE : reduce dt\n"; return 1; }
    if (NGRAIN > NP){ cout << "!! NGRAIN > NP : grains will merge. Raise NP.\n"; }

    // ---------------- Variable Decl ----------------
    int N = nx * ny;
    vector< vector<double> > eta    (NP, vector<double>(N, 0.0));
    vector< vector<double> > eta_new(NP, vector<double>(N, 0.0));
    vector<double>           S(N, 0.0);   // per-cell sum_q eta_q^2 (precomputed each step)

    mt19937 rng(seed);
    uniform_real_distribution<double> U(0.0, 1.0);

    // ---------------- I.C. : Voronoi grains + noise ----------------
    // 1. random seeds
    vector<double> seedx(NGRAIN), seedy(NGRAIN);
    vector<int>    seedOP(NGRAIN);
    for (int g = 0; g < NGRAIN; g++) {
        seedx[g]  = U(rng) * nx;
        seedy[g]  = U(rng) * ny;
        seedOP[g] = g % NP;      // grain g -> order parameter (g%NP); with NGRAIN<=NP, one each
    }

    // 2. assign each cell to nearest seed (Voronoi), set that field to 1
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
        {
            int k = j * nx + i;
            int    best = 0;
            double bestd = 1e300;
            for (int g = 0; g < NGRAIN; g++) {
                double d = pdist2(i, j, seedx[g], seedy[g], nx, ny);
                if (d < bestd) { bestd = d; best = g; }
            }
            eta[seedOP[best]][k] = 1.0;

            // 3. "(1) + 0.01 RN" : small noise on every field at this cell
            for (int p = 0; p < NP; p++)
                eta[p][k] += noise_amp * U(rng);
        }

    // dump the initial grain map so you can see the starting Voronoi structure
    auto write_grain = [&](long t) {
        char fn[256]; snprintf(fn, sizeof fn, "%s/grain_%05ld.csv", outdir.c_str(), t);
        ofstream f(fn);
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int k = j * nx + i;
                int best = 0; double bv = eta[0][k];
                for (int p = 1; p < NP; p++) if (eta[p][k] > bv) { bv = eta[p][k]; best = p; }
                f << best << (i + 1 < nx ? ' ' : '\n');
            }
        }
    };
    auto write_phi = [&](long t) {
        char fn[256]; snprintf(fn, sizeof fn, "%s/phi_%05ld.csv", outdir.c_str(), t);
        ofstream f(fn);
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int k = j * nx + i;
                double phi = 0.0;
                for (int p = 0; p < NP; p++) phi += eta[p][k] * eta[p][k];
                f << phi << (i + 1 < nx ? ' ' : '\n');
            }
        }
    };
    write_grain(0);
    write_phi(0);

    // ---------------- Time loop ----------------
    for (long t = 1; t <= total_timesteps; t++)
    {
        // precompute S[k] = sum_q eta_q(k)^2   -> sum_other = S[k] - eta_p(k)^2
        for (int k = 0; k < N; k++) {
            double s = 0.0;
            for (int p = 0; p < NP; p++) s += eta[p][k] * eta[p][k];
            S[k] = s;
        }

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

                    double sum_other = S[k] - e * e;   // = sum_{q!=p} eta_q^2

                    double d1f = -alpha * e + beta * e * e * e
                                 + 2.0 * gam * e * sum_other;

                    eta_new[p][k] = e + dt * ( -Lmob * (d1f - kappa * lap) );
                }
        for (int p = 0; p < NP; p++) eta[p].swap(eta_new[p]);

        if (t % nprint == 0) {
            write_grain(t);
            write_phi(t);
            // count "alive" order parameters (roughly, surviving grains)
            int alive = 0;
            for (int p = 0; p < NP; p++) {
                double mx = 0.0;
                for (int k = 0; k < N; k++) if (eta[p][k] > mx) mx = eta[p][k];
                if (mx > 0.5) alive++;
            }
            cout << "step " << t << " / " << total_timesteps
                 << "   surviving orientations = " << alive << "\n";
        }
    }

    // ---------------- Output : line scan across the mid row (the "X" crossing) ----------------
    int jmid = ny / 2;
    ofstream out((outdir + "/eta_profiles.dat").c_str());
    out << "# x";
    for (int p = 0; p < NP; p++) out << "  eta" << p;
    out << "\n";
    for (int i = 0; i < nx; i++) {
        int k = jmid * nx + i;
        out << i * dx;
        for (int p = 0; p < NP; p++) out << " " << eta[p][k];
        out << "\n";
    }
    out.close();

    cout << "files written to " << outdir << "/\n";
    cout << "visualise:  python plot_grain_growth.py " << outdir << "\n";
    return 0;
}