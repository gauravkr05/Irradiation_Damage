// ============================================================
//  Shrinking of a circular grain in an infinite matrix (Chen 1995, Fig. 5 IC)
//  Multi-order-parameter Allen-Cahn, 9-point Laplacian, explicit Euler.
//
//  Local free energy density:
//      f = sum_i ( -alpha/2 eta_i^2 + beta/4 eta_i^4 )
//          + gamma sum_i sum_{j>i} eta_i^2 eta_j^2      (2p wells if gamma>beta/2)
//  Kinetics:
//      d(eta_i)/dt = -L [ -alpha eta_i + beta eta_i^3
//                         + 2 gamma eta_i sum_{j!=i} eta_j^2 - kappa lap(eta_i) ]
//
//  I.C.: matrix eta_1=1, circle eta_2=1 (radius R0, centred). Periodic BCs.
//  Curvature-driven: v = L kappa / R  ->  R^2(t)=R0^2-2 L kappa t (area linear).
//  Writes f5_area.dat (step time area) and eta snapshots at selected steps.
//
//  Build : g++ -O3 -o f5 f5_circle_shrink.cpp
//  Run   : ./f5
//  Plot  : python plot_area.py f5_area.dat
// ============================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/stat.h>
using namespace std;

const int P  = 2;
const int NX = 200;
const int NY = 200;

inline int id(int i, int j) { return i * NY + j; }

// ---- 9-point Laplacian, periodic B.C. ----
inline double laplacian(const vector<double> &e, int i, int j, double dx)
{
    int ip = (i + 1) % NX, im = (i - 1 + NX) % NX;
    int jp = (j + 1) % NY, jm = (j - 1 + NY) % NY;
    double c = e[id(i, j)];

    double nn  = e[id(ip, j )] + e[id(im, j )]
               + e[id(i , jp)] + e[id(i , jm)] - 4.0 * c;
    double nnn = e[id(ip, jp)] + e[id(ip, jm)]
               + e[id(im, jp)] + e[id(im, jm)] - 4.0 * c;

    return (0.5 * nn + 0.25 * nnn) / (dx * dx);
}

int main()
{
    // ---------------- Input ----------------
    double alpha = 1.0;
    double beta  = 1.0;
    double gam   = 1.0;
    double kappa = 2.0;
    double L     = 1.0;

    double dx = 2.0;
    double dt = 0.1;
    long   total_timesteps = 5200;                 // carried past the disappearance

    const int n_saves = 6;
    long   save_steps[n_saves] = {200, 2000, 4000, 4600, 4800, 5000};

    double R0 = 41.0;              // initial radius (= 20.5 dx)
    double xc = 0.5 * NX * dx;     // centre of the cell
    double yc = 0.5 * NY * dx;

    const string tag = "f5";

    // ---------------- Stability check ----------------
    double r = 3.0 * kappa * L * dt / (dx * dx);
    cout << "Fig.5 : alpha=" << alpha << " beta=" << beta << " gamma=" << gam
         << " kappa=" << kappa << " L=" << L << "\n";
    cout << "dx=" << dx << " dt=" << dt << " grid=" << NX << "x" << NY
         << " R0=" << R0 << " (" << R0 / dx << " lattice spacings)\n";
    cout << "3*kappa*L*dt/dx^2 = " << r << "  (needs < 1)\n";
    if (gam <= 0.5 * beta) { cout << "!! gamma <= beta/2 : no 2p wells\n"; return 1; }
    if (r >= 1.0)          { cout << "!! UNSTABLE : reduce dt\n";          return 1; }

    // ---------------- Variable Decl ----------------
    vector< vector<double> > eta    (P, vector<double>(NX * NY, 0.0));
    vector< vector<double> > eta_new(P, vector<double>(NX * NY, 0.0));

    // ---------------- I.C. : circular grain in a matrix ----------------
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
        {
            double X = i * dx - xc;
            double Y = j * dx - yc;
            bool inside = (X * X + Y * Y <= R0 * R0);
            if (inside) { eta[0][id(i, j)] = 0.0; eta[1][id(i, j)] = 1.0; }  // grain
            else        { eta[0][id(i, j)] = 1.0; eta[1][id(i, j)] = 0.0; }  // matrix
        }

    ofstream area_out((tag + "_area.dat").c_str());
    area_out << "step time area\n";

    long vanish_step = -1;

    // ---------------- Time loop ----------------
    for (long t = 0; t <= total_timesteps; t++)
    {
        if (t > 0)
        {
            for (int p = 0; p < P; p++)
                for (int i = 0; i < NX; i++)
                    for (int j = 0; j < NY; j++)
                    {
                        int    k  = id(i, j);
                        double ei = eta[p][k];

                        double sum_other = 0.0;
                        for (int q = 0; q < P; q++)
                            if (q != p) sum_other += eta[q][k] * eta[q][k];

                        double d1f = -alpha * ei + beta * ei * ei * ei
                                   + 2.0 * gam * ei * sum_other;
                        double lap = laplacian(eta[p], i, j, dx);
                        eta_new[p][k] = ei + dt * ( -L * (d1f - kappa * lap) );
                    }
            for (int p = 0; p < P; p++) eta[p].swap(eta_new[p]);
        }

        // ---- area of the shrinking grain : cells where eta_2 dominates ----
        if (t % 20 == 0)
        {
            long ncell = 0;
            for (int i = 0; i < NX; i++)
                for (int j = 0; j < NY; j++)
                    if (eta[1][id(i, j)] > eta[0][id(i, j)]) ncell++;

            area_out << t << " " << t * dt << " " << ncell * dx * dx << "\n";

            if (ncell == 0 && vanish_step < 0)
            {
                vanish_step = t;
                cout << "grain disappeared at step " << t
                     << " (time " << t * dt << ")\n";
            }
        }

        // ---------------- snapshots ----------------
        for (int k = 0; k < n_saves; k++)
            if (t == save_steps[k])
            {
                ostringstream folder;
                folder << tag << "_step_" << save_steps[k];
                mkdir(folder.str().c_str(), 0777);

                ostringstream fname;
                fname << folder.str() << "/eta_step_" << save_steps[k] << ".dat";
                ofstream out(fname.str().c_str());
                out << "# x y eta1 eta2\n";
                for (int i = 0; i < NX; i++)
                {
                    for (int j = 0; j < NY; j++)
                        out << i << " " << j << " "
                            << eta[0][id(i, j)] << " " << eta[1][id(i, j)] << "\n";
                    out << "\n";
                }
                out.close();
                cout << "saved step " << t << " (time " << t * dt << ")\n";
            }
    }

    area_out.close();
    if (vanish_step < 0)
        cout << "grain still present at last step -- increase total_timesteps\n";
    cout << "area history written to " << tag << "_area.dat\n";
    return 0;
}