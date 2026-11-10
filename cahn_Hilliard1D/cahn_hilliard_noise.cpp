// ============================================================
//  1D Cahn-Hilliard equation  (explicit / forward Euler)
//
//      dc/dt = M d2/dx2 ( df/dc - 2*kappa*d2c/dx2 )
//
//  Expanded exactly as on the sheet:
//
//      dc/dt = M (d3f/dc3)(dc/dx)^2 + M (d2f/dc2)(d2c/dx2)
//              - 2*M*kappa*(d4c/dx4)
//
//  Bulk free energy :  f(c) = A c^2 (1-c)^2
//      df/dc   = A ( 2c - 6c^2 + 4c^3 )
//      d2f/dc2 = A ( 2 - 12c + 12c^2 )
//      d3f/dc3 = A ( -12 + 24c )
//
//  I.C. :  c = 0.5 + noise,   noise = 0.5*(1/100)*(2R-1),  R in [0,1]
//          -> spinodal decomposition
// ============================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>          // POSIX mkdir (NOT <direct.h>, that is Windows-only)
using namespace std;

int main()
{
    // ---------------- Input ----------------
    double M     = 1.0;        // mobility
    double kappa = 1.0;        // gradient energy coefficient
    double A     = 1.0;        // free energy barrier height

    double dx = 1.0;
    double dt = 0.002;
    int    nx = 512;
    long   total_timesteps = 400000;              // 4 lakh -> total time = 800
    int    n_saves = 7;
    int    save_times[7] = {0, 10, 50, 100, 200, 400, 800};

    unsigned int seed = 42;    // fixed seed = reproducible run; use time(NULL) for a new one
    const string tag = "ch_noise";

    // set true to use the conservative form instead (see note at bottom)
    const bool USE_CONSERVATIVE = false;

    cout << "nx = " << nx << ", M = " << M << ", kappa = " << kappa << ", A = " << A
         << ", dt = " << dt << ", total time = " << total_timesteps*dt << "\n";

    // ---------------- Variable Decl ----------------
    vector<double> c(nx, 0.0);
    vector<double> c_new(nx, 0.0);
    vector<double> mu(nx, 0.0);          // only used by the conservative form

    // ---------------- I.C. : 0.5 + noise ----------------
    srand(seed);
    for (int i = 0; i < nx; i++)
    {
        double R = (double)rand() / (double)RAND_MAX;      // R in [0,1]
        double noise = 0.5 * (1.0 / 100.0) * (2.0 * R - 1.0);   // +/- 0.005
        c[i] = 0.5 + noise;
    }

    // ---------------- Time loop ----------------
    for (long t = 0; t <= total_timesteps; t++)
    {
        if (t > 0)
        {
            if (!USE_CONSERVATIVE)
            {
                for (int i = 0; i < nx; i++)
                {
                    // B.C. : periodic, need i-2 .. i+2 for the 4th derivative
                    int im2 = (i - 2 + nx) % nx;
                    int im1 = (i - 1 + nx) % nx;
                    int ip1 = (i + 1) % nx;
                    int ip2 = (i + 2) % nx;

                    double ci = c[i];

                    // free energy derivatives
                    double d2f = A * (2.0 - 12.0*ci + 12.0*ci*ci);      // d2f/dc2
                    double d3f = A * (-12.0 + 24.0*ci);                 // d3f/dc3

                    // central differences
                    double grad = (c[ip1] - c[im1]) / (2.0*dx);
                    double lap  = (c[ip1] - 2.0*ci + c[im1]) / (dx*dx);
                    double bih  = (c[ip2] - 4.0*c[ip1] + 6.0*ci
                                   - 4.0*c[im1] + c[im2]) / (dx*dx*dx*dx);

                    // D.E.
                    c_new[i] = ci + dt * ( M*d3f*grad*grad
                                         + M*d2f*lap
                                         - 2.0*M*kappa*bih );
                }
            }
            else
            {
                // ---- conservative form : mu first, then laplacian of mu ----
                for (int i = 0; i < nx; i++)
                {
                    int im1 = (i - 1 + nx) % nx;
                    int ip1 = (i + 1) % nx;
                    double ci = c[i];
                    double d1f = A * (2.0*ci - 6.0*ci*ci + 4.0*ci*ci*ci);   // df/dc
                    double lap = (c[ip1] - 2.0*ci + c[im1]) / (dx*dx);
                    mu[i] = d1f - 2.0*kappa*lap;
                }
                for (int i = 0; i < nx; i++)
                {
                    int im1 = (i - 1 + nx) % nx;
                    int ip1 = (i + 1) % nx;
                    c_new[i] = c[i] + dt * M * (mu[ip1] - 2.0*mu[i] + mu[im1]) / (dx*dx);
                }
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

                cout << "saved time = " << time << "\n";
            }
        }
    }
    return 0;
}
