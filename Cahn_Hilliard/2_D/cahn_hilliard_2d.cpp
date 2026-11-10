// ============================================================
//  2D Cahn-Hilliard equation  (explicit / forward Euler)
//
//      dc/dt = M grad^2 ( df/dc - 2*kappa*grad^2 c )
//
//  Expanded, exactly the same idea as the 1D version:
//
//      dc/dt = M (d3f/dc3) |grad c|^2 + M (d2f/dc2) grad^2 c
//              - 2*M*kappa*grad^4 c
//
//  where   |grad c|^2 = (dc/dx)^2 + (dc/dy)^2
//          grad^2 c   = d2c/dx2 + d2c/dy2                 (5-point)
//          grad^4 c   = grad^2 ( grad^2 c )               (laplacian applied twice)
//
//  Bulk free energy :  f(c) = A c^2 (1-c)^2
//      df/dc   = A ( 2c - 6c^2 + 4c^3 )
//      d2f/dc2 = A ( 2 - 12c + 12c^2 )
//      d3f/dc3 = A ( -12 + 24c )
//
//  I.C. :  c = 0.5 + noise  at EVERY grid point (i,j),
//          noise = 0.5*(1/100)*(2R-1),  R random in [0,1]
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
    double dy = 1.0;
    double dt = 0.002;
    int    nx = 100;
    int    ny = 100;
    long   total_timesteps = 400000;               // total time = 800
    int    n_saves = 4;
    int    save_times[4] = {0, 50, 200, 800};

    unsigned int seed = 42;    // fixed seed = reproducible; use time(NULL) for a new one
    const string tag = "ch2d";

    // set true to use the conservative form instead
    const bool USE_CONSERVATIVE = false;

    cout << "2D Cahn-Hilliard : M = " << M << ", kappa = " << kappa
         << ", A = " << A << ", dt = " << dt
         << ", grid = " << nx << "x" << ny
         << ", total time = " << total_timesteps*dt << "\n";

    // ---------------- Variable Decl ----------------
    vector<vector<double> > c    (nx, vector<double>(ny, 0.0));
    vector<vector<double> > c_new(nx, vector<double>(ny, 0.0));
    vector<vector<double> > lap  (nx, vector<double>(ny, 0.0));   // grad^2 c
    vector<vector<double> > mu   (nx, vector<double>(ny, 0.0));   // conservative form only

    // ---------------- I.C. : 0.5 + noise on every grid point ----------------
    srand(seed);
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
        {
            double R = (double)rand() / (double)RAND_MAX;             // R in [0,1]
            double noise = 0.5 * (1.0 / 100.0) * (2.0 * R - 1.0);     // +/- 0.005
            c[i][j] = 0.5 + noise;
        }

    // ---------------- Time loop ----------------
    for (long t = 0; t <= total_timesteps; t++)
    {
        if (t > 0)
        {
            if (!USE_CONSERVATIVE)
            {
                // ---- PASS 1 : laplacian of c at every point ----
                for (int i = 0; i < nx; i++)
                {
                    int im1 = (i - 1 + nx) % nx;
                    int ip1 = (i + 1) % nx;
                    for (int j = 0; j < ny; j++)
                    {
                        int jm1 = (j - 1 + ny) % ny;
                        int jp1 = (j + 1) % ny;
                        lap[i][j] = (c[ip1][j] - 2.0*c[i][j] + c[im1][j]) / (dx*dx)
                                  + (c[i][jp1] - 2.0*c[i][j] + c[i][jm1]) / (dy*dy);
                    }
                }
                // ---- PASS 2 : update ----
                for (int i = 0; i < nx; i++)
                {
                    int im1 = (i - 1 + nx) % nx;
                    int ip1 = (i + 1) % nx;
                    for (int j = 0; j < ny; j++)
                    {
                        int jm1 = (j - 1 + ny) % ny;
                        int jp1 = (j + 1) % ny;

                        double ci = c[i][j];

                        // free energy derivatives
                        double d2f = A * (2.0 - 12.0*ci + 12.0*ci*ci);
                        double d3f = A * (-12.0 + 24.0*ci);

                        // central differences
                        double gx = (c[ip1][j] - c[im1][j]) / (2.0*dx);
                        double gy = (c[i][jp1] - c[i][jm1]) / (2.0*dy);
                        double grad2 = gx*gx + gy*gy;                 // |grad c|^2

                        // biharmonic = laplacian of the laplacian
                        double bih = (lap[ip1][j] - 2.0*lap[i][j] + lap[im1][j]) / (dx*dx)
                                   + (lap[i][jp1] - 2.0*lap[i][j] + lap[i][jm1]) / (dy*dy);

                        // D.E.
                        c_new[i][j] = ci + dt * ( M*d3f*grad2
                                                + M*d2f*lap[i][j]
                                                - 2.0*M*kappa*bih );
                    }
                }
            }
            else
            {
                // ---- PASS 1 : chemical potential mu ----
                for (int i = 0; i < nx; i++)
                {
                    int im1 = (i - 1 + nx) % nx;
                    int ip1 = (i + 1) % nx;
                    for (int j = 0; j < ny; j++)
                    {
                        int jm1 = (j - 1 + ny) % ny;
                        int jp1 = (j + 1) % ny;
                        double ci = c[i][j];
                        double d1f = A * (2.0*ci - 6.0*ci*ci + 4.0*ci*ci*ci);
                        double L = (c[ip1][j] - 2.0*ci + c[im1][j]) / (dx*dx)
                                 + (c[i][jp1] - 2.0*ci + c[i][jm1]) / (dy*dy);
                        mu[i][j] = d1f - 2.0*kappa*L;
                    }
                }
                // ---- PASS 2 : laplacian of mu ----
                for (int i = 0; i < nx; i++)
                {
                    int im1 = (i - 1 + nx) % nx;
                    int ip1 = (i + 1) % nx;
                    for (int j = 0; j < ny; j++)
                    {
                        int jm1 = (j - 1 + ny) % ny;
                        int jp1 = (j + 1) % ny;
                        double Lmu = (mu[ip1][j] - 2.0*mu[i][j] + mu[im1][j]) / (dx*dx)
                                   + (mu[i][jp1] - 2.0*mu[i][j] + mu[i][jm1]) / (dy*dy);
                        c_new[i][j] = c[i][j] + dt * M * Lmu;
                    }
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
                out << "x y c\n";
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        out << i*dx << " " << j*dy << " " << c[i][j] << "\n";
                out.close();

                cout << "saved time = " << time << "\n";
            }
        }
    }
    return 0;
}
