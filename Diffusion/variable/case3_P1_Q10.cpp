// ============================================================
//  1D Non-linear Diffusion  :  dc/dt = d/dx [ D(c) dc/dx ]
//  with  D(c) = P + Q*c^2
//
//  CASE 3 :  P = 1.0 ,  Q = 10.0    (D goes from 1 at c=0  to  11 at c=1)
//
//  Expanded (non-conservative) form used, exactly as derived by hand:
//      dc/dt = D(c)*d2c/dx2 + (dD/dc)*(dc/dx)^2
//      D(c)    = P + Q c^2
//      dD/dc   = 2 Q c
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
    const double P = 1.0;              // D(c) = P + Q*c^2
    const double Q = 10.0;
    const string tag = "P1_Q10";       // prefix for output folders

    double dx = 1.0;
    double dt = 0.02;                  // dt = 0.25 BLOWS UP here: Dmax = 11 needs dt < 0.0454
    int    nx = 100;
    double total_time = 1000.0;
    int    total_timesteps = (int)(total_time / dt + 0.5);
    int    save_times[4] = {0, 200, 500, 1000};

    // set to true to use the conservative flux form instead (mass-exact)
    const bool USE_CONSERVATIVE = false;

    // ---------------- Stability check ----------------
    double Dmax = P + Q;               // c is in [0,1], so D is largest at c = 1
    double r = Dmax * dt / (dx * dx);
    cout << "D(c) = " << P << " + " << Q << "*c^2 ,  dt = " << dt
         << " ,  Dmax*dt/dx^2 = " << r << "  (needs <= 0.5)\n";
    if (r > 0.5) {
        cout << "!! UNSTABLE : reduce dt below " << 0.5*dx*dx/Dmax << "\n";
        return 1;
    }

    // ---------------- Variable Decl ----------------
    vector<double> c(nx, 0.0);
    vector<double> c_new(nx, 0.0);

    // ---------------- I.C. ----------------
    for (int i = 0; i < nx; i++)
    {
        if (i < nx / 4)            c[i] = 0.0;
        else if (i < 3 * nx / 4)   c[i] = 1.0;
        else                       c[i] = 0.0;
    }

    // ---------------- Time loop ----------------
    for (int t = 0; t <= total_timesteps; t++)
    {
        if (t > 0)
        {
            for (int i = 0; i < nx; i++)
            {
                // B.C. : periodic
                int im1 = (i == 0)      ? nx - 1 : i - 1;
                int ip1 = (i == nx - 1) ? 0      : i + 1;

                if (!USE_CONSERVATIVE)
                {
                    // ---- D.E. : expanded form ----
                    double D    = P + Q * c[i] * c[i];          // D(c)
                    double dDdc = 2.0 * Q * c[i];               // dD/dc

                    double lap  = (c[ip1] - 2.0 * c[i] + c[im1]) / (dx * dx);
                    double grad = (c[ip1] - c[im1]) / (2.0 * dx);   // central diff

                    c_new[i] = c[i] + dt * (D * lap + dDdc * grad * grad);
                }
                else
                {
                    // ---- D.E. : conservative flux form ----
                    double c_p = 0.5 * (c[i] + c[ip1]);         // c at i+1/2
                    double c_m = 0.5 * (c[i] + c[im1]);         // c at i-1/2
                    double D_p = P + Q * c_p * c_p;
                    double D_m = P + Q * c_m * c_m;

                    c_new[i] = c[i] + (dt / (dx * dx)) *
                               (D_p * (c[ip1] - c[i]) - D_m * (c[i] - c[im1]));
                }
            }
            c.swap(c_new);
        }

        // ---------------- Output ----------------
        double time = t * dt;
        for (int k = 0; k < 4; k++)
        {
            if (t == (int)(save_times[k] / dt + 0.5))
            {
                ostringstream folder;
                folder << tag << "_time_" << save_times[k];
                mkdir(folder.str().c_str(), 0777);

                ostringstream fname;
                fname << folder.str() << "/c_time_" << save_times[k] << ".dat";
                ofstream out(fname.str().c_str());
                out << "x c\n";
                for (int i = 0; i < nx; i++)
                    out << i * dx << " " << c[i] << "\n";
                out.close();

                cout << "saved time = " << time << "\n";
            }
        }
    }
    return 0;
}
