#include <iostream>
#include <fstream>
#include <sstream>
#include <direct.h>
#include <vector>
#include <cmath>
#include <sys/stat.h>
using namespace std;

int main()
{
    // Input
    double D  = 1.0;
    double dx = 1.0;
    double dt = 1.0;
    int nx = 100;
    int total_timesteps = 1000;
    int save_times[4] = {0, 200, 500, 1000};

    double r = D * dt / (dx * dx);

    // Gauss-Seidel parameters
    int    max_iter = 10000;
    double tol      = 1e-8;

    // Variable Decl
    vector<double> c(nx, 0.0);
    vector<double> c_new(nx, 0.0);
    vector<double> b(nx, 0.0);

    // Build A (nx x nx) -- tridiagonal with periodic BC corners
    vector<vector<double>> A(nx, vector<double>(nx, 0.0));
    for (int i = 0; i < nx; i++)
    {
        int im1 = (i == 0)    ? nx-1 : i-1;
        int ip1 = (i == nx-1) ? 0    : i+1;
        A[i][i]   =  1.0 + r;
        A[i][im1] = -r / 2.0;
        A[i][ip1] = -r / 2.0;
    }

    // I.C.
    for (int i = 0; i < nx; i++)
    {
        if (i < nx / 4)      c[i] = 0.0;
        else if (i < 3*nx/4) c[i] = 1.0;
        else                  c[i] = 0.0;
    }

    for (int t = 0; t <= total_timesteps; t++)
    {
        if (t > 0)
        {
            // Build b = c[i] + r/2*(c[im1] - 2*c[i] + c[ip1])
            for (int i = 0; i < nx; i++)
            {
                int im1 = (i == 0)    ? nx-1 : i-1;
                int ip1 = (i == nx-1) ? 0    : i+1;
                b[i] = c[i] + (r/2.0) * (c[im1] - 2.0*c[i] + c[ip1]);
            }

            //  solve A * c_new = b
            c_new = c;
            for (int iter = 0; iter < max_iter; iter++)
            {
                double res = 0.0;
                for (int i = 0; i < nx; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < nx; j++)
                        if (j != i) sum += A[i][j] * c_new[j];
                    double c_new_i = (b[i] - sum) / A[i][i];
                    res += (c_new_i - c_new[i]) * (c_new_i - c_new[i]);
                    c_new[i] = c_new_i;
                }
                if (sqrt(res) < tol) break;
            }
            c.swap(c_new);
        }

        double time = t * dt;
        for (int k = 0; k < 4; k++)
        {
            if (t == (int)(save_times[k] / dt + 0.5))
            {
                ostringstream folder;
                folder << "semi_implicit_time_" << save_times[k];
                mkdir(folder.str().c_str());
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
