#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
using namespace std;

int main()
{
    // Input
    double D = 1.0;
    double dx = 1.0;
    double dy = 1.0;
    double dt = 0.25;
    int nx = 100;
    int ny = 100;
    int total_timesteps = 4000;
    int save_times[4] = {0, 200, 500, 1000};

    // Variable Decl
    vector<vector<double> > c(nx, vector<double>(ny, 0.0));
    vector<vector<double> > c_new(nx, vector<double>(ny, 0.0));

    // I.C.
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
        {
            if (i >= nx/4 && i < 3*nx/4)
                c[i][j] = 1.0;
            else
                c[i][j] = 0.0;
        }

    for (int t = 0; t <= total_timesteps; t++)
    {
        if (t > 0)
        {
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    // B.C.
                    int im1 = (i == 0)    ? nx-1 : i-1;
                    int ip1 = (i == nx-1) ? 0    : i+1;
                    int jm1 = (j == 0)    ? ny-1 : j-1;
                    int jp1 = (j == ny-1) ? 0    : j+1;
                    // D.E
                    c_new[i][j] = c[i][j]
                        + D*dt*((c[im1][j] - 2.0*c[i][j] + c[ip1][j]) / (dx*dx)
                               +(c[i][jm1] - 2.0*c[i][j] + c[i][jp1]) / (dy*dy));
                }
            }
            c.swap(c_new);
        }

        double time = t * dt;
        for (int k = 0; k < 4; k++)
        {
            if (t == (int)(save_times[k] / dt + 0.5))
            {
                // save to time2d_X folder to avoid clash with 1D output
                ostringstream folder;
                folder << "time2d_" << save_times[k];
                mkdir(folder.str().c_str(),0777);
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