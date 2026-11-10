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
    double dt = 0.25;
    int nx = 100;
    int total_timesteps = 4000;
    int save_times[4] = {0, 200, 500, 1000};
 
    // Variable Decl
    vector<double> c(nx, 0.0);
    vector<double> c_new(nx, 0.0);
 
    // I.C.
    for (int i = 0; i < nx; i++)
    {
        if (i < nx / 4)
            c[i] = 0.0;
        else if (i < 3 * nx / 4)
            c[i] = 1.0;
        else
            c[i] = 0.0;
    }
 
    for (int t = 0; t <= total_timesteps; t++)
    {
        if (t > 0)
        {
            for (int i = 0; i < nx; i++)
            {
                // B.C.
                int im1 = (i == 0) ? nx - 1 : i - 1;
                int ip1 = (i == nx - 1) ? 0 : i + 1;
                // D.E
                c_new[i] = c[i] + (D * dt / (dx * dx)) * (c[im1] - 2.0 * c[i] + c[ip1]);
            }
            c.swap(c_new);
        }
 
        double time = t * dt;
        for (int k = 0; k < 4; k++) {
             if (t == (int)(save_times[k] / dt + 0.5)) {
                 ostringstream folder;
                  folder << "time_" << save_times[k]; 
                  mkdir(folder.str().c_str(),0777); 
                  ostringstream fname; 
                  fname << folder.str() << "/c_time_" << save_times[k] << ".dat"; 
                  ofstream out(fname.str().c_str()); 
                  out << "x c\n"; 
                  for (int i = 0; i < nx; i++) out << i * dx << " " << c[i] << "\n"; 
                  out.close(); cout << "saved time = " << time << "\n";
            }
        }
    }
 
    return 0;
}