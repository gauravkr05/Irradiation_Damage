#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <sys/stat.h>
using namespace std;

void run(double dt, int steps)
{
    double D  = 1.0;
    double dx = 1.0;
    int    nx = 100;
    double r  = D * dt / (dx * dx);

    vector<double> c(nx, 0.0);
    vector<double> c_new(nx, 0.0);

    // I.C.
    for (int i = 0; i < nx; i++)
    {
        if (i < nx / 4)       c[i] = 0.0;
        else if (i < 3*nx/4)  c[i] = 1.0;
        else                   c[i] = 0.0;
    }

    for (int t = 1; t <= steps; t++)
    {
        for (int i = 0; i < nx; i++)
        {
            // B.C.
            int im1 = (i == 0)      ? nx-1 : i-1;
            int ip1 = (i == nx-1)   ? 0    : i+1;
            // D.E
            c_new[i] = c[i] + r * (c[im1] - 2.0*c[i] + c[ip1]);
        }
        c.swap(c_new);
    }

    // save result
    ostringstream folder;
    folder << "stability_dt" << dt;
    mkdir(folder.str().c_str(),0777);

    ostringstream fname;
    fname << folder.str() << "/c.dat";
    ofstream out(fname.str().c_str());
    out << "x c\n";
    for (int i = 0; i < nx; i++)
        out << i*dx << " " << c[i] << "\n";
    out.close();

    cout << "dt=" << dt << "  r=" << r << "  steps=" << steps
         << "  physical_time=" << steps*dt << "\n";
}

int main()
{
    // all run to same physical time = 10
    run(0.25, 40);   // r=0.25  stable
    run(0.40, 25);   // r=0.40  stable but close
    run(0.50, 20);   // r=0.50  borderline
    run(0.60, 17);   // r=0.60  unstable
    run(1.00, 10);   // r=1.00  very unstable

    return 0;
}
