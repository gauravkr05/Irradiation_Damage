// ============================================================
//  2D Allen-Cahn  --  EXPLICIT FINITE DIFFERENCE
//
//  Converted from the FFTW semi-implicit spectral version.
//
//      d(eta)/dt = -L ( g(eta) - 2*kappa*grad^2 eta )
//
//      f(eta) = A eta^2 (1-eta)^2
//      g(eta) = df/d(eta) = 2*A*eta*(1-eta)*(1-2*eta)
//
//  The spectral code treated -2*kappa*grad^2 implicitly (that is where
//  the  1 + 2*dt*L*kappa*k^2  denominator comes from) and g(eta)
//  explicitly. Here EVERYTHING is explicit, so:
//
//      eta_new = eta + dt * L * ( 2*kappa*grad^2 eta  -  g(eta) )
//
//  grad^2 is the 5-point central difference, periodic in x and y
//  (the FFT version is periodic, so this matches it).
//
//  I.C. : circular grain of radius 25 at eta = 1, matrix at eta = 0.2
//         (same as setInitialConditions in the spectral code)
//
//  Reads the SAME input.dat as the spectral version.
// ============================================================
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <sys/stat.h>
using namespace std;

int    NX, NY, NT;
double Lparam, kappa, A, dx, dy, dt;

// ---------------- Input ----------------
void get_input_parameters()
{
    FILE *fp = fopen("input.dat", "r");
    if (!fp) { printf("Error opening input.dat\n"); exit(EXIT_FAILURE); }
    if (fscanf(fp,
               "NX = %d\nNY = %d\nNT = %d\nL = %lf\nk = %lf\nA = %lf\ndx = %lf\ndy = %lf\ndt = %lf\n",
               &NX, &NY, &NT, &Lparam, &kappa, &A, &dx, &dy, &dt) != 9) {
        printf("Error reading input.dat\n");
        fclose(fp); exit(EXIT_FAILURE);
    }
    fclose(fp);
}

// ---------------- I.C. ----------------
void setInitialConditions(vector<double> &c)
{
    double centerX = (NX * dx) / 2.0;
    double centerY = (NY * dy) / 2.0;
    double radius  = 25.0;

    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
        {
            double x = i * dx, y = j * dy;
            double d = sqrt((x-centerX)*(x-centerX) + (y-centerY)*(y-centerY));
            c[i*NY + j] = (d <= radius) ? 1.0 : 0.2;
        }
}

// ---------------- df/d(eta) ----------------
static inline double g_of(double eta)
{
    return 2.0 * A * eta * (1.0 - eta) * (1.0 - 2.0*eta);
}

// ---------------- Output ----------------
void printConcentration(const char *folder, const vector<double> &c, double current_time)
{
    char filename[200];
    sprintf(filename, "%s/time%f.dat", folder, current_time);
    FILE *fp = fopen(filename, "w");
    if (!fp) { printf("Error opening %s\n", filename); return; }
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            fprintf(fp, "%.6f\t%.6f\t%.6f\n", i*dx, j*dy, c[i*NY + j]);
    fclose(fp);
}

void printInputParameters(const char *folder)
{
    char filename[200];
    sprintf(filename, "%s/inputs_taken.dat", folder);
    FILE *fp = fopen(filename, "w");
    if (!fp) return;
    fprintf(fp, "NX = %d\nNY = %d\nNT = %d\n", NX, NY, NT);
    fprintf(fp, "L = %lf\nk = %lf\nA = %lf\n", Lparam, kappa, A);
    fprintf(fp, "dx = %lf\ndy = %lf\ndt = %lf\n", dx, dy, dt);
    fclose(fp);
}

// ---------------- Main ----------------
int main()
{
    get_input_parameters();

    // ---- stability check (this is the price of going explicit) ----
    double r_diff = 2.0*kappa*Lparam*dt * (1.0/(dx*dx) + 1.0/(dy*dy));
    double r_reac = Lparam * dt * 2.0 * A;      // |dg/d(eta)| is largest at eta = 0 and 1
    printf("2*kappa*L*dt*(1/dx^2+1/dy^2) = %g   (needs <= 0.5)\n", r_diff);
    printf("L*dt*|g'|max                 = %g   (needs <= 1)\n", r_reac);
    if (r_diff > 0.5 || r_reac > 1.0) {
        printf("!! UNSTABLE : reduce dt in input.dat.  Try dt <= %g\n",
               0.5 / (2.0*kappa*Lparam*(1.0/(dx*dx) + 1.0/(dy*dy))));
        return 1;
    }

    vector<double> c    (NX*NY, 0.0);
    vector<double> c_new(NX*NY, 0.0);
    setInitialConditions(c);

    const char *folder = "output_ac_fd_2d";
    mkdir(folder, 0777);

    printf("Writing Initial conditions\n");
    printConcentration(folder, c, 0.0);
    printInputParameters(folder);

    int out_freq = NT / 10;
    if (out_freq < 1) out_freq = 1;

    for (int t = 1; t <= NT; t++)
    {
        for (int i = 0; i < NX; i++)
        {
            int im1 = (i - 1 + NX) % NX;      
            int ip1 = (i + 1) % NX;
            for (int j = 0; j < NY; j++)
            {
                int jm1 = (j - 1 + NY) % NY;
                int jp1 = (j + 1) % NY;

                double e = c[i*NY + j];

                double lap = (c[ip1*NY + j] - 2.0*e + c[im1*NY + j]) / (dx*dx)
                           + (c[i*NY + jp1] - 2.0*e + c[i*NY + jm1]) / (dy*dy);

                c_new[i*NY + j] = e + dt * Lparam * (2.0*kappa*lap - g_of(e));
            }
        }
        c.swap(c_new);

        if (t % out_freq == 0 || t == NT)
        {
            double current_time = t * dt;
            printf("Time = %f\n", current_time);
            printConcentration(folder, c, current_time);
        }
    }
    return 0;
}
