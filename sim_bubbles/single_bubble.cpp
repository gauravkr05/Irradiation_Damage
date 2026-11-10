// 2D phase-field simulation of a SINGLE gas bubble growing in a matrix that is
// supersaturated with point defects / gas atoms.
//
// Same free energy as the grain-boundary model, with one grain only:
//   c   : conserved gas / point-defect concentration  -> Cahn-Hilliard  (c = 1 bubble, c = 0 matrix)
//   phi : non-conserved solid order parameter         -> Allen-Cahn     (phi = 1 matrix, 0 bubble)
//
// F = int [ f(c,phi) + kap/2 |grad c|^2 + kap/2 |grad phi|^2 ] dA
// f = A0 [ c^4/4 - c^2/2 + phi^4/4 - phi^2/2 + bS c^2 phi^2 ]
//
//   dc/dt   = M lap( dF/dc ),   M = Dg / (A0 (2 bS - 1))      [Eq. 13]
//   dphi/dt = -L dF/dphi
//
// In the matrix (phi = 1) the c-part is A0[c^4/4 + (bS - 1/2)c^2], minimum at c = 0, so a
// matrix held at c = c_sup > 0 is supersaturated. That excess diffuses to the bubble and the
// bubble grows. Inside the bubble (phi = 0) the minimum is c = 1.
//
// Initial state: one circular bubble of radius R0 at the centre of a uniformly
// supersaturated matrix. Periodic boundaries.
//
// Outputs
//   bubble_radius.dat : step, t_sim, t_hours, area, radius  (radius from area, R = sqrt(A/pi))
//   snap_XXXXXX.bin   : raw float32 c field, for the panel figure
//   meta.txt          : scaling information used to convert to real time
//
// Real-time scaling. With grid spacing dx (metres) and real diffusivity D_real (m^2/s),
// matching the sim equation to Fick's law gives   t_real = t_sim * Dg * dx^2 / D_real.
// The defaults (dx = 7.5 nm, D_real = 3.4e-18 m^2/s, Dg = 2) put the last frame near 11 h.
//
// Build:  g++ -O3 -fopenmp -std=c++17 single_bubble.cpp -o sb
// Run:    ./sb [Nx Ny nsteps dump_every]
// Override anything from the environment, e.g.  CSUP=0.2 R0=10 DXNM=5 ./sb 160 160 240000 24000

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

using std::vector;

// ---------------- model parameters ----------------
static double A0    = 4.0;      // scaling coefficient of the bulk free energy
static double kap   = 2.0;      // gradient coefficient
static double bS    = 0.8;      // bubble surface barrier coefficient (needs 2 bS - 1 > 0)
static double Dg    = 2.0;      // gas / defect diffusivity in the matrix (sim units)
static double Lphi  = 1.0;      // order-parameter mobility
static double dt    = 0.005;    // time step
static double c_sup = 0.15;     // supersaturation held in the matrix at t = 0
static double R0    = 12.0;     // initial bubble radius, in cells
static double cthr  = 0.5;      // threshold defining the bubble

// ---------------- real-unit scaling ----------------
static double dx_nm  = 7.5;     // grid spacing in nm
static double D_real = 3.4e-18; // real diffusivity in m^2/s

int Nx = 160, Ny = 160;
inline int idx(int i, int j) { return j * Nx + i; }
inline int wx(int i) { return (i + Nx) % Nx; }
inline int wy(int j) { return (j + Ny) % Ny; }

inline double lap(const vector<double>& f, int i, int j) {
    return f[idx(wx(i + 1), j)] + f[idx(wx(i - 1), j)]
         + f[idx(i, wy(j + 1))] + f[idx(i, wy(j - 1))]
         - 4.0 * f[idx(i, j)];
}

int main(int argc, char** argv) {
    long nsteps = 240000, dumpEvery = 24000;
    if (argc > 2) { Nx = atoi(argv[1]); Ny = atoi(argv[2]); }
    if (argc > 3) nsteps    = atol(argv[3]);
    if (argc > 4) dumpEvery = atol(argv[4]);

    auto envd = [](const char* n, double& v) { if (const char* s = getenv(n)) v = atof(s); };
    envd("A0", A0); envd("KAP", kap); envd("BS", bS); envd("DG", Dg); envd("L", Lphi);
    envd("DT", dt); envd("CSUP", c_sup); envd("R0", R0); envd("CTHR", cthr);
    envd("DXNM", dx_nm); envd("DREAL", D_real);

    const int N = Nx * Ny;
    vector<double> c(N), phi(N), mu(N), cn(N), phin(N);

    // seconds of real time per unit of simulation time
    const double dx_m  = dx_nm * 1e-9;
    const double tscale = Dg * dx_m * dx_m / D_real;

    // ---- initial condition: one circular bubble in a supersaturated matrix ----
    double cx = Nx * 0.5, cy = Ny * 0.5;
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++) {
            double r = std::sqrt((i - cx) * (i - cx) + (j - cy) * (j - cy));
            double b = 0.5 * (1.0 - std::tanh((r - R0) / 1.5));   // 1 inside bubble, 0 outside
            int id = idx(i, j);
            c[id]   = c_sup + b * (1.0 - c_sup);
            phi[id] = 1.0 - b;
        }

    FILE* fo = fopen("bubble_radius.dat", "w");
    fprintf(fo, "# step t_sim t_hours area radius\n");
    FILE* fm = fopen("meta.txt", "w");
    fprintf(fm, "dt %g\nDg %g\ndx_nm %g\nD_real %g\ntscale_s_per_tsim %g\n",
            dt, Dg, dx_nm, D_real, tscale);
    fclose(fm);

    const double M = Dg / (A0 * (2.0 * bS - 1.0));

    for (long step = 0; step <= nsteps; step++) {
        if (step % (dumpEvery / 10 > 0 ? dumpEvery / 10 : 1) == 0) {
            long area = 0;
            for (int n = 0; n < N; n++) if (c[n] > cthr) area++;
            double tsim = step * dt;
            double R = std::sqrt(area / M_PI);
            fprintf(fo, "%ld %.4f %.4f %ld %.4f\n", step, tsim, tsim * tscale / 3600.0, area, R);
            fflush(fo);
        }
        if (step % dumpEvery == 0) {
            char fn[128];
            snprintf(fn, sizeof fn, "snap_%06ld.bin", step);
            FILE* fs = fopen(fn, "wb");
            int hdr[2] = { Nx, Ny };
            fwrite(hdr, sizeof(int), 2, fs);
            vector<float> buf(N);
            for (int n = 0; n < N; n++) buf[n] = (float)c[n];
            fwrite(buf.data(), sizeof(float), N, fs);
            fclose(fs);
            long area = 0;
            for (int n = 0; n < N; n++) if (c[n] > cthr) area++;
            printf("step %8ld  t = %7.2f sim = %6.2f h   R = %6.2f cells = %6.1f nm\n",
                   step, step * dt, step * dt * tscale / 3600.0,
                   std::sqrt(area / M_PI), std::sqrt(area / M_PI) * dx_nm);
        }
        if (step == nsteps) break;

        #pragma omp parallel for
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++) {
                int id = idx(i, j);
                double cc = c[id], p = phi[id];
                mu[id] = A0 * (cc * cc * cc - cc + 2.0 * bS * cc * p * p) - kap * lap(c, i, j);
            }

        #pragma omp parallel for
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++) {
                int id = idx(i, j);
                double cc = c[id], p = phi[id];
                double dFdp = A0 * (p * p * p - p + 2.0 * bS * cc * cc * p) - kap * lap(phi, i, j);
                cn[id]   = cc + dt * M * lap(mu, i, j);
                phin[id] = p  - dt * Lphi * dFdp;
            }
        c.swap(cn); phi.swap(phin);
    }
    fclose(fo);
    printf("done. wrote bubble_radius.dat, meta.txt and the snapshots\n");
    return 0;
}
