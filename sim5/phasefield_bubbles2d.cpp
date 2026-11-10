// 2D phase-field simulation of gas bubbles ON a grain boundary plane.
//
// The domain IS the interface between the two grains (cubes), seen face-on.
// Two coupled fields live on that plane:
//   c   : conserved gas concentration          -> Cahn-Hilliard   (c = 1 bubble, c = 0 solid)
//   phi : non-conserved solid-contact order parameter -> Allen-Cahn (phi = 1 intact contact,
//                                                                    phi = 0 inside a bubble)
//
// F = int [ f(c,phi) + kap/2 |grad c|^2 + kap/2 |grad phi|^2 ] dA
// f = A0 [ c^4/4 - c^2/2 + phi^4/4 - phi^2/2 + bS c^2 phi^2 ]
//
// In the solid (phi = 1) the c-part becomes A0[c^4/4 + (bS - 1/2) c^2], which has its
// minimum at c = 0 provided 2 bS - 1 > 0. Inside a bubble (phi = 0) the minimum is c = 1.
//
//   dc/dt   = M lap( dF/dc ) + g        M = Dg / (A0 (2 bS - 1))    [Eq. 13 of the paper]
//   dphi/dt = -L dF/dphi
//   dF/dc   = A0 ( c^3 - c + 2 bS c phi^2 ) - kap lap(c)
//   dF/dphi = A0 ( phi^3 - phi + 2 bS c^2 phi ) - kap lap(phi)
//
// g = src * phi^2 * (1 - c) is a gas generation term acting only in the solid. It plays the
// role of fission-gas production: gas is fed into the lattice, diffuses to the bubbles, and
// they grow and coalesce, so coverage sweeps from small values up to 1.
//
// Initial state: many small circular nuclei at random positions on the plane, plus noise.
// Periodic boundary conditions in both directions.
//
// Measured every outEvery steps (bubble = connected set of cells with c > cthr):
//   coverage        = bubble area / domain area          (= A_total / A_domain in 2D)
//   A_largest/A_tot = area of the biggest bubble / total bubble area
//
// Output: bubble_data.dat with columns
//   step  time  coverage  Alargest_over_Atotal  Atotal  Nbubbles
// Optional field snapshots (DUMP=n) as snap_XXXXXX.bin for plotting the morphology.
//
// Build:  g++ -O3 -fopenmp -std=c++17 phasefield_bubbles2d.cpp -o pf2
// Run:    ./pf2 [Nx Ny nsteps out_every]
// Parameters can be overridden from the environment, e.g.
//         NSEED=800 R0=1.6 SRC=0.001 ./pf2 256 256 200000 200

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

using std::vector;

// ---------------- model parameters ----------------
static double A0   = 4.0;     // scaling coefficient of the bulk free energy
static double kap  = 2.0;     // gradient coefficient
static double bS   = 0.8;     // solid-bubble barrier coefficient (needs 2 bS - 1 > 0)
static double Dg   = 1.0;     // gas diffusivity in the solid
static double Lphi = 1.0;     // order-parameter mobility
static double dt   = 0.005;   // time step
static double c_bulk = 0.15;  // initial gas concentration dissolved in the solid
static double src  = 8.0e-4;  // gas generation rate in the solid
static double noise = 0.02;   // amplitude of random noise in the initial c field
static double cthr = 0.5;     // threshold that defines a bubble cell

int Nx = 256, Ny = 256;
inline int idx(int i, int j) { return j * Nx + i; }
inline int wx(int i) { return (i + Nx) % Nx; }
inline int wy(int j) { return (j + Ny) % Ny; }

// 5-point Laplacian, dx = 1, periodic both directions
inline double lap(const vector<double>& f, int i, int j) {
    return f[idx(wx(i + 1), j)] + f[idx(wx(i - 1), j)]
         + f[idx(i, wy(j + 1))] + f[idx(i, wy(j - 1))]
         - 4.0 * f[idx(i, j)];
}

// ---------------- bubble analysis ----------------
struct Stats { double coverage, ratio, atot; int nb; };

Stats analyse(const vector<double>& c) {
    const int N = Nx * Ny;
    vector<int> label(N, -1);
    vector<long> sizes;
    vector<int> stack;
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++) {
            int id = idx(i, j);
            if (c[id] <= cthr || label[id] >= 0) continue;
            int lab = (int)sizes.size();
            sizes.push_back(0);
            stack.clear();
            stack.push_back(id);
            label[id] = lab;
            while (!stack.empty()) {
                int cur = stack.back(); stack.pop_back();
                sizes[lab]++;
                int ci = cur % Nx, cj = cur / Nx;
                int nbr[4] = { idx(wx(ci + 1), cj), idx(wx(ci - 1), cj),
                               idx(ci, wy(cj + 1)), idx(ci, wy(cj - 1)) };
                for (int n = 0; n < 4; n++) {
                    int q = nbr[n];
                    if (c[q] <= cthr || label[q] >= 0) continue;
                    label[q] = lab;
                    stack.push_back(q);
                }
            }
        }
    long atot = 0;
    for (long s : sizes) atot += s;
    long amax = sizes.empty() ? 0 : *std::max_element(sizes.begin(), sizes.end());

    Stats s;
    s.coverage = (double)atot / (double)N;          // in 2D coverage IS the area fraction
    s.ratio    = atot > 0 ? (double)amax / (double)atot : 0.0;
    s.atot     = (double)atot;
    s.nb       = (int)sizes.size();
    return s;
}

int main(int argc, char** argv) {
    long nsteps = 600000, outEvery = 200;
    if (argc > 2) { Nx = atoi(argv[1]); Ny = atoi(argv[2]); }
    if (argc > 3) nsteps = atol(argv[3]);
    if (argc > 4) outEvery = atol(argv[4]);

    auto envd = [](const char* n, double& v) { if (const char* s = getenv(n)) v = atof(s); };
    envd("A0", A0); envd("KAP", kap); envd("BS", bS); envd("DG", Dg); envd("L", Lphi);
    envd("DT", dt); envd("CBULK", c_bulk); envd("SRC", src); envd("NOISE", noise);
    envd("CTHR", cthr);

    int nseed = 800; double r0 = 1.6; unsigned seed = 12345;
    if (const char* s = getenv("NSEED")) nseed = atoi(s);
    if (const char* s = getenv("SEED"))  seed  = (unsigned)atoi(s);
    envd("R0", r0);
    long dumpEvery = 0;
    if (const char* s = getenv("DUMP")) dumpEvery = atol(s);

    const int N = Nx * Ny;
    vector<double> c(N), phi(N), mu(N), cn(N), phin(N);

    // ---- initial condition: random small nuclei on the boundary plane ----
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U01(0.0, 1.0), U(-1.0, 1.0);
    vector<double> sx(nseed), sy(nseed), sr(nseed);
    for (int s = 0; s < nseed; s++) {
        sx[s] = U01(rng) * Nx;
        sy[s] = U01(rng) * Ny;
        sr[s] = r0 * (0.75 + 0.5 * U01(rng));
    }
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++) {
            double b = 0.0;                       // bubble indicator profile
            for (int s = 0; s < nseed; s++) {
                double dx = std::fabs(i - sx[s]); dx = std::min(dx, Nx - dx);
                double dy = std::fabs(j - sy[s]); dy = std::min(dy, Ny - dy);
                double r = std::sqrt(dx * dx + dy * dy);
                b = std::max(b, 0.5 * (1.0 - std::tanh((r - sr[s]) / 1.2)));
            }
            int id = idx(i, j);
            double cc = c_bulk + b * (1.0 - c_bulk) + noise * U(rng);
            c[id]   = std::min(1.0, std::max(0.0, cc));
            phi[id] = std::min(1.0, std::max(0.0, 1.0 - b));
        }

    FILE* fo = fopen("bubble_data.dat", "w");
    fprintf(fo, "# step time coverage Alargest_over_Atotal Atotal Nbubbles\n");

    const double M = Dg / (A0 * (2.0 * bS - 1.0));

    for (long step = 0; step <= nsteps; step++) {
        if (step % outEvery == 0) {
            Stats s = analyse(c);
            fprintf(fo, "%ld %.4f %.6f %.6f %.0f %d\n",
                    step, step * dt, s.coverage, s.ratio, s.atot, s.nb);
            fflush(fo);
            if ((step / outEvery) % 20 == 0)
                printf("step %8ld  cov %.3f  A1/Atot %.3f  nb %d\n",
                       step, s.coverage, s.ratio, s.nb);
            if (s.coverage > 0.995) {   // boundary fully covered, nothing left to follow
                printf("step %8ld  coverage saturated, stopping\n", step);
                break;
            }
        }
        if (dumpEvery > 0 && step % dumpEvery == 0) {
            char fn[128];
            snprintf(fn, sizeof fn, "snap_%06ld.bin", step);
            FILE* fs = fopen(fn, "wb");
            int hdr[2] = { Nx, Ny };
            fwrite(hdr, sizeof(int), 2, fs);
            vector<float> buf(N);
            for (int n = 0; n < N; n++) buf[n] = (float)c[n];
            fwrite(buf.data(), sizeof(float), N, fs);
            fclose(fs);
        }
        if (step == nsteps) break;

        // chemical potential  dF/dc
        #pragma omp parallel for
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++) {
                int id = idx(i, j);
                double cc = c[id], p = phi[id];
                mu[id] = A0 * (cc * cc * cc - cc + 2.0 * bS * cc * p * p) - kap * lap(c, i, j);
            }

        // Cahn-Hilliard for c, Allen-Cahn for phi, explicit Euler
        #pragma omp parallel for
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++) {
                int id = idx(i, j);
                double cc = c[id], p = phi[id];
                double dFdp = A0 * (p * p * p - p + 2.0 * bS * cc * cc * p) - kap * lap(phi, i, j);
                double gen = src * p * p * (1.0 - cc);     // gas generated in the solid only
                if (gen < 0.0) gen = 0.0;
                cn[id]   = cc + dt * (M * lap(mu, i, j) + gen);
                phin[id] = p  - dt * Lphi * dFdp;
            }
        c.swap(cn); phi.swap(phin);
    }
    fclose(fo);
    printf("done. wrote bubble_data.dat\n");
    return 0;
}
