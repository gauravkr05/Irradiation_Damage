/*
    single_bubble_gb.cpp

    3-D phase-field simulation of ONE gas bubble sitting on a flat
    grain boundary (GB). Same model as phasefield_gb_bubbles3d.cpp,
    reduced to a single bubble so that growth rate can be measured
    quickly and parameters can be tuned before going to many bubbles.

    ------------------------------------------------------------------
    EQUATIONS (unchanged from the multi-bubble code, NO p(1-c) term)
    ------------------------------------------------------------------

        F = integral [ A0 f(c,n1,n2) + j/2 |grad c|^2
                       + jGB/2 (|grad n1|^2 + |grad n2|^2) ] dV

        f = (c^4/4 - c^2/2) + (n1^4/4 - n1^2/2) + (n2^4/4 - n2^2/2)
            + bGB/2 n1^2 n2^2 + bS c^2 (n1^2 + n2^2)

        mu      = A0 [ c^3 - c + 2 bS c (n1^2+n2^2) ] - j lap(c)
        dc/dt   = div( M grad mu )                    (conserved)
        M       = Ms*16*c^2*(1-c^2) + Mgb*n1^2*n2^2 + Mb

        dn1/dt  = -LGB { A0 [n1^3 - n1 + bGB n1 n2^2 + 2 bS c^2 n1] - jGB lap(n1) }
        dn2/dt  = -LGB { A0 [n2^3 - n2 + bGB n2 n1^2 + 2 bS c^2 n2] - jGB lap(n2) }

    ------------------------------------------------------------------
    WHAT IS DIFFERENT FROM THE MULTI-BUBBLE CODE (numerics only)
    ------------------------------------------------------------------

    1. One bubble, centred on the GB, in a smaller box.
    2. All parameters are read from the command line, so you do NOT
       recompile to change them:
           ./single_bubble_gb c0=0.15 R0=6 Mb=0.1 dt=0.001 out=run1
    3. c is NOT clamped to [0,1] after each step (clampC=0 default).
       Clamping deletes or creates gas and breaks the conservation that
       Cahn-Hilliard guarantees. Total gas is written to the CSV so you
       can check it stays constant. (Clamp is still used inside the
       mobility, exactly as before, to keep M >= Mb.)
    4. Mobility is computed once per cell per step (it was computed six
       times per cell before). Same values, about 2-3x faster.
    5. The program prints the explicit stability limit for dt. The old
       dt = 5e-5 was ~40x below the limit, i.e. ~40x wasted steps.
    6. Optional: n1, n2 are set to ~0 inside the initial bubble
       (initN=1, default), which is what the free energy wants anyway
       (the 2 bS c^2 n term drives n -> 0 where c = 1). This only
       removes an initial transient, the equations are the same.
    7. Physical time uses the dimensionless BULK diffusivity
       D* = Mb A0 (3 c0^2 - 1 + 2 bS), so that changing Mb rescales
       physical time correctly instead of looking like "faster growth".

    ------------------------------------------------------------------
    OUTPUT  (in folder given by out=..., default "run")
    ------------------------------------------------------------------

        growth.csv     time series of bubble size and gas balance
        c_XXXXXX.raw   float32 c field, shape [NZ][NY][NX] (x fastest)
                       written every saveEvery steps (0 = start/end only)
        params.txt     the parameters used for this run

    Compile:
        g++ -O3 -march=native -fopenmp -std=c++17 single_bubble_gb.cpp -o single_bubble_gb
*/

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// =====================================================================
// Parameters (defaults; override any of them as key=value on the
// command line)
// =====================================================================
struct Params
{
    // grid
    int NX = 64, NY = 64, NZ = 64;
    double dx = 1.0;
    double dxPhysical_m = 10.0e-9;     // 1 grid cell = 10 nm

    // bubble
    double R0 = 6.0;                   // initial radius in cells (60 nm)
    double width = 1.5;                // tanh width of initial profile

    // gas
    double c0 = 0.12;                  // matrix (supersaturated) concentration
    double cBubble = 1.0;
    double noise = 0.001;              // only at t = 0
    unsigned seed = 2024;

    // free energy
    double A0 = 4.0;
    int theta = 60;                    // 40, 50, 60, 80 -> bS, bGB below

    // kinetics
    double LGB = 0.10;
    double Mb = 0.01;
    double Mgb = 2.0;
    double Ms = 0.50;

    // time stepping
    double dt = -1.0;                  // <= 0 means: use 0.5 * stability limit
    long steps = 20000;
    long outputEvery = 200;
    long saveEvery = 0;                // 0 = save c only at start and end

    // misc
    int clampC = 0;                    // 1 = old behaviour (not conservative)
    int initN = 1;                     // 1 = n1,n2 -> 0 inside initial bubble
    double DPhysical = 3.4e-18;        // m^2/s, from the paper
    std::string out = "run";

    // derived
    double bS = 1.0, bGB = 1.0, j = 2.0, jGB = 2.0;

    void derive()
    {
        switch (theta)
        {
            case 40: bS = 1.0; bGB = 1.6; break;
            case 50: bS = 1.0; bGB = 1.3; break;
            case 60: bS = 1.0; bGB = 1.0; break;
            case 80: bS = 1.7; bGB = 0.7; break;
            default:
                std::cerr << "theta must be 40, 50, 60 or 80\n";
                std::exit(1);
        }
        j = 2.0 * bS;     // j/bS = 2 as in the paper
        jGB = 2.0 * bGB;
    }
};

static void parseArgs(int argc, char** argv, Params& P)
{
    for (int a = 1; a < argc; ++a)
    {
        std::string s(argv[a]);
        auto eq = s.find('=');
        if (eq == std::string::npos)
        {
            std::cerr << "Ignoring argument '" << s << "' (use key=value)\n";
            continue;
        }
        std::string k = s.substr(0, eq), v = s.substr(eq + 1);

        if      (k == "NX") P.NX = std::stoi(v);
        else if (k == "NY") P.NY = std::stoi(v);
        else if (k == "NZ") P.NZ = std::stoi(v);
        else if (k == "N")  { P.NX = P.NY = P.NZ = std::stoi(v); }
        else if (k == "R0") P.R0 = std::stod(v);
        else if (k == "width") P.width = std::stod(v);
        else if (k == "c0") P.c0 = std::stod(v);
        else if (k == "noise") P.noise = std::stod(v);
        else if (k == "seed") P.seed = (unsigned)std::stoul(v);
        else if (k == "A0") P.A0 = std::stod(v);
        else if (k == "theta") P.theta = std::stoi(v);
        else if (k == "LGB") P.LGB = std::stod(v);
        else if (k == "Mb") P.Mb = std::stod(v);
        else if (k == "Mgb") P.Mgb = std::stod(v);
        else if (k == "Ms") P.Ms = std::stod(v);
        else if (k == "dt") P.dt = std::stod(v);
        else if (k == "steps") P.steps = std::stol(v);
        else if (k == "outputEvery") P.outputEvery = std::stol(v);
        else if (k == "saveEvery") P.saveEvery = std::stol(v);
        else if (k == "clampC") P.clampC = std::stoi(v);
        else if (k == "initN") P.initN = std::stoi(v);
        else if (k == "out") P.out = v;
        else std::cerr << "Unknown parameter '" << k << "'\n";
    }
    P.derive();
}

// =====================================================================
// Simulation
// =====================================================================
class SingleBubble
{
public:
    explicit SingleBubble(const Params& p) : P(p)
    {
        NX = P.NX; NY = P.NY; NZ = P.NZ;
        N = (size_t)NX * NY * NZ;
        zGB = 0.5 * NZ - 0.5;   // GB plane between the two middle layers
        xC = 0.5 * NX;          // bubble centre
        yC = 0.5 * NY;

        for (auto* v : {&c, &n1, &n2, &cNew, &n1New, &n2New, &mu,
                        &lapC, &lapN1, &lapN2, &mob, &fx, &fy, &fz})
            v->assign(N, 0.0);
    }

    // ------------------------------------------------------------
    // Explicit stability estimate (largest eigenvalue of the
    // linearised operators on the 7-point stencil)
    // ------------------------------------------------------------
    double stabilityLimit() const
    {
        double k2max = 12.0 / (P.dx * P.dx);              // 4*d/dx^2, d=3
        double Mmax  = 4.0 * P.Ms + P.Mgb / 16.0 + P.Mb;  // max of M
        double fppC  = P.A0 * (2.0 + 2.0 * P.bS);         // max d2f/dc2
        double lamCH = Mmax * (P.j * k2max * k2max + fppC * k2max);
        double fppN  = P.A0 * (2.0 + P.bGB + 2.0 * P.bS);
        double lamAC = P.LGB * (P.jGB * k2max + fppN);
        return std::min(2.0 / lamCH, 2.0 / lamAC);
    }

    double bulkDiffusivity() const   // dimensionless D* in the matrix
    {
        return P.Mb * P.A0 * (3.0 * P.c0 * P.c0 - 1.0 + 2.0 * P.bS);
    }

    void run(double dt)
    {
        namespace fs = std::filesystem;
        fs::create_directories(P.out);

        initialize();
        writeParams(dt);

        std::ofstream csv(P.out + "/growth.csv");
        csv << "step,t,t_phys_s,V,R_eq,A_gb,a_gb,h_up,h_down,"
               "theta_est_deg,c_far,total_gas,gas_drift_rel,c_min,c_max\n";

        gas0 = totalGas();
        writeRow(csv, 0, dt);
        if (P.saveEvery >= 0) saveC(0);

        auto t0 = std::chrono::steady_clock::now();
        for (long s = 1; s <= P.steps; ++s)
        {
            step(dt);

            if (s % P.outputEvery == 0 || s == P.steps)
            {
                bool ok = writeRow(csv, s, dt);
                double el = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t0).count();
                std::cout << "step " << s << "/" << P.steps
                          << "  R_eq=" << std::fixed << std::setprecision(3) << lastReq
                          << "  c_far=" << std::setprecision(4) << lastCfar
                          << "  gas drift=" << std::scientific << std::setprecision(2) << lastDrift
                          << std::defaultfloat
                          << "  (" << std::setprecision(3) << el << " s)\n";
                if (!ok)
                {
                    std::cerr << "NaN/Inf detected: dt too large. Stopping.\n";
                    break;
                }
            }
            if (P.saveEvery > 0 && s % P.saveEvery == 0) saveC(s);
        }
        saveC(P.steps);
        std::cout << "Done. Output in " << P.out << "/\n";
    }

private:
    Params P;
    int NX, NY, NZ;
    size_t N;
    double zGB, xC, yC;
    double gas0 = 0.0, lastReq = 0.0, lastCfar = 0.0, lastDrift = 0.0;

    std::vector<double> c, n1, n2, cNew, n1New, n2New, mu,
                        lapC, lapN1, lapN2, mob, fx, fy, fz;

    inline size_t id(int x, int y, int z) const
    {
        return ((size_t)z * NY + y) * NX + x;
    }
    inline int px(int x) const { return x < 0 ? x + NX : (x >= NX ? x - NX : x); }
    inline int py(int y) const { return y < 0 ? y + NY : (y >= NY ? y - NY : y); }
    inline int cz(int z) const { return z < 0 ? 0 : (z >= NZ ? NZ - 1 : z); }

    // periodic x,y ; zero-gradient (no-flux) z  -- same as original
    inline double lap(const std::vector<double>& f, int x, int y, int z) const
    {
        double ctr = f[id(x, y, z)];
        return (f[id(px(x - 1), y, z)] + f[id(px(x + 1), y, z)]
              + f[id(x, py(y - 1), z)] + f[id(x, py(y + 1), z)]
              + f[id(x, y, cz(z - 1))] + f[id(x, y, cz(z + 1))]
              - 6.0 * ctr) / (P.dx * P.dx);
    }

    // ------------------------------------------------------------
    // Initial condition
    // ------------------------------------------------------------
    void initialize()
    {
        std::mt19937 gen(P.seed);
        std::uniform_real_distribution<double> nz(-P.noise, P.noise);

        for (int z = 0; z < NZ; ++z)
            for (int y = 0; y < NY; ++y)
                for (int x = 0; x < NX; ++x)
                {
                    size_t p = id(x, y, z);
                    double d = z - zGB;

                    double g1 = 0.5 * (1.0 - std::tanh(d / P.width));
                    double g2 = 0.5 * (1.0 + std::tanh(d / P.width));

                    double rx = x - xC, ry = y - yC, rz = d;
                    double r = std::sqrt(rx * rx + ry * ry + rz * rz);
                    double b = 0.5 * (1.0 - std::tanh((r - P.R0) / P.width));

                    double cv = P.c0 + (P.cBubble - P.c0) * b + nz(gen);
                    c[p] = std::clamp(cv, 0.0, 1.0);   // only at t = 0

                    if (P.initN)
                    {
                        g1 *= (1.0 - b);
                        g2 *= (1.0 - b);
                    }
                    n1[p] = g1;
                    n2[p] = g2;
                }
    }

    // ------------------------------------------------------------
    // One timestep: CH for c, AC for n1, n2  (explicit Euler)
    // ------------------------------------------------------------
    void step(double dt)
    {
        const double A0 = P.A0, bS = P.bS, bGB = P.bGB;
        const double j = P.j, jGB = P.jGB, dx = P.dx;

        // laplacians
        #pragma omp parallel for collapse(2) schedule(static)
        for (int z = 0; z < NZ; ++z)
            for (int y = 0; y < NY; ++y)
                for (int x = 0; x < NX; ++x)
                {
                    size_t p = id(x, y, z);
                    lapC[p]  = lap(c, x, y, z);
                    lapN1[p] = lap(n1, x, y, z);
                    lapN2[p] = lap(n2, x, y, z);
                }

        // chemical potential and mobility (once per cell)
        #pragma omp parallel for schedule(static)
        for (long long q = 0; q < (long long)N; ++q)
        {
            size_t p = (size_t)q;
            double cv = c[p];
            double g = n1[p] * n1[p] + n2[p] * n2[p];
            mu[p] = A0 * (cv * cv * cv - cv + 2.0 * bS * cv * g) - j * lapC[p];

            double cm = std::clamp(cv, 0.0, 1.0);
            mob[p] = P.Ms * 16.0 * cm * cm * (1.0 - cm * cm)
                   + P.Mgb * n1[p] * n1[p] * n2[p] * n2[p]
                   + P.Mb;
        }

        // face fluxes  F = M_face * grad(mu)   (dc/dt = div F)
        #pragma omp parallel for collapse(2) schedule(static)
        for (int z = 0; z < NZ; ++z)
            for (int y = 0; y < NY; ++y)
                for (int x = 0; x < NX; ++x)
                {
                    size_t p = id(x, y, z);
                    size_t qx = id(px(x + 1), y, z);
                    size_t qy = id(x, py(y + 1), z);
                    fx[p] = 0.5 * (mob[p] + mob[qx]) * (mu[qx] - mu[p]) / dx;
                    fy[p] = 0.5 * (mob[p] + mob[qy]) * (mu[qy] - mu[p]) / dx;
                    if (z == NZ - 1) fz[p] = 0.0;      // no flux through top
                    else
                    {
                        size_t qz = id(x, y, z + 1);
                        fz[p] = 0.5 * (mob[p] + mob[qz]) * (mu[qz] - mu[p]) / dx;
                    }
                }

        // update
        #pragma omp parallel for collapse(2) schedule(static)
        for (int z = 0; z < NZ; ++z)
            for (int y = 0; y < NY; ++y)
                for (int x = 0; x < NX; ++x)
                {
                    size_t p = id(x, y, z);
                    double fzm = (z > 0) ? fz[id(x, y, z - 1)] : 0.0;  // no flux bottom
                    double div = (fx[p] - fx[id(px(x - 1), y, z)]) / dx
                               + (fy[p] - fy[id(x, py(y - 1), z)]) / dx
                               + (fz[p] - fzm) / dx;

                    double cn = c[p] + dt * div;
                    if (P.clampC) cn = std::clamp(cn, 0.0, 1.0);
                    cNew[p] = cn;

                    double a = n1[p], b = n2[p], cc = c[p] * c[p];
                    double dF1 = A0 * (a * a * a - a + bGB * a * b * b + 2.0 * bS * cc * a)
                               - jGB * lapN1[p];
                    double dF2 = A0 * (b * b * b - b + bGB * b * a * a + 2.0 * bS * cc * b)
                               - jGB * lapN2[p];

                    n1New[p] = std::clamp(a - dt * P.LGB * dF1, 0.0, 1.0);
                    n2New[p] = std::clamp(b - dt * P.LGB * dF2, 0.0, 1.0);
                }

        c.swap(cNew);
        n1.swap(n1New);
        n2.swap(n2New);
    }

    // ------------------------------------------------------------
    // Diagnostics
    // ------------------------------------------------------------
    double totalGas() const
    {
        double s = 0.0;
        #pragma omp parallel for reduction(+ : s)
        for (long long q = 0; q < (long long)N; ++q) s += c[(size_t)q];
        return s;
    }

    bool writeRow(std::ofstream& csv, long s, double dt)
    {
        const double thr = 0.5;
        double V = 0.0, gas = 0.0, cmin = 1e9, cmax = -1e9;
        double farSum = 0.0, farCnt = 0.0;
        bool finite = true;

        // far-field matrix concentration: cells farther than
        // 0.4*min(NX,NY) from the bubble centre axis, away from the GB
        double rFar = 0.4 * std::min(NX, NY);

        for (int z = 0; z < NZ; ++z)
            for (int y = 0; y < NY; ++y)
                for (int x = 0; x < NX; ++x)
                {
                    double v = c[id(x, y, z)];
                    if (!std::isfinite(v)) finite = false;
                    gas += v;
                    cmin = std::min(cmin, v);
                    cmax = std::max(cmax, v);
                    if (v >= thr) V += 1.0;
                    double rx = x - xC, ry = y - yC;
                    if (std::sqrt(rx * rx + ry * ry) > rFar
                        && std::abs(z - zGB) > 4.0)
                    {
                        farSum += v;
                        farCnt += 1.0;
                    }
                }

        // projected area on the GB plane (column contains c >= 0.5)
        double A = 0.0;
        for (int y = 0; y < NY; ++y)
            for (int x = 0; x < NX; ++x)
                for (int z = 0; z < NZ; ++z)
                    if (c[id(x, y, z)] >= thr) { A += 1.0; break; }

        // half-heights along the column through the centre, measured from
        // the GB plane, with linear interpolation of the c = 0.5 crossing
        int xi = (int)std::floor(xC), yi = (int)std::floor(yC);
        double hUp = 0.0, hDown = 0.0;
        {
            int zs = (int)std::ceil(zGB);           // first layer above GB
            if (c[id(xi, yi, zs)] >= thr)
            {
                int z = zs;
                while (z + 1 < NZ && c[id(xi, yi, z + 1)] >= thr) ++z;
                double c1 = c[id(xi, yi, z)];
                double c2 = (z + 1 < NZ) ? c[id(xi, yi, z + 1)] : c1;
                double fr = (c1 != c2) ? (c1 - thr) / (c1 - c2) : 0.0;
                hUp = (z + fr) - zGB;
            }
            int zb = (int)std::floor(zGB);          // first layer below GB
            if (c[id(xi, yi, zb)] >= thr)
            {
                int z = zb;
                while (z - 1 >= 0 && c[id(xi, yi, z - 1)] >= thr) --z;
                double c1 = c[id(xi, yi, z)];
                double c2 = (z - 1 >= 0) ? c[id(xi, yi, z - 1)] : c1;
                double fr = (c1 != c2) ? (c1 - thr) / (c1 - c2) : 0.0;
                hDown = zGB - (z - fr);
            }
        }

        double Req = std::cbrt(3.0 * V / (4.0 * M_PI));
        double a = std::sqrt(A / M_PI);
        double h = 0.5 * (hUp + hDown);
        // symmetric lens: a = R sin(th), h = R (1 - cos(th))  ->  th = 2 atan(h/a)
        double th = (a > 0.0) ? 2.0 * std::atan(h / a) * 180.0 / M_PI : 0.0;
        double cfar = farCnt > 0 ? farSum / farCnt : 0.0;

        double t = s * dt;
        double tPhys = t * bulkDiffusivity() * P.dxPhysical_m * P.dxPhysical_m / P.DPhysical;
        double drift = (gas - gas0) / gas0;

        lastReq = Req; lastCfar = cfar; lastDrift = drift;

        csv << std::setprecision(10) << s << "," << t << "," << tPhys << ","
            << V << "," << Req << "," << A << "," << a << ","
            << hUp << "," << hDown << "," << th << ","
            << cfar << "," << gas << "," << drift << ","
            << cmin << "," << cmax << "\n";
        csv.flush();
        return finite;
    }

    void saveC(long s) const
    {
        std::ostringstream fn;
        fn << P.out << "/c_" << std::setw(6) << std::setfill('0') << s << ".raw";
        std::ofstream o(fn.str(), std::ios::binary);
        for (double v : c)
        {
            float f = (float)v;
            o.write(reinterpret_cast<const char*>(&f), sizeof(float));
        }
    }

    void writeParams(double dt) const
    {
        std::ofstream o(P.out + "/params.txt");
        o << "NX=" << NX << "\nNY=" << NY << "\nNZ=" << NZ
          << "\nR0=" << P.R0 << "\nc0=" << P.c0 << "\ntheta=" << P.theta
          << "\nbS=" << P.bS << "\nbGB=" << P.bGB << "\nj=" << P.j << "\njGB=" << P.jGB
          << "\nA0=" << P.A0 << "\nLGB=" << P.LGB
          << "\nMb=" << P.Mb << "\nMgb=" << P.Mgb << "\nMs=" << P.Ms
          << "\ndt=" << dt << "\nsteps=" << P.steps
          << "\nclampC=" << P.clampC << "\ninitN=" << P.initN
          << "\nD_bulk_dimensionless=" << bulkDiffusivity()
          << "\ndt_stability_limit=" << stabilityLimit() << "\n";
    }
};

// =====================================================================
int main(int argc, char** argv)
{
    Params P;
    parseArgs(argc, argv, P);

    SingleBubble sim(P);
    double dtMax = sim.stabilityLimit();
    double dt = (P.dt > 0.0) ? P.dt : 0.5 * dtMax;

    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif

    std::cout << "Single bubble on a GB\n"
              << "  grid           " << P.NX << " x " << P.NY << " x " << P.NZ << "\n"
              << "  R0             " << P.R0 << " cells (" << P.R0 * 10 << " nm)\n"
              << "  c0             " << P.c0 << "\n"
              << "  theta          " << P.theta << "  (bS=" << P.bS << ", bGB=" << P.bGB << ")\n"
              << "  Mb, Mgb, Ms    " << P.Mb << ", " << P.Mgb << ", " << P.Ms << "\n"
              << "  LGB            " << P.LGB << "\n"
              << "  D* bulk        " << sim.bulkDiffusivity() << "\n"
              << "  dt             " << dt << "   (stability limit ~ " << dtMax << ")\n"
              << "  steps          " << P.steps << "  -> t_end = " << dt * P.steps << "\n"
              << "  threads        " << threads << "\n\n";

    if (dt > dtMax)
        std::cout << "WARNING: dt is above the estimated stability limit.\n\n";

    sim.run(dt);
    return 0;
}
