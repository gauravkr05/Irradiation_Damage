/*
    phasefield_gb_bubbles3d.cpp   (conserved version)

    3-D phase-field simulation of gas bubble growth and coalescence on a
    flat grain boundary (GB) between two grains.

    Domain   : 128 x 128 x 64, dx = 10 nm
    Grain 1  : z = 0 ... 31      (n1 = 1)
    Grain 2  : z = 32 ... 63     (n2 = 1)
    GB plane : z = 31.5
    x, y periodic, z no-flux.

    Free energy:

        F = integral [ A0 f(c,n1,n2) + j/2 |grad c|^2
                       + jGB/2 (|grad n1|^2 + |grad n2|^2) ] dV

        f = c^4/4 - c^2/2
          + n1^4/4 - n1^2/2 + n2^4/4 - n2^2/2
          + bGB/2 n1^2 n2^2
          + bS c^2 (n1^2 + n2^2)

    Gas (conserved Cahn-Hilliard with variable mobility, no source):

        dc/dt = div( M grad mu )

        mu = A0 [ c^3 - c + 2 bS c (n1^2 + n2^2) ] - j lap(c)

        M  = Ms * 16 c^2 (1 - c^2)        bubble surface
           + Mgb * 16 n1^2 n2^2           grain boundary
           + Mb                           bulk

    Grains (Allen-Cahn):

        dn_i/dt = -LGB dF/dn_i

    ------------------------------------------------------------------
    GAS CONSERVATION
    ------------------------------------------------------------------

    The total gas sum(c) is conserved to round-off:

      * The update is in flux form. The flux through each cell face,
        0.5 (M_p + M_q) (mu_q - mu_p), is added to one cell and
        subtracted from its neighbour, so it cancels in the sum.
        x, y are periodic and the z faces carry zero flux, so no gas
        leaves or enters through the boundaries.

      * There is NO source term. The earlier P (1 - c) generation term
        has been removed completely.

      * c is NOT clamped to [0, 1] after each step. Clamping cut off the
        small negative undershoots that Cahn-Hilliard produces next to
        interfaces (and overshoots above 1), and every cut changed the
        total gas. In the previous version this added about 22% gas
        over a run. c may now dip slightly below 0 or above 1 near
        interfaces. That is normal for Cahn-Hilliard and is reported
        in the output (c_min, c_max). The mobility still uses c clamped
        to [0, 1] so that it stays positive; this does not affect
        conservation, because M only scales the fluxes.

    Every metrics row records the total gas and its relative drift from
    the initial value (column mass_drift), so conservation can be
    checked directly.

    ------------------------------------------------------------------
    TIME STEP (adaptive, from the stability limit)
    ------------------------------------------------------------------

    The update is explicit (forward Euler). For a Fourier mode with
    discrete Laplacian eigenvalue -k^2, k^2 <= 12 for the 7-point
    stencil in 3-D (dx = 1), the linearised Cahn-Hilliard rate is

        lambda = - M k^2 ( A0 f_cc + j k^2 ),   f_cc = 3c^2 - 1 + 2 bS (n1^2 + n2^2)

    and forward Euler is stable when dt * |lambda| <= 2. Taking the
    worst case k^2 = 12:

        dt_CH <= 2 / [ M_max ( 12 A0 max(f_cc) + 144 j ) ]

    M_max and max(f_cc) change as bubbles grow (the GB mobility after
    relaxation can exceed Mgb), so they are measured over the whole grid
    EVERY step and dt is recomputed. The Allen-Cahn limit

        dt_AC <= 2 / [ LGB ( A0 max(f_nn) + 12 jGB ) ]

    is applied as well. The step used is

        dt = min( dtMax, safety * dt_CH, safety * dt_AC )

    with safety = 0.5 by default. The run is controlled by end TIME,
    not by a number of steps, and steps are shortened to land exactly on
    the output times, so runs at different contact angles (which have
    different stable dt) are sampled at the same times.

    ------------------------------------------------------------------
    PARAMETERS FROM THE COMMAND LINE (all optional), e.g.

        ./phasefield_gb_bubbles3d theta=60 tEnd=20

        theta    contact angle: 40, 50, 60 or 80          (default 60)
        tEnd     end time, dimensionless                  (default 20)
        outDt    write a metrics row every outDt          (default 0.1)
        snapDt   write a 3-D morphology file every snapDt (default 3)
        dtMax    largest allowed time step                (default 1e-3)
        safety   fraction of the stability limit used     (default 0.5)
        minThick min. bubble thickness (cells) for a GB pixel to count
                 as covered                               (default 3 = 30 nm)
        Mgb      GB mobility coefficient                  (default 0.5)
        Ms       surface mobility                         (default 0.05)
        Mb       bulk mobility                            (default 5)
        LGB      Allen-Cahn mobility                      (default 1)
        cInit    initial gas in the matrix                (default 0.30)
        r0       initial bubble radius, cells             (default 3.8)
        nBubbles number of initial bubbles                (default 33)
        minSep   min. centre spacing of bubbles, cells    (default 19)
        seed     random seed                              (default 11)

    Output:
        bubble_metrics.csv
        morphology_step_N.raw   (float32 c, shape [NZ][NY][NX]), N = step number

    Compile:
        g++ -O3 -march=native -fopenmp -std=c++17 phasefield_gb_bubbles3d.cpp -o phasefield_gb_bubbles3d
*/

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// ------------------------------------------------------------
// Grid
// ------------------------------------------------------------
constexpr int NX = 128;
constexpr int NY = 128;
constexpr int NZ = 64;
constexpr int N  = NX * NY * NZ;
constexpr double ZGB = 31.5;

constexpr double DX_M = 10.0e-9;      // one grid cell = 10 nm
constexpr double D_PHYS = 3.4e-18;    // m^2/s, from the paper (time conversion only)

// ------------------------------------------------------------
// Parameters (defaults; can be overridden on the command line)
// ------------------------------------------------------------
struct Params
{
    int theta = 60;

    double A0 = 4.0;
    double bS = 1.0;
    double bGB = 1.0;
    double j = 2.0;     // = 2 bS
    double jGB = 2.0;   // = 2 bGB

    double Ms = 0.05;
    double Mgb = 0.5;
    double Mb = 5.0;
    double LGB = 1.0;

    double tEnd = 20.0;     // end time (dimensionless)
    double outDt = 0.1;     // metrics every this much time
    double snapDt = 3.0;    // 3-D morphology file every this much time
    double dtMax = 1.0e-3;  // upper cap on the time step
    double safety = 0.5;    // fraction of the stability limit actually used
    int minThick = 3;       // GB pixel counts as bubble if >= this many cells in its column have c >= 0.5

    int nBubbles = 33;
    double r0 = 3.8;          // 38 nm
    double minSep = 19.0;     // 190 nm (near the random-packing limit for 33 bubbles)
    double cInit = 0.30;     // above ~0.35 new bubbles nucleate on the GB
    double noise = 0.001;
    unsigned seed = 11;       // places all 33 bubbles at minSep = 19

    void setTheta(int t)
    {
        theta = t;
        if (t == 40)      { bS = 1.0; bGB = 1.6; }
        else if (t == 50) { bS = 1.0; bGB = 1.3; }
        else if (t == 60) { bS = 1.0; bGB = 1.0; }
        else if (t == 80) { bS = 1.7; bGB = 0.7; }
        else { std::cerr << "theta must be 40, 50, 60 or 80\n"; std::exit(1); }
        j = 2.0 * bS;
        jGB = 2.0 * bGB;
    }
};

static inline int idx(int x, int y, int z) { return (z * NY + y) * NX + x; }

class Simulation
{
public:
    explicit Simulation(const Params& p)
        : P(p), c(N), n1(N), n2(N), cN(N), n1N(N), n2N(N),
          mu(N), mob(N), lapN1(N), lapN2(N), rng(p.seed)
    {
        // neighbour tables: periodic x/y, mirrored (no-flux) z
        for (int i = 0; i < NX; i++) { xm[i] = (i + NX - 1) % NX; xp[i] = (i + 1) % NX; }
        for (int i = 0; i < NY; i++) { ym[i] = (i + NY - 1) % NY; yp[i] = (i + 1) % NY; }
        for (int i = 0; i < NZ; i++) { zm[i] = std::max(i - 1, 0); zp[i] = std::min(i + 1, NZ - 1); }
    }

    void initialize()
    {
        // ---------------- bubble centres on the GB ----------------
        std::uniform_real_distribution<double> ux(0.0, NX), uy(0.0, NY);
        std::vector<std::pair<double, double>> centres;
        int attempts = 0;
        while ((int)centres.size() < P.nBubbles && attempts++ < 1000000)
        {
            double x = ux(rng), y = uy(rng);
            bool ok = true;
            for (auto& q : centres)
            {
                double dx = std::abs(x - q.first);  dx = std::min(dx, NX - dx);
                double dy = std::abs(y - q.second); dy = std::min(dy, NY - dy);
                if (dx * dx + dy * dy < P.minSep * P.minSep) { ok = false; break; }
            }
            if (ok) centres.push_back({x, y});
        }
        std::cout << "Initial bubbles placed = " << centres.size() << "\n";

        std::uniform_real_distribution<double> un(-P.noise, P.noise);
        const double w = 1.5;   // interface half-width (cells)

        for (int z = 0; z < NZ; z++)
            for (int y = 0; y < NY; y++)
                for (int x = 0; x < NX; x++)
                {
                    int p = idx(x, y, z);
                    double dz = z - ZGB;

                    // bubble indicator b: 1 inside, 0 outside
                    double b = 0.0;
                    for (auto& q : centres)
                    {
                        double dx = x - q.first;
                        double dy = y - q.second;
                        if (dx >  NX / 2.0) dx -= NX;
                        if (dx < -NX / 2.0) dx += NX;
                        if (dy >  NY / 2.0) dy -= NY;
                        if (dy < -NY / 2.0) dy += NY;
                        double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                        b = std::max(b, 0.5 * (1.0 - std::tanh((r - P.r0) / w)));
                    }

                    // grains: tanh profile across the GB, removed inside bubbles
                    double g1 = 0.5 * (1.0 - std::tanh(dz / w));
                    double g2 = 0.5 * (1.0 + std::tanh(dz / w));
                    n1[p] = g1 * (1.0 - b);
                    n2[p] = g2 * (1.0 - b);

                    // gas (clamped once here only; this sets the initial inventory)
                    double v = P.cInit + (1.0 - P.cInit) * b + un(rng);
                    c[p] = std::clamp(v, 0.0, 1.0);
                }

        mass0 = totalGas();
    }

    double totalGas() const
    {
        double s = 0.0;
        #pragma omp parallel for reduction(+:s) schedule(static)
        for (int p = 0; p < N; p++) s += c[p];
        return s;
    }

    // one explicit time step; dtLimit = time left to the next output.
    // Returns the dt actually used.
    double step(double dtLimit)
    {
        const double A0 = P.A0, bS = P.bS, bGB = P.bGB, j = P.j, jGB = P.jGB;
        const double Ms = P.Ms, Mgb = P.Mgb, Mb = P.Mb;

        // 1. chemical potential, mobility, grain laplacians,
        //    and the grid maxima needed for the stability limit
        double Mmax = 0.0, fccMax = 0.0, fnnMax = 0.0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(max:Mmax, fccMax, fnnMax)
        for (int z = 0; z < NZ; z++)
            for (int y = 0; y < NY; y++)
                for (int x = 0; x < NX; x++)
                {
                    int p = idx(x, y, z);
                    int pxm = idx(xm[x], y, z), pxp = idx(xp[x], y, z);
                    int pym = idx(x, ym[y], z), pyp = idx(x, yp[y], z);
                    int pzm = idx(x, y, zm[z]), pzp = idx(x, y, zp[z]);

                    double cv = c[p], a = n1[p], b = n2[p];
                    double lapC = c[pxm] + c[pxp] + c[pym] + c[pyp] + c[pzm] + c[pzp] - 6.0 * cv;
                    lapN1[p] = n1[pxm] + n1[pxp] + n1[pym] + n1[pyp] + n1[pzm] + n1[pzp] - 6.0 * a;
                    lapN2[p] = n2[pxm] + n2[pxp] + n2[pym] + n2[pyp] + n2[pzm] + n2[pzp] - 6.0 * b;

                    mu[p] = A0 * (cv * cv * cv - cv + 2.0 * bS * cv * (a * a + b * b)) - j * lapC;

                    // mobility uses c clamped to [0,1] so M stays positive.
                    // This does not affect conservation: M only scales fluxes.
                    double cc = std::clamp(cv, 0.0, 1.0);
                    double m = Ms * 16.0 * cc * cc * (1.0 - cc * cc)
                             + Mgb * 16.0 * a * a * b * b
                             + Mb;
                    mob[p] = m;

                    // second derivatives of f, for the stability limit
                    double fcc = 3.0 * cv * cv - 1.0 + 2.0 * bS * (a * a + b * b);
                    double f11 = 3.0 * a * a - 1.0 + bGB * b * b + 2.0 * bS * cv * cv;
                    double f22 = 3.0 * b * b - 1.0 + bGB * a * a + 2.0 * bS * cv * cv;

                    Mmax = std::max(Mmax, m);
                    fccMax = std::max(fccMax, fcc);
                    fnnMax = std::max(fnnMax, std::max(f11, f22));
                }

        // 2. stable time step
        const double L = P.LGB;
        double dtCH = 2.0 / (Mmax * (12.0 * A0 * fccMax + 144.0 * j));
        double dtAC = 2.0 / (L * (A0 * fnnMax + 12.0 * jGB));
        if (!std::isfinite(dtCH) || !std::isfinite(dtAC) || !(dtCH > 0.0) || !(dtAC > 0.0))
        {
            std::cerr << "Fields became non-finite (numerical blow-up). Reduce safety or dtMax.\n";
            std::exit(2);
        }
        double dt = std::min({P.dtMax, P.safety * dtCH, P.safety * dtAC, dtLimit});

        // 3. update c (flux divergence, conservative, no clamp) and n1, n2
        #pragma omp parallel for collapse(2) schedule(static)
        for (int z = 0; z < NZ; z++)
            for (int y = 0; y < NY; y++)
                for (int x = 0; x < NX; x++)
                {
                    int p = idx(x, y, z);
                    int pxm = idx(xm[x], y, z), pxp = idx(xp[x], y, z);
                    int pym = idx(x, ym[y], z), pyp = idx(x, yp[y], z);

                    double m = mob[p], u = mu[p];
                    double div =
                          0.5 * (m + mob[pxp]) * (mu[pxp] - u) - 0.5 * (m + mob[pxm]) * (u - mu[pxm])
                        + 0.5 * (m + mob[pyp]) * (mu[pyp] - u) - 0.5 * (m + mob[pym]) * (u - mu[pym]);
                    // z faces: zero flux through the top and bottom walls
                    if (z < NZ - 1) { int q = idx(x, y, z + 1); div += 0.5 * (m + mob[q]) * (mu[q] - u); }
                    if (z > 0)      { int q = idx(x, y, z - 1); div -= 0.5 * (m + mob[q]) * (u - mu[q]); }

                    double cv = c[p], a = n1[p], b = n2[p];
                    cN[p] = cv + dt * div;        // no clamp: gas is conserved exactly

                    double dF1 = A0 * (a * a * a - a + bGB * a * b * b + 2.0 * bS * cv * cv * a) - jGB * lapN1[p];
                    double dF2 = A0 * (b * b * b - b + bGB * b * a * a + 2.0 * bS * cv * cv * b) - jGB * lapN2[p];
                    // grains are not conserved, so keeping them in [0,1] is harmless
                    n1N[p] = std::clamp(a - dt * L * dF1, 0.0, 1.0);
                    n2N[p] = std::clamp(b - dt * L * dF2, 0.0, 1.0);
                }

        c.swap(cN);
        n1.swap(n1N);
        n2.swap(n2N);
        return dt;
    }

    // projected bubble metrics on the GB plane
    //   coverage   = covered GB pixels / (NX*NY)
    //   A1/Atot    = largest connected covered region / all covered pixels
    //   porosity   = cells with c >= 0.5 / all cells
    //   mass_drift = (sum c - sum c at t=0) / (sum c at t=0)
    void writeMetrics(std::ofstream& f, long s, double t, double dt, bool print)
    {
        std::vector<unsigned char> mask(NX * NY, 0);
        long poreCells = 0;
        #pragma omp parallel for reduction(+:poreCells) collapse(2)
        for (int y = 0; y < NY; y++)
            for (int x = 0; x < NX; x++)
            {
                // A thin gas-rich film (1-2 cells, c ~ 0.5-0.6) forms along the GB
                // between bubbles. Counting it as bubble makes separate bubbles
                // flicker in and out of contact. Require a minimum thickness.
                int thick = 0;
                for (int z = 0; z < NZ; z++)
                    if (c[idx(x, y, z)] >= 0.5) thick++;
                poreCells += thick;
                if (thick >= P.minThick) mask[y * NX + x] = 1;
            }

        // connected components, 8-neighbour, periodic
        std::vector<unsigned char> seen(NX * NY, 0);
        std::vector<int> queue(NX * NY);
        long total = 0, largest = 0; int count = 0;
        for (int s0 = 0; s0 < NX * NY; s0++)
        {
            if (!mask[s0] || seen[s0]) continue;
            int head = 0, tail = 0;
            queue[tail++] = s0; seen[s0] = 1;
            long area = 0;
            while (head < tail)
            {
                int q = queue[head++]; area++;
                int qx = q % NX, qy = q / NX;
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++)
                    {
                        if (!dx && !dy) continue;
                        int nx = (qx + dx + NX) % NX, ny = (qy + dy + NY) % NY;
                        int r = ny * NX + nx;
                        if (mask[r] && !seen[r]) { seen[r] = 1; queue[tail++] = r; }
                    }
            }
            count++; total += area; largest = std::max(largest, area);
        }

        // gas inventory and range of c
        double mass = 0.0, cmin = 1e30, cmax = -1e30;
        #pragma omp parallel for reduction(+:mass) reduction(min:cmin) reduction(max:cmax) schedule(static)
        for (int p = 0; p < N; p++) { mass += c[p]; cmin = std::min(cmin, c[p]); cmax = std::max(cmax, c[p]); }
        double drift = (mass - mass0) / mass0;

        if (!std::isfinite(mass))
        {
            std::cerr << "Field became non-finite at step " << s << ", t = " << t
                      << ". Reduce safety or dtMax.\n";
            std::exit(2);
        }

        double tSec = t * DX_M * DX_M / D_PHYS;
        double coverage = double(total) / (NX * NY);
        double frac = total > 0 ? double(largest) / total : 0.0;
        double porosity = double(poreCells) / N;

        f << s << "," << std::setprecision(10) << t << "," << tSec << "," << tSec / 3600.0 << ","
          << count << "," << total << "," << largest << "," << frac << ","
          << coverage << "," << porosity << ","
          << dt << "," << mass << "," << drift << "," << cmin << "," << cmax << std::endl;

        if (print)
            std::cout << "t " << std::fixed << std::setprecision(2) << std::setw(6) << t
                      << "  step " << std::setw(7) << s
                      << "  dt " << std::scientific << std::setprecision(2) << dt
                      << "  bubbles " << std::setw(3) << count
                      << std::fixed << std::setprecision(3)
                      << "  coverage " << coverage
                      << "  A1/Atot " << frac
                      << "  porosity " << std::setprecision(4) << porosity
                      << "  mass_drift " << std::scientific << std::setprecision(1) << drift
                      << std::fixed << std::setprecision(3)
                      << "  c in [" << cmin << ", " << cmax << "]"
                      << std::defaultfloat << std::endl;
    }

    void saveMorphology(long s)
    {
        std::ofstream f("morphology_step_" + std::to_string(s) + ".raw", std::ios::binary);
        std::vector<float> buf(c.begin(), c.end());
        f.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(float));
    }

    void run()
    {
        initialize();
        std::ofstream f("bubble_metrics.csv");
        f << "step,time_dimensionless,time_seconds,time_hours,number_of_projected_bubbles,"
             "total_projected_area,largest_projected_area,largest_area_fraction,coverage,porosity,"
             "dt,gas_total,mass_drift,c_min,c_max\n";

        std::cout << "theta = " << P.theta << "  bS = " << P.bS << "  bGB = " << P.bGB
                  << "  Ms = " << P.Ms << "  Mgb = " << P.Mgb << "  Mb = " << P.Mb
                  << "  LGB = " << P.LGB
                  << "  tEnd = " << P.tEnd << "  dtMax = " << P.dtMax << "  safety = " << P.safety << "\n"
                  << "Gas: conserved Cahn-Hilliard, no source, no clamp. Initial total gas = "
                  << std::setprecision(10) << mass0 << std::defaultfloat << "\n";

        writeMetrics(f, 0, 0.0, 0.0, true);
        saveMorphology(0);

        // output times are k * outDt; snapshots every snapEvery outputs
        const long nOut = std::lround(P.tEnd / P.outDt);
        const long snapEvery = std::max(1L, std::lround(P.snapDt / P.outDt));
        const long printEvery = std::max(1L, std::lround(1.0 / P.outDt));

        double t = 0.0;
        long s = 0;
        double dtSum = 0.0, dtMin = 1e30; long dtCount = 0;
        for (long k = 1; k <= nOut; k++)
        {
            const double tOut = k * P.outDt;
            double dtReport = 0.0;   // last step not shortened to hit the output time
            while (t < tOut)
            {
                double left = tOut - t;
                if (left <= 1e-9 * P.outDt) { t = tOut; break; }   // round-off only, no step needed
                double dt = step(left);
                s++;
                if (dt >= left) t = tOut;  // landed exactly on the output time
                else { t += dt; dtReport = dt; dtSum += dt; dtMin = std::min(dtMin, dt); dtCount++; }
            }
            writeMetrics(f, s, t, dtReport, k % printEvery == 0);
            if (k % snapEvery == 0) saveMorphology(s);
        }
        std::cout << "Done. " << s << " steps. dt (excluding output-aligned steps): mean "
                  << std::scientific << std::setprecision(3) << (dtCount ? dtSum / dtCount : 0.0)
                  << ", min " << dtMin << std::defaultfloat << "\n"
                  << "Output: bubble_metrics.csv, morphology_step_*.raw\n";
    }

private:
    Params P;
    std::vector<double> c, n1, n2, cN, n1N, n2N, mu, mob, lapN1, lapN2;
    int xm[NX], xp[NX], ym[NY], yp[NY], zm[NZ], zp[NZ];
    std::mt19937 rng;
    double mass0 = 0.0;
};

int main(int argc, char** argv)
{
    Params p;
    p.setTheta(60);
    for (int i = 1; i < argc; i++)
    {
        std::string a = argv[i];
        auto eq = a.find('=');
        if (eq == std::string::npos) continue;
        std::string k = a.substr(0, eq);
        double v = std::atof(a.c_str() + eq + 1);
        if (k == "theta") p.setTheta((int)v);
        else if (k == "tEnd") p.tEnd = v;
        else if (k == "outDt") p.outDt = v;
        else if (k == "snapDt") p.snapDt = v;
        else if (k == "dtMax") p.dtMax = v;
        else if (k == "safety") p.safety = v;
        else if (k == "minThick") p.minThick = (int)v;
        else if (k == "Mgb") p.Mgb = v;
        else if (k == "Ms") p.Ms = v;
        else if (k == "Mb") p.Mb = v;
        else if (k == "LGB") p.LGB = v;
        else if (k == "seed") p.seed = (unsigned)v;
        else if (k == "cInit") p.cInit = v;
        else if (k == "r0") p.r0 = v;
        else if (k == "nBubbles") p.nBubbles = (int)v;
        else if (k == "minSep") p.minSep = v;
        else { std::cerr << "unknown parameter " << k << "\n"; return 1; }
    }
    Simulation sim(p);
    sim.run();
    return 0;
}