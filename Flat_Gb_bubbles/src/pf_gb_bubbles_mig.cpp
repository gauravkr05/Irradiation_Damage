/*
    phasefield_gb_bubbles3d.cpp   (conserved version, GB-energy study)

    3-D phase-field simulation of gas bubble growth and coalescence on a
    flat grain boundary (GB) between two grains.

    Domain   : 128 x 128 x 64, dx = 10 nm   (compile-time; see GRID below)
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
    CHANGES FOR THE GRAIN-BOUNDARY-ENERGY STUDY (everything else unchanged)
    ------------------------------------------------------------------

    1. bGB is now an independent input. j, jGB, bS are NOT re-derived from it.
       In the old setTheta() a change of bGB also changed jGB = 2 bGB, and
       jGB is the gradient coefficient of n1, n2 in the bubble/grain
       interface as well (n goes 1 -> 0 there). So the old parameter sets
       changed the bubble/grain energy together with the GB energy.
       With j, jGB, bS, A0 fixed and only bGB varied:
         - bubble/grain energy is exactly unchanged (n2 = 0 on that
           interface, so the bGB n1^2 n2^2 term is identically zero there),
         - bulk states (grain: c=0,n=1; bubble: c=1,n=0; f = -1/4) are
           unchanged,
         - only the GB energy (and GB width) changes.
       bGB must be > 1: for bGB = 1 the (n1,n2) potential depends only on
       n1^2+n2^2, so the GB has no barrier and spreads without limit; for
       bGB < 1 a mixed n1 = n2 state is lower in energy than either grain.
       A warning is printed in those cases. setTheta() is kept only for
       reproducing old runs (theta=...) and prints a warning.

    2. Common initial condition: initOut=FILE writes the generated
       (c, n1, n2) fields and bubble centres to a binary file; initIn=FILE
       loads them, so every case starts from bit-identical fields.
       initOnly=1 writes the file and exits. centresIn=FILE (csv "x,y")
       places bubbles at given centres instead of random ones.

    3. outDir=DIR : all output goes to DIR. A params.txt and a
       snapshots.csv (time <-> step <-> file) are written there. Snapshot
       files keep their original names (morphology_step_N.raw); since
       step numbers differ between cases, snapshots.csv gives the times,
       which are identical across cases (multiples of snapDt).
       grains_step_N.raw (float32 n1 then n2, each [NZ][NY][NX]) is also
       written so GB integrity can be inspected.

    4. Extra columns are APPENDED to bubble_metrics.csv (old columns are
       unchanged and in the same order):
         free_energy   total F (should not increase)
         gb_z_mean     mean height of the n1 = n2 surface over bubble-free
                       columns (initially 31.5)
         gb_z_std      its standard deviation
         gb_delta      mean of sum_z 16 n1^2 n2^2 over bubble-free
                       columns (effective GB width seen by the Mgb term)
         n_min, n_max  range of n1, n2
         nonfinite     number of non-finite values in c, n1, n2

    5. Checkpoint / restart (bit-exact): checkpoint.bin is written at every
       snapshot time and at the end (fields in double precision, time,
       step, initial gas inventory and run statistics). restartIn=FILE
       continues the run to the new tEnd, appending to bubble_metrics.csv
       and snapshots.csv. Because dt depends only on the fields and the
       next output time, a restarted run is identical to an uninterrupted
       one. mass_drift stays relative to the original t = 0 inventory.

    6. Grid size can be overridden at compile time (-DGRID_NX=.. etc.),
       used only for the contact-angle and critical-radius test runs.
       Defaults are unchanged. ZGB = NZ/2 - 0.5.

    ------------------------------------------------------------------
    CHANGES FOR THE GRAIN-BOUNDARY MIGRATION EXPERIMENT (gb_migration/)
    ------------------------------------------------------------------
    7. Bulk free-energy bias between the grains (parameter dfBias, default 0):

           F  ->  F + integral dfBias * h(n1,n2) dV,   h = n1^2 / (n1^2 + n2^2)

       h = 1 in grain I (n1=1, n2=0) and 0 in grain II (n1=0, n2=1), so
           f_I - f_II = dfBias   exactly  (per unit volume, same units as A0 f).
       dh/dn1 = 2 n1 n2^2 / s^2, dh/dn2 = -2 n1^2 n2 / s^2 (s = n1^2+n2^2) both
       vanish at the pure grain states, which therefore remain exact
       stationary states with unchanged order-parameter values. The term
       enters only the Allen-Cahn equations (dF/dn_i += dfBias dh/dn_i), the
       free-energy diagnostic, and (as |dfBias d2h/dn_i^2|) the dt_AC bound.
       h is guarded for s < 1e-12 (only possible inside a bubble; none here).
       With dfBias = 0 every code path is the original one.

    8. gbOnly=1: bubble-free diagnostics. bubble_metrics.csv (coverage,
       components) is NOT computed; gb_position.csv is written instead with,
       per output, the n1 = n2 crossing of every (x,y) column (either grain
       below), its mean, std, min, max, the number of columns with none or
       several crossings, the GB width integral of 16 n1^2 n2^2, min of
       n1^2 + n2^2, dfBias, F, gas total, c/n ranges and non-finite count.
       Snapshots contain only the grain fields (grains_step_N.raw).

    ------------------------------------------------------------------
    CHANGES FOR THE COUPLED GB-MIGRATION + BUBBLE EXPERIMENT (Flat_Gb_bubbles/)
    ------------------------------------------------------------------
    9. zgb0=Z : height of the planar GB AND of the bubble centres used by the
       generator (default NZ/2 - 0.5 = 31.5, i.e. unchanged). Only generate()
       uses it. Grain I (n1) stays below, grain II (n2) above, as in the original.

   10. biasMode (default 0 = item 7, unchanged).
       biasMode=1: polynomial bias on grain II,
           F -> F + integral dfBias * g(n2) dV,   g(n) = n^2 (3 - 2 n)
       g(0)=0, g(1)=1, g'(0)=g'(1)=0, so grain I (n2=0), grain II (n2=1) and the
       bubble (n1=n2=0) all stay exact stationary states, and
           f_II - f_I = dfBias,   f_bubble - f_I = 0 .
       Grain II (upper) is the high-energy grain, so the GB migrates UP.
       Why not item 7: h = n1^2/(n1^2+n2^2) is undefined at n1=n2=0, which is
       exactly the bubble interior; dh/dn ~ 1/|n| and d2h/dn2 ~ 1/|n|^2 there,
       so the explicit AC update becomes singular and dt -> 0. Any bias that is
       continuous must give the bubble a single value; g(n2) gives it the value
       of the growing grain I, so a bubble that ends up in grain I behind the GB
       feels exactly the original (unbiased) thermodynamics. A bubble surface in
       contact with grain II gains dfBias per unit volume converted, the same
       amount the GB gains: converting grain II into grain I through a pore is
       physically the same driving force.
       Derivative: dF/dn2 += dfBias * 6 n2 (1 - n2); dt_AC bound uses
       |dfBias * (6 - 12 n2)|.

   11. Coupled diagnostics (bubble mode, every output): gb_bubble.csv with
         gb_real_*   : GB height over all columns that contain a REAL GB
                       crossing (n1-n2 changes sign from + to - between two
                       cells with n1+n2 > 0.5 in both), mean/std/min/max/count
         cols_*      : covered columns (>= minThick cells with c >= 0.5) in
                       which the GB passes through the bubble (no real crossing),
                       lies above the topmost bubble cell (bubble left behind in
                       grain I), or lies below the lowest bubble cell
         gb_coverage : covered columns through which the GB passes / NX*NY
         bub_z_*     : mean z of all c >= 0.5 cells, mean bottom/top z of the
                       bubble over covered columns
         vol_II      : number of cells with n2 > n1 and n1+n2 > 0.5
       and, at every snapshot time, maps/tTTTT.TT.bin (compact):
         uint8 thick[NY][NX], float32 gbz[NY][NX] (real crossing, NaN if none),
         uint8 zbot[NY][NX], uint8 ztop[NY][NX] (255 if no bubble cell),
         float32 c[NZ][NX] and float32 (n1-n2)[NZ][NX] on the plane y = sliceY.
       Full raw fields (morphology/grains .raw) are written only every
       rawEvery-th snapshot (default 1 = every snapshot, as before).

    ------------------------------------------------------------------
    SPHERICAL GRAIN (sphereR > 0)   (Flat_Gb_bubbles/sphere study)
    ------------------------------------------------------------------
   12. Geometry: an inner spherical grain of radius sphereR, centred at the box
       centre ((NX-1)/2, (NY-1)/2, (NZ-1)/2), surrounded by the outer grain.
       The INNER grain is grain II (n2), the OUTER grain is grain I (n1), so with
       biasMode=1 and dfBias > 0 the inner grain is the high-energy one and
       shrinks faster; with dfBias = 0 it shrinks by curvature only.
       Same tanh profile (w = 1.5). Bubbles (nBubbles, r0) are placed at random
       points ON the sphere surface with chord distance >= minSep (seed), as
       lenses straddling the GB, exactly like the flat generator.
   13. Diagnostics in sphere mode (sphere.csv, every output), measured along
       nDir Fibonacci-distributed rays from the centre (trilinear interpolation,
       step 0.25):
         r_gb_mean/std/min/max, n_dir_gb : first REAL GB crossing (n2-n1 from
                       + to - with n1+n2 > 0.5 on both sides) per ray
         coverage    : fraction of rays with >= minThick cells of c >= 0.5
                       (solid-angle fraction of the GB covered by gas)
         n_patches, largest_fraction : connected components of covered rays
                       (6-nearest-neighbour graph on the direction set)
         dirs_gb_through_gas : covered rays without a real crossing
         V_inner, R_vol : cells with n2 > n1 and n1+n2 > 0.5, and the radius of
                       the sphere with that volume
         centroid_offset : distance of the inner-grain centroid from the box
                       centre; gas_cells, gas_r_mean; gas_total, mass_drift,
                       free_energy, c/n range, nonfinite
       maps/tTTTT.TT.bin at snapshot times: float32 c[NZ][NX] and
       (n2-n1)[NZ][NX] on the plane y = (NY-1)/2 (rounded down), float32
       r_gb[nDir] (NaN if none), uint8 covered[nDir]; dirs.bin = float32
       xyz[nDir]. stopAfterVanish = T (>= 0): stop T time units after the inner
       grain has vanished (V_inner = 0). bubble_metrics.csv is not written.

    ------------------------------------------------------------------
    GAS CONSERVATION   (unchanged)
    ------------------------------------------------------------------
      * flux-form update, fluxes cancel pairwise, zero flux through z walls
      * no source term
      * c is NOT clamped after a step (mobility uses c clamped to [0,1],
        which only scales fluxes)
    Every metrics row records the total gas and its relative drift.

    ------------------------------------------------------------------
    TIME STEP (adaptive, from the stability limit)   (unchanged)
    ------------------------------------------------------------------
        dt_CH <= 2 / [ M_max ( 12 A0 max(f_cc) + 144 j ) ]
        dt_AC <= 2 / [ LGB ( A0 max(f_nn) + 12 jGB ) ]
        dt = min( dtMax, safety * dt_CH, safety * dt_AC, time to next output )

    ------------------------------------------------------------------
    PARAMETERS FROM THE COMMAND LINE (all optional)
    ------------------------------------------------------------------
        bGB, bS, j, jGB, A0       free-energy coefficients
                                  (defaults: 1, 1, 2, 2, 4 = old theta=60)
        theta    LEGACY: sets bS, bGB, j=2bS, jGB=2bGB as before
        tEnd, outDt, snapDt, dtMax, safety, minThick,
        Mgb, Ms, Mb, LGB, cInit, r0, nBubbles, minSep, seed, noise
        outDir, initIn, initOut, initOnly, centresIn, saveGrains (0/1),
        restartIn, ckpt (0/1)

    Compile:
        g++ -O3 -march=native -fopenmp -std=c++17 phasefield_gb_bubbles3d.cpp -o phasefield_gb_bubbles3d
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// ------------------------------------------------------------
// Grid
// ------------------------------------------------------------
#ifndef GRID_NX
#define GRID_NX 128
#endif
#ifndef GRID_NY
#define GRID_NY 128
#endif
#ifndef GRID_NZ
#define GRID_NZ 64
#endif
constexpr int NX = GRID_NX;
constexpr int NY = GRID_NY;
constexpr int NZ = GRID_NZ;
constexpr int N  = NX * NY * NZ;
constexpr double ZGB = NZ / 2 - 0.5;   // 31.5 for NZ = 64

constexpr double DX_M = 10.0e-9;      // one grid cell = 10 nm
constexpr double D_PHYS = 3.4e-18;    // m^2/s, from the paper (time conversion only)

// ------------------------------------------------------------
// Parameters (defaults; can be overridden on the command line)
// ------------------------------------------------------------
struct Params
{
    int theta = 0;          // 0 = not set (legacy label only)

    double A0 = 4.0;
    double bS = 1.0;
    double bGB = 1.0;
    double j = 2.0;
    double jGB = 2.0;

    double Ms = 0.05;
    double Mgb = 0.5;
    double Mb = 5.0;
    double LGB = 1.0;

    double tEnd = 20.0;
    double outDt = 0.1;
    double snapDt = 3.0;
    double dtMax = 1.0e-3;
    double safety = 0.5;
    int minThick = 3;

    int nBubbles = 33;
    double r0 = 3.8;
    double minSep = 19.0;
    double cInit = 0.30;
    double noise = 0.001;
    unsigned seed = 11;

    std::string outDir = ".";
    std::string initIn, initOut, centresIn, restartIn;
    bool ckpt = true;       // write checkpoint.bin at every snapshot time and at the end
    bool initOnly = false;
    bool saveGrains = true;
    double dfBias = 0.0;     // bulk bias: f_I - f_II = dfBias (GB-migration experiment)
    bool gbOnly = false;     // bubble-free GB diagnostics only
    double zgb0 = NZ / 2 - 0.5;  // generator: GB plane and bubble-centre height
    int biasMode = 0;        // 0: Moelans h on grain I (item 7); 1: dfBias*g(n2) (item 10)
    int rawEvery = 1;        // write full raw fields every rawEvery-th snapshot
    int sliceY = NY / 2;     // y index of the xz slice stored in maps/
    double sphereR = 0.0;    // > 0: spherical inner grain (item 12)
    int nDir = 6000;         // rays for the sphere diagnostics (item 13)
    double stopAfterVanish = -1.0;  // sphere mode: stop this long after the inner grain is gone

    // LEGACY: reproduces the old coupled parameter sets. Not used in the study.
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
        std::cerr << "WARNING: theta=" << t << " is a legacy label. It changes bS, bGB, j and jGB together;\n"
                     "         the label is NOT a validated contact angle.\n";
    }

    std::string describe() const
    {
        std::ostringstream o;
        o << std::setprecision(10)
          << "grid NX NY NZ        = " << NX << " " << NY << " " << NZ << "\n"
          << "A0 bS bGB j jGB      = " << A0 << " " << bS << " " << bGB << " " << j << " " << jGB << "\n"
          << "Ms Mgb Mb LGB        = " << Ms << " " << Mgb << " " << Mb << " " << LGB << "\n"
          << "tEnd outDt snapDt    = " << tEnd << " " << outDt << " " << snapDt << "\n"
          << "dtMax safety         = " << dtMax << " " << safety << "\n"
          << "minThick             = " << minThick << "\n"
          << "nBubbles r0 minSep   = " << nBubbles << " " << r0 << " " << minSep << "\n"
          << "cInit noise seed     = " << cInit << " " << noise << " " << seed << "\n"
          << "initIn               = " << (initIn.empty() ? "(generated)" : initIn) << "\n"
          << "centresIn            = " << (centresIn.empty() ? "(random)" : centresIn) << "\n"
          << "restart from         = " << (restartIn.empty() ? "(none)" : restartIn) << "\n"
          << "legacy theta label   = " << (theta ? std::to_string(theta) : "none") << "\n"
          << "dfBias gbOnly        = " << dfBias << " " << gbOnly << "\n"
          << "zgb0 biasMode        = " << zgb0 << " " << biasMode << "\n"
          << "rawEvery sliceY      = " << rawEvery << " " << sliceY << "\n"
          << "sphereR nDir stopAV  = " << sphereR << " " << nDir << " " << stopAfterVanish << "\n";
        return o.str();
    }
};

static inline int idx(int x, int y, int z) { return (z * NY + y) * NX + x; }

// bulk bias interpolation h = n1^2/(n1^2+n2^2): 1 in grain I, 0 in grain II
static inline double biasH(double a, double b)
{
    double s = a * a + b * b;
    return s > 1e-12 ? a * a / s : 0.5;
}
static inline void biasDerivs(double a, double b, double& h1, double& h2, double& h11, double& h22)
{
    double s = a * a + b * b;
    if (s <= 1e-12) { h1 = h2 = h11 = h22 = 0.0; return; }
    double s2 = s * s, s3 = s2 * s;
    h1 = 2.0 * a * b * b / s2;
    h2 = -2.0 * a * a * b / s2;
    h11 = 2.0 * b * b * (b * b - 3.0 * a * a) / s3;
    h22 = -2.0 * a * a * (a * a - 3.0 * b * b) / s3;
}

// item 10: polynomial bias on grain II, g(n) = n^2 (3 - 2n)
static inline double polyG(double n) { return n * n * (3.0 - 2.0 * n); }
static inline double polyG1(double n) { return 6.0 * n * (1.0 - n); }
static inline double polyG2(double n) { return 6.0 - 12.0 * n; }

class Simulation
{
public:
    explicit Simulation(const Params& p)
        : P(p), c(N), n1(N), n2(N), cN(N), n1N(N), n2N(N),
          mu(N), mob(N), lapN1(N), lapN2(N), rng(p.seed)
    {
        for (int i = 0; i < NX; i++) { xm[i] = (i + NX - 1) % NX; xp[i] = (i + 1) % NX; }
        for (int i = 0; i < NY; i++) { ym[i] = (i + NY - 1) % NY; yp[i] = (i + 1) % NY; }
        for (int i = 0; i < NZ; i++) { zm[i] = std::max(i - 1, 0); zp[i] = std::min(i + 1, NZ - 1); }
    }

    std::string path(const std::string& name) const { return P.outDir + "/" + name; }

    // ---------------- initial condition: generate (unchanged algorithm) ----------------
    void generateSphere()
    {
        const double cx = (NX - 1) / 2.0, cy = (NY - 1) / 2.0, cz = (NZ - 1) / 2.0, R0 = P.sphereR;
        std::uniform_real_distribution<double> uz(-1.0, 1.0), uphi(0.0, 2.0 * M_PI);
        centres.clear(); centres3d.clear();
        int attempts = 0;
        while ((int)centres3d.size() < P.nBubbles && attempts++ < 1000000)
        {
            double z = uz(rng), phi = uphi(rng), rr = std::sqrt(std::max(0.0, 1.0 - z * z));
            double px = cx + R0 * rr * std::cos(phi), py = cy + R0 * rr * std::sin(phi), pz = cz + R0 * z;
            bool ok = true;
            for (auto& q : centres3d)
            {
                double dx = px - q[0], dy = py - q[1], dz = pz - q[2];
                if (dx * dx + dy * dy + dz * dz < P.minSep * P.minSep) { ok = false; break; }
            }
            if (ok) { centres3d.push_back({px, py, pz}); centres.push_back({px, py}); }
        }
        std::cout << "Initial bubbles placed on the sphere = " << centres3d.size() << " (requested " << P.nBubbles << ")\n";
        std::uniform_real_distribution<double> un(-P.noise, P.noise);
        const double w = 1.5;
        for (int z = 0; z < NZ; z++)
            for (int y = 0; y < NY; y++)
                for (int x = 0; x < NX; x++)
                {
                    int p = idx(x, y, z);
                    double r = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz));
                    double b = 0.0;
                    for (auto& q : centres3d)
                    {
                        double d = std::sqrt((x - q[0]) * (x - q[0]) + (y - q[1]) * (y - q[1]) + (z - q[2]) * (z - q[2]));
                        b = std::max(b, 0.5 * (1.0 - std::tanh((d - P.r0) / w)));
                    }
                    double g2 = 0.5 * (1.0 - std::tanh((r - R0) / w));   // inner grain = n2
                    double g1 = 0.5 * (1.0 + std::tanh((r - R0) / w));   // outer grain = n1
                    n1[p] = g1 * (1.0 - b);
                    n2[p] = g2 * (1.0 - b);
                    double v = P.cInit + (1.0 - P.cInit) * b + (P.noise > 0 ? un(rng) : 0.0);
                    c[p] = std::clamp(v, 0.0, 1.0);
                }
    }

    void generate()
    {
        if (P.sphereR > 0.0) { generateSphere(); return; }
        std::uniform_real_distribution<double> ux(0.0, NX), uy(0.0, NY);
        centres.clear();
        if (!P.centresIn.empty())
        {
            std::ifstream f(P.centresIn);
            std::string line;
            while (std::getline(f, line))
            {
                if (line.empty() || line[0] == '#' || line[0] == 'x') continue;
                std::replace(line.begin(), line.end(), ',', ' ');
                std::istringstream s(line);
                double x, y;
                if (s >> x >> y) centres.push_back({x, y});
            }
        }
        else
        {
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
        }
        std::cout << "Initial bubbles placed = " << centres.size() << "\n";

        std::uniform_real_distribution<double> un(-P.noise, P.noise);
        const double w = 1.5;

        for (int z = 0; z < NZ; z++)
            for (int y = 0; y < NY; y++)
                for (int x = 0; x < NX; x++)
                {
                    int p = idx(x, y, z);
                    double dz = z - P.zgb0;
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
                    double g1 = 0.5 * (1.0 - std::tanh(dz / w));
                    double g2 = 0.5 * (1.0 + std::tanh(dz / w));
                    n1[p] = g1 * (1.0 - b);
                    n2[p] = g2 * (1.0 - b);
                    double v = P.cInit + (1.0 - P.cInit) * b + (P.noise > 0 ? un(rng) : 0.0);
                    c[p] = std::clamp(v, 0.0, 1.0);
                }
    }

    // binary IC file: header, centres, then c, n1, n2 as double
    void saveInitial(const std::string& fn) const
    {
        std::ofstream f(fn, std::ios::binary);
        const char magic[8] = {'P','F','G','B','I','C','0','1'};
        int32_t dims[3] = {NX, NY, NZ};
        int32_t nc = (int32_t)centres.size();
        double meta[4] = {P.r0, P.cInit, P.noise, (double)P.seed};
        f.write(magic, 8);
        f.write((const char*)dims, sizeof dims);
        f.write((const char*)meta, sizeof meta);
        f.write((const char*)&nc, sizeof nc);
        for (auto& q : centres) { double xy[2] = {q.first, q.second}; f.write((const char*)xy, sizeof xy); }
        f.write((const char*)c.data(), N * sizeof(double));
        f.write((const char*)n1.data(), N * sizeof(double));
        f.write((const char*)n2.data(), N * sizeof(double));
        std::ofstream t(fn + ".centres.csv");
        t << "x,y\n" << std::setprecision(17);
        for (auto& q : centres) t << q.first << "," << q.second << "\n";
        if (!centres3d.empty())
        {
            std::ofstream t3(fn + ".centres3d.csv");
            t3 << "x,y,z\n" << std::setprecision(17);
            for (auto& q : centres3d) t3 << q[0] << "," << q[1] << "," << q[2] << "\n";
        }
    }

    void loadInitial(const std::string& fn)
    {
        std::ifstream f(fn, std::ios::binary);
        if (!f) { std::cerr << "cannot open initial condition " << fn << "\n"; std::exit(1); }
        char magic[8]; int32_t dims[3]; double meta[4]; int32_t nc;
        f.read(magic, 8);
        f.read((char*)dims, sizeof dims);
        if (std::memcmp(magic, "PFGBIC01", 8) || dims[0] != NX || dims[1] != NY || dims[2] != NZ)
        { std::cerr << "initial condition " << fn << " does not match this grid\n"; std::exit(1); }
        f.read((char*)meta, sizeof meta);
        f.read((char*)&nc, sizeof nc);
        centres.resize(nc);
        for (auto& q : centres) { double xy[2]; f.read((char*)xy, sizeof xy); q = {xy[0], xy[1]}; }
        f.read((char*)c.data(), N * sizeof(double));
        f.read((char*)n1.data(), N * sizeof(double));
        f.read((char*)n2.data(), N * sizeof(double));
        if (!f) { std::cerr << "initial condition " << fn << " is truncated\n"; std::exit(1); }
        P.r0 = meta[0]; P.cInit = meta[1]; P.noise = meta[2]; P.seed = (unsigned)meta[3];
        P.nBubbles = nc;
        std::cout << "Loaded initial condition " << fn << " (" << nc << " bubbles, r0 = " << P.r0
                  << ", cInit = " << P.cInit << ")\n";
    }

    void initialize()
    {
        if (!P.initIn.empty()) loadInitial(P.initIn);
        else generate();
        if (!P.initOut.empty()) saveInitial(P.initOut);
        mass0 = totalGas();
    }

    void saveCheckpoint(double t, long s) const
    {
        std::string fn = path("checkpoint.bin"), tmp = fn + ".tmp";
        {
            std::ofstream f(tmp, std::ios::binary);
            const char magic[8] = {'P','F','G','B','C','K','0','1'};
            int32_t dims[3] = {NX, NY, NZ};
            int32_t nc = (int32_t)centres.size();
            double st[8] = {t, (double)s, mass0, lastF, maxAbsDrift, runCmin, runCmax, (double)Fincreases};
            double meta[4] = {P.r0, P.cInit, P.noise, (double)P.seed};
            int32_t mc = maxCount;
            f.write(magic, 8); f.write((const char*)dims, sizeof dims);
            f.write((const char*)st, sizeof st); f.write((const char*)meta, sizeof meta);
            f.write((const char*)&mc, sizeof mc); f.write((const char*)&nc, sizeof nc);
            for (auto& q : centres) { double xy[2] = {q.first, q.second}; f.write((const char*)xy, sizeof xy); }
            f.write((const char*)c.data(), N * sizeof(double));
            f.write((const char*)n1.data(), N * sizeof(double));
            f.write((const char*)n2.data(), N * sizeof(double));
        }
        std::filesystem::rename(tmp, fn);
    }

    void loadCheckpoint(const std::string& fn, double& t, long& s)
    {
        std::ifstream f(fn, std::ios::binary);
        char magic[8]; int32_t dims[3], mc, nc; double st[8], meta[4];
        f.read(magic, 8); f.read((char*)dims, sizeof dims);
        if (!f || std::memcmp(magic, "PFGBCK01", 8) || dims[0] != NX || dims[1] != NY || dims[2] != NZ)
        { std::cerr << "checkpoint " << fn << " is missing or does not match this grid\n"; std::exit(1); }
        f.read((char*)st, sizeof st); f.read((char*)meta, sizeof meta);
        f.read((char*)&mc, sizeof mc); f.read((char*)&nc, sizeof nc);
        centres.resize(nc);
        for (auto& q : centres) { double xy[2]; f.read((char*)xy, sizeof xy); q = {xy[0], xy[1]}; }
        f.read((char*)c.data(), N * sizeof(double));
        f.read((char*)n1.data(), N * sizeof(double));
        f.read((char*)n2.data(), N * sizeof(double));
        if (!f) { std::cerr << "checkpoint " << fn << " is truncated\n"; std::exit(1); }
        t = st[0]; s = std::lround(st[1]); mass0 = st[2]; lastF = st[3]; maxAbsDrift = st[4];
        runCmin = st[5]; runCmax = st[6]; Fincreases = std::lround(st[7]); maxCount = mc;
        P.r0 = meta[0]; P.cInit = meta[1]; P.noise = meta[2]; P.seed = (unsigned)meta[3]; P.nBubbles = nc;
        std::cout << "Restarting from " << fn << " at t = " << t << ", step " << s << "\n";
    }

    double totalGas() const
    {
        double s = 0.0;
        #pragma omp parallel for reduction(+:s) schedule(static)
        for (int p = 0; p < N; p++) s += c[p];
        return s;
    }

    // total free energy, same discretisation as the evolution equations
    // (face differences; zero-flux z walls contribute no face)
    double freeEnergy() const
    {
        const double A0 = P.A0, bS = P.bS, bGB = P.bGB, j = P.j, jGB = P.jGB;
        double F = 0.0;
        #pragma omp parallel for collapse(2) reduction(+:F) schedule(static)
        for (int z = 0; z < NZ; z++)
            for (int y = 0; y < NY; y++)
                for (int x = 0; x < NX; x++)
                {
                    int p = idx(x, y, z);
                    double cv = c[p], a = n1[p], b = n2[p];
                    double fb = 0.25 * cv * cv * cv * cv - 0.5 * cv * cv
                              + 0.25 * a * a * a * a - 0.5 * a * a + 0.25 * b * b * b * b - 0.5 * b * b
                              + 0.5 * bGB * a * a * b * b + bS * cv * cv * (a * a + b * b);
                    double g = 0.0;
                    int q[3] = {idx(xp[x], y, z), idx(x, yp[y], z), z < NZ - 1 ? idx(x, y, z + 1) : -1};
                    for (int k = 0; k < 3; k++)
                    {
                        if (q[k] < 0) continue;
                        if (k == 1 && NY == 1) continue;
                        if (k == 0 && NX == 1) continue;
                        double dc = c[q[k]] - cv, d1 = n1[q[k]] - a, d2 = n2[q[k]] - b;
                        g += 0.5 * j * dc * dc + 0.5 * jGB * (d1 * d1 + d2 * d2);
                    }
                    double hb = (P.dfBias == 0.0) ? 0.0
                              : (P.biasMode == 1 ? P.dfBias * polyG(b) : P.dfBias * biasH(a, b));
                    F += A0 * fb + g + hb;
                }
        return F;
    }

    double step(double dtLimit)
    {
        const double A0 = P.A0, bS = P.bS, bGB = P.bGB, j = P.j, jGB = P.jGB;
        const double Ms = P.Ms, Mgb = P.Mgb, Mb = P.Mb;
        const double dfB = P.dfBias;

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

                    double cc = std::clamp(cv, 0.0, 1.0);
                    double m = Ms * 16.0 * cc * cc * (1.0 - cc * cc)
                             + Mgb * 16.0 * a * a * b * b
                             + Mb;
                    mob[p] = m;

                    double fcc = 3.0 * cv * cv - 1.0 + 2.0 * bS * (a * a + b * b);
                    double f11 = 3.0 * a * a - 1.0 + bGB * b * b + 2.0 * bS * cv * cv;
                    double f22 = 3.0 * b * b - 1.0 + bGB * a * a + 2.0 * bS * cv * cv;
                    if (dfB != 0.0 && P.biasMode == 1)
                        f22 += std::abs(dfB * polyG2(b)) / A0;
                    else if (dfB != 0.0)
                    {
                        double h1, h2, h11, h22; biasDerivs(a, b, h1, h2, h11, h22);
                        f11 += std::abs(dfB * h11) / A0; f22 += std::abs(dfB * h22) / A0;
                    }

                    Mmax = std::max(Mmax, m);
                    fccMax = std::max(fccMax, fcc);
                    fnnMax = std::max(fnnMax, std::max(f11, f22));
                }

        const double L = P.LGB;
        double dtCH = 2.0 / (Mmax * (12.0 * A0 * fccMax + 144.0 * j));
        double dtAC = 2.0 / (L * (A0 * fnnMax + 12.0 * jGB));
        if (!std::isfinite(dtCH) || !std::isfinite(dtAC) || !(dtCH > 0.0) || !(dtAC > 0.0))
        {
            std::cerr << "Fields became non-finite (numerical blow-up). Reduce safety or dtMax.\n";
            std::exit(2);
        }
        lastDtCH = dtCH; lastDtAC = dtAC;
        double dt = std::min({P.dtMax, P.safety * dtCH, P.safety * dtAC, dtLimit});

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
                    if (z < NZ - 1) { int q = idx(x, y, z + 1); div += 0.5 * (m + mob[q]) * (mu[q] - u); }
                    if (z > 0)      { int q = idx(x, y, z - 1); div -= 0.5 * (m + mob[q]) * (u - mu[q]); }

                    double cv = c[p], a = n1[p], b = n2[p];
                    cN[p] = cv + dt * div;        // no clamp: gas is conserved exactly

                    double dF1 = A0 * (a * a * a - a + bGB * a * b * b + 2.0 * bS * cv * cv * a) - jGB * lapN1[p];
                    double dF2 = A0 * (b * b * b - b + bGB * b * a * a + 2.0 * bS * cv * cv * b) - jGB * lapN2[p];
                    if (dfB != 0.0 && P.biasMode == 1)
                        dF2 += dfB * polyG1(b);
                    else if (dfB != 0.0)
                    {
                        double h1, h2, h11, h22; biasDerivs(a, b, h1, h2, h11, h22);
                        dF1 += dfB * h1; dF2 += dfB * h2;
                    }
                    n1N[p] = std::clamp(a - dt * L * dF1, 0.0, 1.0);
                    n2N[p] = std::clamp(b - dt * L * dF2, 0.0, 1.0);
                }

        c.swap(cN);
        n1.swap(n1N);
        n2.swap(n2N);
        return dt;
    }

    void writeMetrics(std::ofstream& f, long s, double t, double dt, bool print)
    {
        std::vector<unsigned char> mask(NX * NY, 0);
        long poreCells = 0;
        // GB diagnostics over bubble-free columns
        double gbZ = 0.0, gbZ2 = 0.0, gbDelta = 0.0; long gbCols = 0;
        #pragma omp parallel for reduction(+:poreCells, gbZ, gbZ2, gbDelta, gbCols) collapse(2)
        for (int y = 0; y < NY; y++)
            for (int x = 0; x < NX; x++)
            {
                int thick = 0;
                for (int z = 0; z < NZ; z++)
                    if (c[idx(x, y, z)] >= 0.5) thick++;
                poreCells += thick;
                if (thick >= P.minThick) mask[y * NX + x] = 1;
                if (thick == 0)
                {
                    double zc = -1.0, del = 0.0;
                    for (int z = 0; z < NZ; z++)
                    {
                        int p = idx(x, y, z);
                        del += 16.0 * n1[p] * n1[p] * n2[p] * n2[p];
                        if (z + 1 < NZ && zc < 0.0)
                        {
                            double d0 = n1[p] - n2[p], d1 = n1[idx(x, y, z + 1)] - n2[idx(x, y, z + 1)];
                            if (d0 > 0.0 && d1 <= 0.0) zc = z + d0 / (d0 - d1);
                        }
                    }
                    if (zc >= 0.0) { gbZ += zc; gbZ2 += zc * zc; gbDelta += del; gbCols++; }
                }
            }

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

        double mass = 0.0, cmin = 1e30, cmax = -1e30, nmin = 1e30, nmax = -1e30;
        long bad = 0;
        #pragma omp parallel for reduction(+:mass, bad) reduction(min:cmin, nmin) reduction(max:cmax, nmax) schedule(static)
        for (int p = 0; p < N; p++)
        {
            mass += c[p]; cmin = std::min(cmin, c[p]); cmax = std::max(cmax, c[p]);
            nmin = std::min({nmin, n1[p], n2[p]}); nmax = std::max({nmax, n1[p], n2[p]});
            bad += !std::isfinite(c[p]) + !std::isfinite(n1[p]) + !std::isfinite(n2[p]);
        }
        double drift = (mass - mass0) / mass0;

        if (!std::isfinite(mass) || bad)
        {
            std::cerr << "Field became non-finite at step " << s << ", t = " << t
                      << ". Reduce safety or dtMax.\n";
            std::exit(2);
        }

        double F = freeEnergy();
        if (s > 0 && F > lastF + 1e-9 * std::abs(lastF)) Fincreases++;
        lastF = F;
        maxAbsDrift = std::max(maxAbsDrift, std::abs(drift));
        runCmin = std::min(runCmin, cmin); runCmax = std::max(runCmax, cmax);
        maxCount = std::max(maxCount, count);

        double tSec = t * DX_M * DX_M / D_PHYS;
        double coverage = double(total) / (NX * NY);
        double frac = total > 0 ? double(largest) / total : 0.0;
        double porosity = double(poreCells) / N;
        double gzm = gbCols ? gbZ / gbCols : NAN;
        double gzs = gbCols ? std::sqrt(std::max(0.0, gbZ2 / gbCols - gzm * gzm)) : NAN;
        double gdel = gbCols ? gbDelta / gbCols : NAN;

        f << s << "," << std::setprecision(10) << t << "," << tSec << "," << tSec / 3600.0 << ","
          << count << "," << total << "," << largest << "," << frac << ","
          << coverage << "," << porosity << ","
          << dt << "," << mass << "," << drift << "," << cmin << "," << cmax << ","
          << F << "," << gzm << "," << gzs << "," << gdel << "," << nmin << "," << nmax << "," << bad
          << std::endl;

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
                      << "  GB z " << std::setprecision(2) << gzm
                      << "  F " << std::setprecision(1) << F
                      << std::defaultfloat << std::endl;
    }

    // item 11: coupled GB/bubble diagnostics (bubble mode)
    void writeCoupled(std::ofstream& f, double t, bool writeMap)
    {
        std::vector<unsigned char> thick(NX * NY, 0), zb(NX * NY, 255), zt(NX * NY, 255);
        std::vector<float> gbz(NX * NY, NAN);
        double s = 0.0, s2 = 0.0, zmin = 1e30, zmax = -1e30, bz = 0.0, bbot = 0.0, btop = 0.0;
        long cnt = 0, multi = 0, cov = 0, through = 0, above = 0, below = 0, bcells = 0, volII = 0;
        #pragma omp parallel for collapse(2) reduction(+:s, s2, cnt, multi, cov, through, above, below, bz, bbot, btop, bcells, volII) reduction(min:zmin) reduction(max:zmax)
        for (int y = 0; y < NY; y++)
            for (int x = 0; x < NX; x++)
            {
                int th = 0, lo = -1, hi = -1, nc = 0; double zc = -1.0;
                for (int z = 0; z < NZ; z++)
                {
                    int p = idx(x, y, z);
                    double a = n1[p], b = n2[p];
                    if (c[p] >= 0.5) { th++; if (lo < 0) lo = z; hi = z; bz += z; bcells++; }
                    if (b > a && a + b > 0.5) volII++;
                    if (z + 1 < NZ)
                    {
                        int q = idx(x, y, z + 1);
                        double d0 = a - b, d1 = n1[q] - n2[q];
                        if (d0 > 0.0 && d1 <= 0.0 && a + b > 0.5 && n1[q] + n2[q] > 0.5)
                        { if (nc == 0) zc = z + d0 / (d0 - d1); nc++; }
                    }
                }
                int m = y * NX + x;
                thick[m] = (unsigned char)std::min(th, 255);
                if (lo >= 0) { zb[m] = (unsigned char)lo; zt[m] = (unsigned char)hi; }
                if (nc > 0) { gbz[m] = (float)zc; s += zc; s2 += zc * zc; cnt++; zmin = std::min(zmin, zc); zmax = std::max(zmax, zc); }
                if (nc > 1) multi++;
                if (th >= P.minThick)
                {
                    cov++; bbot += lo; btop += hi;
                    if (nc == 0) through++;
                    else if (zc > hi) above++;
                    else if (zc < lo) below++;
                }
            }
        double gm = cnt ? s / cnt : NAN;
        double gs = cnt ? std::sqrt(std::max(0.0, s2 / cnt - gm * gm)) : NAN;
        f << std::setprecision(10) << t << "," << gm << "," << gs << "," << (cnt ? zmin : NAN) << "," << (cnt ? zmax : NAN) << ","
          << cnt << "," << multi << "," << cov << "," << through << "," << above << "," << below << ","
          << double(through) / (NX * NY) << "," << (bcells ? bz / bcells : NAN) << ","
          << (cov ? bbot / cov : NAN) << "," << (cov ? btop / cov : NAN) << "," << volII << std::endl;

        if (writeMap)
        {
            std::filesystem::create_directories(path("maps"));
            std::ostringstream nm; nm << "maps/t" << std::fixed << std::setprecision(2) << std::setw(7) << std::setfill('0') << t << ".bin";
            std::ofstream g(path(nm.str()), std::ios::binary);
            g.write((const char*)thick.data(), NX * NY);
            g.write((const char*)gbz.data(), NX * NY * sizeof(float));
            g.write((const char*)zb.data(), NX * NY);
            g.write((const char*)zt.data(), NX * NY);
            std::vector<float> sc(NZ * NX), sn(NZ * NX);
            int y = std::clamp(P.sliceY, 0, NY - 1);
            for (int z = 0; z < NZ; z++)
                for (int x = 0; x < NX; x++)
                { int p = idx(x, y, z); sc[z * NX + x] = (float)c[p]; sn[z * NX + x] = (float)(n1[p] - n2[p]); }
            g.write((const char*)sc.data(), sc.size() * sizeof(float));
            g.write((const char*)sn.data(), sn.size() * sizeof(float));
        }
    }

    // ---------------- item 13: sphere diagnostics ----------------
    std::vector<std::array<double, 3>> dirs;
    std::vector<std::array<int, 6>> dirNb;
    bool innerVanished = false; double tVanish = -1.0;

    void setupDirs()
    {
        const int n = P.nDir; dirs.resize(n); dirNb.resize(n);
        const double ga = M_PI * (3.0 - std::sqrt(5.0));
        for (int i = 0; i < n; i++)
        {
            double z = 1.0 - 2.0 * (i + 0.5) / n, r = std::sqrt(std::max(0.0, 1.0 - z * z)), ph = ga * i;
            dirs[i] = {r * std::cos(ph), r * std::sin(ph), z};
        }
        #pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < n; i++)
        {
            double bd[6]; int bi[6];
            for (int k = 0; k < 6; k++) { bd[k] = 1e30; bi[k] = -1; }
            for (int j = 0; j < n; j++)
            {
                if (j == i) continue;
                double d = 0.0; for (int k = 0; k < 3; k++) d += (dirs[i][k] - dirs[j][k]) * (dirs[i][k] - dirs[j][k]);
                int k = 5; if (d >= bd[5]) continue;
                while (k > 0 && bd[k - 1] > d) { bd[k] = bd[k - 1]; bi[k] = bi[k - 1]; k--; }
                bd[k] = d; bi[k] = j;
            }
            for (int k = 0; k < 6; k++) dirNb[i][k] = bi[k];
        }
        std::ofstream g(path("dirs.bin"), std::ios::binary);
        std::vector<float> buf; for (auto& d : dirs) { buf.push_back((float)d[0]); buf.push_back((float)d[1]); buf.push_back((float)d[2]); }
        g.write((const char*)buf.data(), buf.size() * sizeof(float));
    }

    // trilinear interpolation of a field (x, y periodic; z clamped)
    double interp(const std::vector<double>& f, double x, double y, double z) const
    {
        int x0 = (int)std::floor(x), y0 = (int)std::floor(y), z0 = (int)std::floor(z);
        double fx = x - x0, fy = y - y0, fz = z - z0;
        int xa = ((x0 % NX) + NX) % NX, xb = (xa + 1) % NX, ya = ((y0 % NY) + NY) % NY, yb = (ya + 1) % NY;
        int za = std::clamp(z0, 0, NZ - 1), zb = std::clamp(z0 + 1, 0, NZ - 1);
        auto F = [&](int X, int Y, int Z) { return f[idx(X, Y, Z)]; };
        return (1 - fz) * ((1 - fy) * ((1 - fx) * F(xa, ya, za) + fx * F(xb, ya, za)) + fy * ((1 - fx) * F(xa, yb, za) + fx * F(xb, yb, za)))
             + fz * ((1 - fy) * ((1 - fx) * F(xa, ya, zb) + fx * F(xb, ya, zb)) + fy * ((1 - fx) * F(xa, yb, zb) + fx * F(xb, yb, zb)));
    }

    void writeSphere(std::ofstream& f, long s, double t, double dt, bool print, bool writeMap)
    {
        const double cx = (NX - 1) / 2.0, cy = (NY - 1) / 2.0, cz = (NZ - 1) / 2.0;
        const double rmax = std::min({NX, NY, NZ}) / 2.0 - 2.0, dr = 0.25;
        const int n = P.nDir, ns = (int)(rmax / dr);
        std::vector<float> rgb(n, NAN); std::vector<unsigned char> cov(n, 0), through(n, 0);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++)
        {
            double dPrev = 0.0, sPrev = 0.0; bool first = true; double gas = 0.0;
            for (int k = 0; k <= ns; k++)
            {
                double r = k * dr, x = cx + r * dirs[i][0], y = cy + r * dirs[i][1], z = cz + r * dirs[i][2];
                double a = interp(n1, x, y, z), b = interp(n2, x, y, z), cv = interp(c, x, y, z);
                if (cv >= 0.5) gas += dr;
                double d = b - a, sm = a + b;
                if (!first && std::isnan(rgb[i]) && dPrev > 0.0 && d <= 0.0 && sPrev > 0.5 && sm > 0.5)
                    rgb[i] = (float)(r - dr + dr * dPrev / (dPrev - d));
                dPrev = d; sPrev = sm; first = false;
            }
            if (gas >= P.minThick) cov[i] = 1;
            if (cov[i] && std::isnan(rgb[i])) through[i] = 1;
        }
        // components of covered directions
        std::vector<int> comp(n, -1); int ncomp = 0; long covered = 0, largest = 0, thr = 0;
        std::vector<int> q(n);
        for (int i0 = 0; i0 < n; i0++)
        {
            if (cov[i0]) covered++; if (through[i0]) thr++;
            if (!cov[i0] || comp[i0] >= 0) continue;
            int h = 0, tl = 0; q[tl++] = i0; comp[i0] = ncomp; long sz = 0;
            while (h < tl) { int u = q[h++]; sz++; for (int k = 0; k < 6; k++) { int v = dirNb[u][k]; if (v >= 0 && cov[v] && comp[v] < 0) { comp[v] = ncomp; q[tl++] = v; } } }
            ncomp++; largest = std::max(largest, sz);
        }
        double sr = 0, sr2 = 0, rmin = 1e30, rmx = -1e30; long ngb = 0;
        for (int i = 0; i < n; i++) if (!std::isnan(rgb[i])) { sr += rgb[i]; sr2 += (double)rgb[i] * rgb[i]; ngb++; rmin = std::min(rmin, (double)rgb[i]); rmx = std::max(rmx, (double)rgb[i]); }
        double rm = ngb ? sr / ngb : NAN, rs = ngb ? std::sqrt(std::max(0.0, sr2 / ngb - rm * rm)) : NAN;

        // volumes, centroid, gas, ranges
        long vin = 0, gcells = 0, bad = 0; double gx = 0, gy = 0, gz = 0, gr = 0, mass = 0, cmin = 1e30, cmax = -1e30, nmin = 1e30, nmax = -1e30;
        #pragma omp parallel for collapse(2) reduction(+:vin, gcells, bad, gx, gy, gz, gr, mass) reduction(min:cmin, nmin) reduction(max:cmax, nmax)
        for (int z = 0; z < NZ; z++)
            for (int y = 0; y < NY; y++)
                for (int x = 0; x < NX; x++)
                {
                    int p = idx(x, y, z); double a = n1[p], b = n2[p], cv = c[p];
                    if (b > a && a + b > 0.5) { vin++; gx += x; gy += y; gz += z; }
                    if (cv >= 0.5) { gcells++; gr += std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz)); }
                    mass += cv; cmin = std::min(cmin, cv); cmax = std::max(cmax, cv);
                    nmin = std::min({nmin, a, b}); nmax = std::max({nmax, a, b});
                    bad += !std::isfinite(cv) + !std::isfinite(a) + !std::isfinite(b);
                }
        if (bad) { std::cerr << "Field became non-finite at step " << s << ", t = " << t << "\n"; std::exit(2); }
        double rvol = std::cbrt(3.0 * vin / (4.0 * M_PI));
        double coff = vin ? std::sqrt((gx / vin - cx) * (gx / vin - cx) + (gy / vin - cy) * (gy / vin - cy) + (gz / vin - cz) * (gz / vin - cz)) : NAN;
        double drift = mass0 != 0.0 ? (mass - mass0) / mass0 : 0.0;
        double F = freeEnergy();
        if (s > 0 && F > lastF + 1e-9 * std::abs(lastF)) Fincreases++;
        lastF = F; maxAbsDrift = std::max(maxAbsDrift, std::abs(drift)); runCmin = std::min(runCmin, cmin); runCmax = std::max(runCmax, cmax);
        maxCount = std::max(maxCount, ncomp);
        if (vin == 0 && !innerVanished) { innerVanished = true; tVanish = t; }

        f << s << "," << std::setprecision(10) << t << "," << dt << "," << rm << "," << rs << "," << (ngb ? rmin : NAN) << "," << (ngb ? rmx : NAN) << ","
          << ngb << "," << double(covered) / n << "," << ncomp << "," << (covered ? double(largest) / covered : 0.0) << "," << thr << ","
          << vin << "," << rvol << "," << coff << "," << gcells << "," << (gcells ? gr / gcells : NAN) << ","
          << mass << "," << drift << "," << F << "," << cmin << "," << cmax << "," << nmin << "," << nmax << "," << bad << std::endl;
        if (print)
            std::cout << "t " << std::fixed << std::setprecision(2) << std::setw(7) << t << "  step " << std::setw(7) << s
                      << "  R_gb " << std::setprecision(2) << rm << "  R_vol " << rvol << "  cov " << std::setprecision(3) << double(covered) / n
                      << "  patches " << ncomp << "  A1/Atot " << (covered ? double(largest) / covered : 0.0)
                      << "  gas " << gcells << "  F " << std::setprecision(1) << F << std::defaultfloat << std::endl;
        if (writeMap)
        {
            std::filesystem::create_directories(path("maps"));
            std::ostringstream nm; nm << "maps/t" << std::fixed << std::setprecision(2) << std::setw(7) << std::setfill('0') << t << ".bin";
            std::ofstream g(path(nm.str()), std::ios::binary);
            std::vector<float> sc(NZ * NX), sn(NZ * NX); int y = (NY - 1) / 2;
            for (int z = 0; z < NZ; z++) for (int x = 0; x < NX; x++) { int p = idx(x, y, z); sc[z * NX + x] = (float)c[p]; sn[z * NX + x] = (float)(n2[p] - n1[p]); }
            g.write((const char*)sc.data(), sc.size() * sizeof(float));
            g.write((const char*)sn.data(), sn.size() * sizeof(float));
            g.write((const char*)rgb.data(), n * sizeof(float));
            g.write((const char*)cov.data(), n);
        }
    }

    // bubble-free GB diagnostics (gbOnly=1)
    void writeGB(std::ofstream& f, long s, double t, double dt, bool print)
    {
        double sum = 0.0, sum2 = 0.0, zmin = 1e30, zmax = -1e30, width = 0.0, smin = 1e30;
        long found = 0, multi = 0, none = 0, lowerN2 = 0;
        #pragma omp parallel for collapse(2) reduction(+:sum, sum2, found, multi, none, lowerN2, width) reduction(min:zmin, smin) reduction(max:zmax)
        for (int y = 0; y < NY; y++)
            for (int x = 0; x < NX; x++)
            {
                double zc = 0.0, del = 0.0; int cross = 0;
                if (n1[idx(x, y, 0)] < n2[idx(x, y, 0)]) lowerN2++;
                for (int z = 0; z < NZ; z++)
                {
                    int p = idx(x, y, z);
                    del += 16.0 * n1[p] * n1[p] * n2[p] * n2[p];
                    smin = std::min(smin, n1[p] * n1[p] + n2[p] * n2[p]);
                    if (z + 1 < NZ)
                    {
                        int q = idx(x, y, z + 1);
                        double da = n1[p] - n2[p], db = n1[q] - n2[q];
                        if ((da > 0.0 && db <= 0.0) || (da < 0.0 && db >= 0.0))
                        { if (++cross == 1) zc = z + da / (da - db); }
                    }
                }
                width += del;
                if (cross == 0) none++;
                else
                {
                    found++; if (cross > 1) multi++;
                    sum += zc; sum2 += zc * zc; zmin = std::min(zmin, zc); zmax = std::max(zmax, zc);
                }
            }
        double mass = 0.0, cmin = 1e30, cmax = -1e30, nmin = 1e30, nmax = -1e30; long bad = 0;
        #pragma omp parallel for reduction(+:mass, bad) reduction(min:cmin, nmin) reduction(max:cmax, nmax) schedule(static)
        for (int p = 0; p < N; p++)
        {
            mass += c[p]; cmin = std::min(cmin, c[p]); cmax = std::max(cmax, c[p]);
            nmin = std::min({nmin, n1[p], n2[p]}); nmax = std::max({nmax, n1[p], n2[p]});
            bad += !std::isfinite(c[p]) + !std::isfinite(n1[p]) + !std::isfinite(n2[p]);
        }
        if (bad) { std::cerr << "Field became non-finite at step " << s << ", t = " << t << "\n"; std::exit(2); }
        double F = freeEnergy();
        if (s > 0 && F > lastF + 1e-9 * std::abs(lastF)) Fincreases++;
        lastF = F;
        runCmin = std::min(runCmin, cmin); runCmax = std::max(runCmax, cmax);
        double zm = found ? sum / found : NAN;
        double zs = found ? std::sqrt(std::max(0.0, sum2 / found - zm * zm)) : NAN;
        f << s << "," << std::setprecision(12) << t << "," << dt << "," << zm << "," << zs << ","
          << (found ? zmin : NAN) << "," << (found ? zmax : NAN) << "," << found << "," << multi << "," << none << ","
          << lowerN2 << "," << width / (NX * NY) << "," << smin << "," << P.dfBias << "," << F << ","
          << mass << "," << cmin << "," << cmax << "," << nmin << "," << nmax << "," << bad << std::endl;
        if (print)
            std::cout << "t " << std::fixed << std::setprecision(2) << std::setw(7) << t << "  step " << std::setw(8) << s
                      << "  GB z mean " << std::setprecision(4) << zm << "  std " << std::scientific << std::setprecision(2) << zs
                      << "  F " << std::setprecision(8) << F << "  max|c| " << std::setprecision(1) << std::max(std::abs(cmin), std::abs(cmax))
                      << std::defaultfloat << std::endl;
    }

    void saveMorphology(long s, double t, std::ofstream& index)
    {
        long snapNo = std::lround(t / P.snapDt);
        bool raw = P.rawEvery > 0 && snapNo % P.rawEvery == 0;
        std::string name = (P.gbOnly || !raw) ? std::string("") : "morphology_step_" + std::to_string(s) + ".raw";
        if (!P.gbOnly && raw)
        {
            std::ofstream f(path(name), std::ios::binary);
            std::vector<float> buf(c.begin(), c.end());
            f.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(float));
        }
        std::string gname;
        if ((P.saveGrains && raw) || P.gbOnly)
        {
            gname = "grains_step_" + std::to_string(s) + ".raw";
            std::ofstream g(path(gname), std::ios::binary);
            std::vector<float> b1(n1.begin(), n1.end()), b2(n2.begin(), n2.end());
            g.write(reinterpret_cast<const char*>(b1.data()), b1.size() * sizeof(float));
            g.write(reinterpret_cast<const char*>(b2.data()), b2.size() * sizeof(float));
        }
        index << std::setprecision(10) << t << "," << s << "," << name << "," << gname << std::endl;
    }

    void run()
    {
        std::filesystem::create_directories(P.outDir);
        double t = 0.0;
        long s = 0;
        const bool restart = !P.restartIn.empty();
        if (restart) loadCheckpoint(P.restartIn, t, s);
        else initialize();
        if (P.initOnly) { std::cout << "Initial condition written to " << P.initOut << "\n"; return; }

        {
            std::ofstream pf(path(restart ? "params_restart_t" + std::to_string(std::lround(t)) + ".txt" : "params.txt"));
            pf << P.describe();
            pf << "initial gas total    = " << std::setprecision(17) << mass0 << "\n";
            pf << "bubble centres (x y):\n" << std::setprecision(17);
            for (auto& q : centres) pf << "  " << q.first << " " << q.second << "\n";
        }
        if (P.bGB <= 1.0)
            std::cout << "WARNING: bGB = " << P.bGB << " <= 1. The grain-boundary potential has no barrier\n"
                         "         (bGB = 1) or the grains are unstable (bGB < 1); the GB will spread.\n";

        const bool SPH = P.sphereR > 0.0;
        std::ofstream f(path(SPH ? "sphere.csv" : P.gbOnly ? "gb_position.csv" : "bubble_metrics.csv"), restart ? std::ios::app : std::ios::trunc);
        if (SPH) setupDirs();
        if (!restart && SPH) f << "step,time,dt,r_gb_mean,r_gb_std,r_gb_min,r_gb_max,n_dir_gb,coverage,n_patches,largest_fraction,"
             "dirs_gb_through_gas,V_inner,R_vol,centroid_offset,gas_cells,gas_r_mean,gas_total,mass_drift,free_energy,"
             "c_min,c_max,n_min,n_max,nonfinite\n";
        else if (!restart && P.gbOnly) f << "step,time,dt,gb_z_mean,gb_z_std,gb_z_min,gb_z_max,columns_with_gb,"
             "columns_multiple_crossings,columns_no_crossing,columns_lower_grain_n2,gb_width_int16n1n2n2,"
             "min_n1sq_plus_n2sq,df_bias,free_energy,gas_total,c_min,c_max,n_min,n_max,nonfinite\n";
        else if (!restart) f << "step,time_dimensionless,time_seconds,time_hours,number_of_projected_bubbles,"
             "total_projected_area,largest_projected_area,largest_area_fraction,coverage,porosity,"
             "dt,gas_total,mass_drift,c_min,c_max,"
             "free_energy,gb_z_mean,gb_z_std,gb_delta,n_min,n_max,nonfinite\n";
        std::ofstream cf;
        if (!P.gbOnly && !SPH)
        {
            cf.open(path("gb_bubble.csv"), restart ? std::ios::app : std::ios::trunc);
            if (!restart) cf << "time,gb_real_mean,gb_real_std,gb_real_min,gb_real_max,cols_real_gb,cols_multi,"
                                "cols_covered,cols_gb_through_bubble,cols_gb_above_bubble,cols_gb_below_bubble,"
                                "gb_coverage,bub_z_mean,bub_zbot_mean,bub_ztop_mean,vol_II\n";
        }
        std::ofstream index(path("snapshots.csv"), restart ? std::ios::app : std::ios::trunc);
        if (!restart) index << "time_dimensionless,step,c_file,grains_file\n";

        std::cout << P.describe()
                  << "Gas: conserved Cahn-Hilliard, no source, no clamp. Initial total gas = "
                  << std::setprecision(17) << mass0 << std::defaultfloat << "\n";

        if (!restart)
        {
            if (SPH) writeSphere(f, 0, 0.0, 0.0, true, true);
            else if (P.gbOnly) writeGB(f, 0, 0.0, 0.0, true); else { writeMetrics(f, 0, 0.0, 0.0, true); writeCoupled(cf, 0.0, true); }
            saveMorphology(0, 0.0, index);
        }

        const long nOut = std::lround(P.tEnd / P.outDt);
        const long snapEvery = std::max(1L, std::lround(P.snapDt / P.outDt));
        const long printEvery = std::max(1L, std::lround(1.0 / P.outDt));

        const long kStart = std::lround(t / P.outDt) + 1;
        double dtSum = 0.0, dtMin = 1e30, dtMaxSeen = 0.0; long dtCount = 0;
        double dtCHmin = 1e30, dtACmin = 1e30;
        for (long k = kStart; k <= nOut; k++)
        {
            const double tOut = k * P.outDt;
            double dtReport = 0.0;
            while (t < tOut)
            {
                double left = tOut - t;
                if (left <= 1e-9 * P.outDt) { t = tOut; break; }
                double dt = step(left);
                s++;
                dtCHmin = std::min(dtCHmin, lastDtCH); dtACmin = std::min(dtACmin, lastDtAC);
                if (dt >= left) t = tOut;
                else { t += dt; dtReport = dt; dtSum += dt; dtMin = std::min(dtMin, dt);
                       dtMaxSeen = std::max(dtMaxSeen, dt); dtCount++; }
            }
            if (SPH) writeSphere(f, s, t, dtReport, k % printEvery == 0, k % snapEvery == 0);
            else if (P.gbOnly) writeGB(f, s, t, dtReport, k % printEvery == 0);
            else { writeMetrics(f, s, t, dtReport, k % printEvery == 0); writeCoupled(cf, t, k % snapEvery == 0); }
            if (k % snapEvery == 0) { saveMorphology(s, t, index); if (P.ckpt) saveCheckpoint(t, s); }
            if (SPH && innerVanished && P.stopAfterVanish >= 0.0 && t >= tVanish + P.stopAfterVanish - 1e-9)
            { std::cout << "Inner grain vanished at t = " << tVanish << "; stopping at t = " << t << "\n"; break; }
        }
        if (P.ckpt) saveCheckpoint(t, s);

        std::ostringstream d;
        d << std::scientific << std::setprecision(4)
          << "steps                       = " << s << "\n"
          << "dt mean / min / max         = " << (dtCount ? dtSum / dtCount : 0.0) << " / " << dtMin << " / " << dtMaxSeen << "\n"
          << "smallest dt_CH, dt_AC limit = " << dtCHmin << ", " << dtACmin << "  (safety " << P.safety << ")\n"
          << "max |mass drift|            = " << maxAbsDrift << "\n"
          << "c range over run            = [" << runCmin << ", " << runCmax << "]\n"
          << "free-energy increases       = " << Fincreases << " (of " << nOut << " outputs)\n"
          << "max projected bubble count  = " << maxCount << (P.gbOnly ? "  (gbOnly: not computed)" : "") << "\n"
          << "dfBias, gbOnly              = " << P.dfBias << ", " << P.gbOnly << "\n";
        std::cout << "Done.\n" << d.str();
        std::ofstream df(path("diagnostics.txt"));
        df << "end time                    = " << t << (restart ? "  (restarted run; statistics cover the whole run)" : "") << "\n";
        df << d.str();
    }

private:
    Params P;
    std::vector<double> c, n1, n2, cN, n1N, n2N, mu, mob, lapN1, lapN2;
    std::vector<std::pair<double, double>> centres;
    std::vector<std::array<double, 3>> centres3d;
    int xm[NX], xp[NX], ym[NY], yp[NY], zm[NZ], zp[NZ];
    std::mt19937 rng;
    double mass0 = 0.0;
    double lastDtCH = 0.0, lastDtAC = 0.0;
    double lastF = 0.0, maxAbsDrift = 0.0, runCmin = 1e30, runCmax = -1e30;
    long Fincreases = 0; int maxCount = 0;
};

int main(int argc, char** argv)
{
    Params p;
    for (int i = 1; i < argc; i++)
    {
        std::string a = argv[i];
        auto eq = a.find('=');
        if (eq == std::string::npos) continue;
        std::string k = a.substr(0, eq);
        std::string sv = a.substr(eq + 1);
        double v = std::atof(sv.c_str());
        if (k == "theta") p.setTheta((int)v);
        else if (k == "A0") p.A0 = v;
        else if (k == "bS") p.bS = v;
        else if (k == "bGB") p.bGB = v;
        else if (k == "j") p.j = v;
        else if (k == "jGB") p.jGB = v;
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
        else if (k == "noise") p.noise = v;
        else if (k == "r0") p.r0 = v;
        else if (k == "nBubbles") p.nBubbles = (int)v;
        else if (k == "minSep") p.minSep = v;
        else if (k == "outDir") p.outDir = sv;
        else if (k == "initIn") p.initIn = sv;
        else if (k == "initOut") p.initOut = sv;
        else if (k == "initOnly") p.initOnly = v != 0.0;
        else if (k == "centresIn") p.centresIn = sv;
        else if (k == "saveGrains") p.saveGrains = v != 0.0;
        else if (k == "restartIn") p.restartIn = sv;
        else if (k == "ckpt") p.ckpt = v != 0.0;
        else if (k == "dfBias") p.dfBias = v;
        else if (k == "gbOnly") p.gbOnly = v != 0.0;
        else if (k == "zgb0") p.zgb0 = v;
        else if (k == "biasMode") p.biasMode = (int)v;
        else if (k == "rawEvery") p.rawEvery = (int)v;
        else if (k == "sliceY") p.sliceY = (int)v;
        else if (k == "sphereR") p.sphereR = v;
        else if (k == "nDir") p.nDir = (int)v;
        else if (k == "stopAfterVanish") p.stopAfterVanish = v;
        else { std::cerr << "unknown parameter " << k << "\n"; return 1; }
    }
    Simulation sim(p);
    sim.run();
    return 0;
}
