// ============================================================
// 2D Phase-Field Simulation of Gas Bubbles on a Grain Boundary
//
// c          : conserved gas concentration (Cahn-Hilliard)
// phi1, phi2 : non-conserved grain order parameters (Allen-Cahn)
//
// Free energy:
//
// F = int [ f(c,phi)
//         + kap/2 |grad c|^2
//         + kap/2 |grad phi1|^2
//         + kap/2 |grad phi2|^2 ] dA
//
// f = A0 [ c^4/4 - c^2/2
//        + phi1^4/4 - phi1^2/2
//        + phi2^4/4 - phi2^2/2
//        + bGB phi1^2 phi2^2
//        + bS c^2 (phi1^2 + phi2^2) ]
//
// dc/dt     = div( M grad(mu) )
// mu        = dF/dc
//
// dphi/dt   = -Lphi dF/dphi
//
// Geometry:
//
//      Grain 2
//      phi2 = 1
//      --------------------
//          GRAIN BOUNDARY
//      --------------------
//      Grain 1
//      phi1 = 1
//
// The GB is the horizontal line y = Ny/2.
//
// Periodic boundary conditions in x.
// No-flux boundary conditions in y.
//
// Gas bubbles nucleate near the GB and grow/coalesce.
//
// Output:
//
// bubble_data.dat
//
// Columns:
//
// step
// time
// coverage
// Alargest_over_Atotal
// Atotal
// Nbubbles
//
// Build:
//
// g++ -O3 -fopenmp -std=c++17 phasefield_bubbles_2d.cpp -o pf2d
//
// Run:
//
// ./pf2d
//
// or:
//
// ./pf2d Nx Ny nsteps out_every
//
// Example:
//
// ./pf2d 128 128 60000 200
//
// ============================================================

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

using std::vector;


// ============================================================
// PARAMETERS
// ============================================================

static double A0   = 4.0;

static double kap  = 2.0;

static double bGB  = 1.0;

static double bS   = 0.8;

static double Dg   = 1.0;

static double Lphi = 1.0;

static double dt   = 0.005;


// ------------------------------------------------------------
// Gas concentration
// ------------------------------------------------------------

static double c_bulk  = 0.33;

static double c_layer = 0.35;

static double noise   = 0.05;


// Width of noisy region around GB

static int layerHalfWidth = 2;


// Bubble threshold

static double cthr = 0.5;


// ============================================================
// GRID
// ============================================================

int Nx = 128;
int Ny = 128;


// ============================================================
// INDEXING
// ============================================================

inline int idx(int i, int j)
{
    return j * Nx + i;
}


// ============================================================
// PERIODIC X
// ============================================================

inline int wrapx(int i)
{
    if (i < 0)
        return i + Nx;

    if (i >= Nx)
        return i - Nx;

    return i;
}


// ============================================================
// CLAMP Y
//
// Used for no-flux boundary condition.
//
// At y = 0:
//     y-1 -> y=0
//
// At y = Ny-1:
//     y+1 -> y=Ny-1
// ============================================================

inline int clampy(int j)
{
    if (j < 0)
        return 0;

    if (j >= Ny)
        return Ny - 1;

    return j;
}


// ============================================================
// 2D 5-POINT LAPLACIAN
//
// dx = dy = 1
//
// Periodic in x
// No-flux in y
// ============================================================

inline double lap(
    const vector<double>& f,
    int i,
    int j)
{
    double c0 = f[idx(i, j)];

    return
        f[idx(wrapx(i + 1), j)] +
        f[idx(wrapx(i - 1), j)] +
        f[idx(i, clampy(j + 1))] +
        f[idx(i, clampy(j - 1))] -
        4.0 * c0;
}


// ============================================================
// BUBBLE STATISTICS
// ============================================================

struct Stats
{
    double coverage;

    double ratio;

    double atot;

    int nb;
};


// ============================================================
// 2D BUBBLE ANALYSIS
//
// Bubble = connected region where
//
//     c > cthr
//
// Connectivity:
//
//     4-neighbour
//
// Periodic in x.
// No periodicity in y.
//
// Since the simulation is now 2D, every bubble has an AREA,
// not a volume.
// ============================================================

Stats analyse(const vector<double>& c)
{
    const int N = Nx * Ny;

    vector<int> label(N, -1);

    vector<long> areas;

    vector<int> stack;

    long atot = 0;


    // --------------------------------------------------------
    // Find connected bubbles
    // --------------------------------------------------------

    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            int id = idx(i, j);


            // Not bubble or already visited

            if (c[id] <= cthr || label[id] >= 0)
                continue;


            int lab = (int)areas.size();

            areas.push_back(0);

            stack.clear();

            stack.push_back(id);

            label[id] = lab;


            // ------------------------------------------------
            // Flood fill
            // ------------------------------------------------

            while (!stack.empty())
            {
                int cur = stack.back();

                stack.pop_back();


                int ci = cur % Nx;

                int cj = cur / Nx;


                areas[lab]++;


                // ------------------------------------------------
                // Four neighbours
                // ------------------------------------------------

                int ni[4] =
                {
                    wrapx(ci + 1),
                    wrapx(ci - 1),
                    ci,
                    ci
                };


                int nj[4] =
                {
                    cj,
                    cj,
                    clampy(cj + 1),
                    clampy(cj - 1)
                };


                for (int n = 0; n < 4; n++)
                {
                    int ii = ni[n];

                    int jj = nj[n];

                    int q = idx(ii, jj);


                    if (label[q] >= 0)
                        continue;


                    if (c[q] <= cthr)
                        continue;


                    label[q] = lab;

                    stack.push_back(q);
                }
            }
        }
    }


    // --------------------------------------------------------
    // Total bubble area
    // --------------------------------------------------------

    for (long a : areas)
        atot += a;


    // --------------------------------------------------------
    // Largest bubble area
    // --------------------------------------------------------

    long amax = 0;

    if (!areas.empty())
    {
        amax =
            *std::max_element(
                areas.begin(),
                areas.end()
            );
    }


    // --------------------------------------------------------
    // Coverage
    //
    // fraction of entire 2D domain occupied by bubbles
    //
    // coverage = Atotal / (Nx * Ny)
    // --------------------------------------------------------

    double coverage =
        (double)atot /
        (double)(Nx * Ny);


    // --------------------------------------------------------
    // Largest bubble / total bubble area
    // --------------------------------------------------------

    double ratio =
        (atot > 0)
        ? (double)amax / (double)atot
        : 0.0;


    Stats s;

    s.coverage = coverage;

    s.ratio = ratio;

    s.atot = (double)atot;

    s.nb = (int)areas.size();


    return s;
}


// ============================================================
// MAIN
// ============================================================

int main(int argc, char** argv)
{
    // --------------------------------------------------------
    // DEFAULT SIMULATION
    // --------------------------------------------------------

    long nsteps = 60000;

    long outEvery = 200;


    // --------------------------------------------------------
    // COMMAND LINE
    //
    // ./pf2d Nx Ny nsteps out_every
    //
    // Example:
    //
    // ./pf2d 128 128 60000 200
    // --------------------------------------------------------

    if (argc > 2)
    {
        Nx = atoi(argv[1]);

        Ny = atoi(argv[2]);
    }


    if (argc > 3)
        nsteps = atol(argv[3]);


    if (argc > 4)
        outEvery = atol(argv[4]);


    // ========================================================
    // ENVIRONMENT VARIABLE OVERRIDES
    // ========================================================

    auto envd =
        [](const char* name, double& value)
        {
            if (const char* s = getenv(name))
                value = atof(s);
        };


    envd("BS", bS);

    envd("BGB", bGB);

    envd("KAP", kap);

    envd("DT", dt);

    envd("CBULK", c_bulk);

    envd("CLAYER", c_layer);

    envd("NOISE", noise);

    envd("CTHR", cthr);

    envd("DG", Dg);


    if (const char* s = getenv("LAYER"))
        layerHalfWidth = atoi(s);


    // ========================================================
    // ALLOCATE FIELDS
    // ========================================================

    const int N = Nx * Ny;


    vector<double> c(N, 0.0);

    vector<double> p1(N, 0.0);

    vector<double> p2(N, 0.0);

    vector<double> mu(N, 0.0);


    vector<double> cn(N, 0.0);

    vector<double> p1n(N, 0.0);

    vector<double> p2n(N, 0.0);


    // ========================================================
    // RANDOM NUMBER GENERATOR
    // ========================================================

    std::mt19937 rng(12345);

    std::uniform_real_distribution<double> U(-1.0, 1.0);

    std::uniform_real_distribution<double> U01(0.0, 1.0);


    // ========================================================
    // GRAIN BOUNDARY
    //
    // Grain 1:
    //
    //     j < kGB
    //
    // Grain 2:
    //
    //     j >= kGB
    // ========================================================

    int kGB = Ny / 2;


    // ========================================================
    // BUBBLE SEEDS
    // ========================================================

    int nseed = 60;

    double r0 = 2.0;


    if (const char* s = getenv("NSEED"))
        nseed = atoi(s);


    envd("R0", r0);


    vector<double> sx(nseed);

    vector<double> sy(nseed);

    vector<double> sr(nseed);


    // --------------------------------------------------------
    // Random seed positions
    //
    // Seeds are centered around the GB.
    // --------------------------------------------------------

    for (int s = 0; s < nseed; s++)
    {
        sx[s] =
            U01(rng) * Nx;


        // Place seed close to GB

        sy[s] =
            kGB +
            (U01(rng) - 0.5) *
            2.0 *
            layerHalfWidth;


        sr[s] =
            r0 *
            (0.7 + 0.6 * U01(rng));
    }


    // ========================================================
    // INITIAL CONDITION
    // ========================================================

    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            int id = idx(i, j);


            // ------------------------------------------------
            // Grain 1
            // ------------------------------------------------

            if (j < kGB)
            {
                p1[id] = 1.0;

                p2[id] = 0.0;
            }


            // ------------------------------------------------
            // Grain 2
            // ------------------------------------------------

            else
            {
                p1[id] = 0.0;

                p2[id] = 1.0;
            }


            // ------------------------------------------------
            // Background gas concentration
            // ------------------------------------------------

            double cc = c_bulk;


            // ------------------------------------------------
            // Initial circular bubble seeds
            // ------------------------------------------------

            for (int s = 0; s < nseed; s++)
            {
                double dx =
                    std::fabs(
                        i - sx[s]
                    );


                // Periodic x

                dx =
                    std::min(
                        dx,
                        Nx - dx
                    );


                double dy =
                    std::fabs(
                        j - sy[s]
                    );


                double r =
                    std::sqrt(
                        dx * dx +
                        dy * dy
                    );


                // Smooth circular bubble profile

                double prof =
                    0.5 *
                    (
                        1.0 -
                        std::tanh(
                            (r - sr[s]) / 1.5
                        )
                    );


                cc =
                    std::max(
                        cc,
                        prof
                    );
            }


            // ------------------------------------------------
            // Noise near grain boundary
            // ------------------------------------------------

            if (std::abs(j - kGB)
                <= layerHalfWidth)
            {
                cc +=
                    noise *
                    U(rng);
            }


            // ------------------------------------------------
            // Clamp concentration
            // ------------------------------------------------

            c[id] =
                std::min(
                    1.0,
                    std::max(
                        0.0,
                        cc
                    )
                );
        }
    }


    // ========================================================
    // OUTPUT FILE
    // ========================================================

    FILE* fo =
        fopen(
            "bubble_data.dat",
            "w"
        );


    if (!fo)
    {
        printf(
            "ERROR: could not open bubble_data.dat\n"
        );

        return 1;
    }


    fprintf(
        fo,
        "# step time coverage "
        "Alargest_over_Atotal "
        "Atotal Nbubbles\n"
    );


    // ========================================================
    // CH MOBILITY
    // ========================================================

    const double M =
        Dg /
        (
            A0 *
            (2.0 * bS - 1.0)
        );


    // ========================================================
    // TIME LOOP
    // ========================================================

    for (long step = 0;
         step <= nsteps;
         step++)
    {
        // ----------------------------------------------------
        // OUTPUT
        // ----------------------------------------------------

        if (step % outEvery == 0)
        {
            Stats s =
                analyse(c);


            fprintf(
                fo,
                "%ld %.4f %.6f %.6f %.0f %d\n",

                step,

                step * dt,

                s.coverage,

                s.ratio,

                s.atot,

                s.nb
            );


            fflush(fo);


            // Print every 10th output

            if ((step / outEvery) % 10 == 0)
            {
                printf(
                    "step %7ld  "
                    "time %.4f  "
                    "cov %.4f  "
                    "Alargest/Atot %.4f  "
                    "Atot %.0f  "
                    "nb %d\n",

                    step,

                    step * dt,

                    s.coverage,

                    s.ratio,

                    s.atot,

                    s.nb
                );
            }
        }


        // ----------------------------------------------------
        // Stop after final output
        // ----------------------------------------------------

        if (step == nsteps)
            break;


        // ====================================================
        // CHEMICAL POTENTIAL
        //
        // mu = dF/dc
        //
        // mu =
        //
        // A0 [c^3 - c
        //     + 2 bS c (phi1^2 + phi2^2)]
        //
        // - kap lap(c)
        // ====================================================

        #pragma omp parallel for

        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                int id =
                    idx(i, j);


                double cc =
                    c[id];


                double s2 =
                    p1[id] * p1[id] +
                    p2[id] * p2[id];


                mu[id] =
                    A0 *
                    (
                        cc * cc * cc
                        -
                        cc
                        +
                        2.0 *
                        bS *
                        cc *
                        s2
                    )
                    -
                    kap *
                    lap(
                        c,
                        i,
                        j
                    );
            }
        }


        // ====================================================
        // CAHN-HILLIARD + ALLEN-CAHN
        // ====================================================

        #pragma omp parallel for

        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                int id =
                    idx(i, j);


                double cc =
                    c[id];


                double a =
                    p1[id];


                double b =
                    p2[id];


                double c2 =
                    cc * cc;


                // ------------------------------------------------
                // dF/dphi1
                // ------------------------------------------------

                double dF1 =
                    A0 *
                    (
                        a * a * a
                        -
                        a
                        +
                        2.0 *
                        bGB *
                        a *
                        b *
                        b
                        +
                        2.0 *
                        bS *
                        c2 *
                        a
                    )
                    -
                    kap *
                    lap(
                        p1,
                        i,
                        j
                    );


                // ------------------------------------------------
                // dF/dphi2
                // ------------------------------------------------

                double dF2 =
                    A0 *
                    (
                        b * b * b
                        -
                        b
                        +
                        2.0 *
                        bGB *
                        b *
                        a *
                        a
                        +
                        2.0 *
                        bS *
                        c2 *
                        b
                    )
                    -
                    kap *
                    lap(
                        p2,
                        i,
                        j
                    );


                // ------------------------------------------------
                // Cahn-Hilliard
                // ------------------------------------------------

                cn[id] =
                    cc +
                    dt *
                    M *
                    lap(
                        mu,
                        i,
                        j
                    );


                // ------------------------------------------------
                // Allen-Cahn
                // ------------------------------------------------

                p1n[id] =
                    a -
                    dt *
                    Lphi *
                    dF1;


                p2n[id] =
                    b -
                    dt *
                    Lphi *
                    dF2;
            }
        }


        // ====================================================
        // SWAP
        // ====================================================

        c.swap(cn);

        p1.swap(p1n);

        p2.swap(p2n);
    }


    // ========================================================
    // CLOSE
    // ========================================================

    fclose(fo);


    printf(
        "done. wrote bubble_data.dat\n"
    );


    return 0;
}