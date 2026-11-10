
/*
    phasefield_gb_bubbles3d.cpp

    3-D phase-field simulation of gas bubbles on a flat grain boundary.

    Geometry:
        Total domain = 128 x 128 x 64
        Grain 1     = z = 0 ... 31
        Grain 2     = z = 32 ... 63
        Grain boundary is the plane around z = 31.5

    The two grain order parameters are:
        n1 = grain 1 order parameter
        n2 = grain 2 order parameter

    Gas phase field:
        c = gas concentration / bubble phase field

    The model uses the equations discussed from the paper:

        F = integral [
                A0 * f(c,n1,n2)
                + j/2 * |grad c|^2
                + jGB/2 * (|grad n1|^2 + |grad n2|^2)
            ] dV

    f(c,n1,n2) =
        (c^4/4 - c^2/2)
        + (n1^4/4 - n1^2/2)
        + (n2^4/4 - n2^2/2)
        + bGB/2 * n1^2 * n2^2
        + bS * c^2 * (n1^2 + n2^2)

    Gas chemical potential:

        mu = dF/dc

            = A0 * [
                  c^3 - c
                  + 2*bS*c*(n1^2+n2^2)
              ]
              - j * laplacian(c)

    Gas evolution:

        dc/dt = div( M * grad(mu) )

    IMPORTANT:
    M is NOT constant here. We use the mobility form requested:

        M(c,n1,n2)
          = Ms  * 16*c^2*(1-c^2)
          + Mgb * n1^2*n2^2
          + Mb

    This makes:
        - bulk mobility approximately Mb
        - GB mobility high where n1 and n2 overlap
        - bubble-interface mobility high where c is between 0 and 1

    Grain evolution:

        dn1/dt = -LGB * dF/dn1
        dn2/dt = -LGB * dF/dn2

    where

        dF/dn1 =
            A0 * [
                n1^3 - n1
                + bGB*n1*n2^2
                + 2*bS*c^2*n1
            ]
            - jGB*laplacian(n1)

        dF/dn2 =
            A0 * [
                n2^3 - n2
                + bGB*n2*n1^2
                + 2*bS*c^2*n2
            ]
            - jGB*laplacian(n2)

    Contact-angle parameter sets from the paper:

        theta = 40 deg : bS = 1.0, bGB = 1.6
        theta = 50 deg : bS = 1.0, bGB = 1.3
        theta = 60 deg : bS = 1.0, bGB = 1.0
        theta = 80 deg : bS = 1.7, bGB = 0.7

    This file uses THETA = 60 by default.
    Change THETA to 40, 50, 60 or 80 and compile again.

    Physical scaling:
        dx = 10 nm
        Dg from the paper = 3.4e-18 m^2/s

    Initial condition:
        33 bubbles
        radius = 44 nm
        initial matrix concentration = 0.12
        bubble concentration = 1.0
        small random noise is added ONLY at t=0
        random seed = 2024

    Output:
        bubble_metrics.csv
        morphology_step_XXXXXX.raw

    The RAW morphology file contains float values of c with shape:
        [NX, NY, NZ]

    Compile:
        g++ -O3 -fopenmp -std=c++17 phasefield_gb_bubbles3d.cpp -o phasefield_gb_bubbles3d

    Run:
        ./phasefield_gb_bubbles3d
*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

class PhaseField3D
{
public:

    // ------------------------------------------------------------
    // Grid
    // ------------------------------------------------------------

    static constexpr int NX = 128;
    static constexpr int NY = 128;
    static constexpr int NZ = 64;

    static constexpr double dx = 1.0;

    // Physical length represented by one grid spacing.
    static constexpr double dxPhysical_m = 10.0e-9;

    // ------------------------------------------------------------
    // Initial bubbles
    // ------------------------------------------------------------

    static constexpr int numberOfBubbles = 33;

    // 44 nm / 10 nm = 4.4 grid cells.
    static constexpr double initialRadius = 4.4;

    // 160 nm / 10 nm = 16 grid cells.
    // This follows the minimum separation used in the paper.
    static constexpr double minimumBubbleSeparation = 16.0;

    // ------------------------------------------------------------
    // Gas concentration
    // ------------------------------------------------------------

    static constexpr double cInitial = 0.12;
    static constexpr double cBubble = 1.0;

    // Small random perturbation.
    // It is used only when creating the initial field.
    static constexpr double noiseAmplitude = 0.001;

    static constexpr unsigned int randomSeed = 2024;

    // ------------------------------------------------------------
    // Free-energy parameters
    // ------------------------------------------------------------

    static constexpr double A0 = 4.0;

    // Default contact angle.
    // Change to 40, 50, 60 or 80.
    static constexpr int THETA = 60;

    /*
        Contact-angle parameters from the paper.

        theta = 40 : bS = 1.0, bGB = 1.6
        theta = 50 : bS = 1.0, bGB = 1.3
        theta = 60 : bS = 1.0, bGB = 1.0
        theta = 80 : bS = 1.7, bGB = 0.7
    */

    static constexpr double bS =
        (THETA == 80) ? 1.7 : 1.0;

    static constexpr double bGB =
        (THETA == 40) ? 1.6 :
        (THETA == 50) ? 1.3 :
        (THETA == 60) ? 1.0 :
        (THETA == 80) ? 0.7 :
        1.0;

    // Interface gradient coefficients.
    //
    // The paper keeps approximately:
    //     jS / bS = 2
    //
    // We use j = 2*bS for the gas interface.
    // For the GB, jGB = 2*bGB.
    static constexpr double j = 2.0 * bS;
    static constexpr double jGB = 2.0 * bGB;

    // ------------------------------------------------------------
    // Allen-Cahn parameters
    // ------------------------------------------------------------

    // Grain boundary relaxation rate.
    // This is a numerical parameter, not a value taken directly
    // from the paper.
    static constexpr double LGB = 0.10;

    // ------------------------------------------------------------
    // Mobility parameters
    // ------------------------------------------------------------

    /*
        Requested mobility:

        M = Ms*16*c^2*(1-c^2)
            + Mgb*n1^2*n2^2
            + Mb

        These are dimensionless model parameters.

        Mb  : bulk mobility
        Mgb : enhanced GB mobility
        Ms  : enhanced bubble-surface mobility
    */

    static constexpr double Mb = 0.01;
    static constexpr double Mgb = 2.0;
    static constexpr double Ms = 0.50;

    // ------------------------------------------------------------
    // Time stepping
    // ------------------------------------------------------------

    /*
        The CH equation contains fourth-order spatial derivatives.

        Therefore an explicit finite-difference calculation needs
        a small timestep.

        This is a numerical parameter.
    */
    static constexpr double dt = 0.00005;

    static constexpr int totalSteps = 30000;

    // Save metrics and morphology every this many steps.
    static constexpr int outputEvery = 1000;

    // ------------------------------------------------------------
    // Physical diffusivity used only for time conversion.
    // The paper gives approximately 3.4e-18 m^2/s.
    // ------------------------------------------------------------

    static constexpr double DPhysical = 3.4e-18;

private:

    // ------------------------------------------------------------
    // Arrays
    // ------------------------------------------------------------

    std::vector<double> c;
    std::vector<double> n1;
    std::vector<double> n2;

    std::vector<double> cNew;
    std::vector<double> n1New;
    std::vector<double> n2New;

    std::vector<double> mu;

    std::vector<double> lapC;
    std::vector<double> lapN1;
    std::vector<double> lapN2;

    // Fluxes for variable-mobility Cahn-Hilliard.
    //
    // Fx[i,j,k] is the flux through the x-face between
    // (i,j,k) and (i+1,j,k).
    std::vector<double> fluxX;
    std::vector<double> fluxY;
    std::vector<double> fluxZ;

    std::mt19937 generator;

public:

    PhaseField3D()
        : c(NX * NY * NZ, 0.0),
          n1(NX * NY * NZ, 0.0),
          n2(NX * NY * NZ, 0.0),
          cNew(NX * NY * NZ, 0.0),
          n1New(NX * NY * NZ, 0.0),
          n2New(NX * NY * NZ, 0.0),
          mu(NX * NY * NZ, 0.0),
          lapC(NX * NY * NZ, 0.0),
          lapN1(NX * NY * NZ, 0.0),
          lapN2(NX * NY * NZ, 0.0),
          fluxX(NX * NY * NZ, 0.0),
          fluxY(NX * NY * NZ, 0.0),
          fluxZ(NX * NY * NZ, 0.0),
          generator(randomSeed)
    {
    }

    // ------------------------------------------------------------
    // Index
    // ------------------------------------------------------------

    inline int index(int x, int y, int z) const
    {
        return (z * NY + y) * NX + x;
    }

    // ------------------------------------------------------------
    // Periodic x and y.
    // z is NOT periodic because the two grains have physical
    // top and bottom surfaces.
    // ------------------------------------------------------------

    inline int periodicX(int x) const
    {
        if (x < 0)
            return x + NX;

        if (x >= NX)
            return x - NX;

        return x;
    }

    inline int periodicY(int y) const
    {
        if (y < 0)
            return y + NY;

        if (y >= NY)
            return y - NY;

        return y;
    }

    // ------------------------------------------------------------
    // Get a value with:
    // periodic x/y
    // no-flux boundary in z
    // ------------------------------------------------------------

    double getValue(
        const std::vector<double>& field,
        int x,
        int y,
        int z) const
    {
        x = periodicX(x);
        y = periodicY(y);

        if (z < 0)
            z = 0;

        if (z >= NZ)
            z = NZ - 1;

        return field[index(x, y, z)];
    }

    // ------------------------------------------------------------
    // 3-D Laplacian
    // ------------------------------------------------------------

    double laplacian(
        const std::vector<double>& field,
        int x,
        int y,
        int z) const
    {
        double center = getValue(field, x, y, z);

        double left = getValue(field, x - 1, y, z);
        double right = getValue(field, x + 1, y, z);

        double down = getValue(field, x, y - 1, z);
        double up = getValue(field, x, y + 1, z);

        double back = getValue(field, x, y, z - 1);
        double front = getValue(field, x, y, z + 1);

        return (
            left + right +
            down + up +
            back + front -
            6.0 * center
        ) / (dx * dx);
    }

    // ------------------------------------------------------------
    // Distance from the GB.
    //
    // The GB is around z = 31.5.
    // ------------------------------------------------------------

    double distanceFromGB(int z) const
    {
        return std::abs(static_cast<double>(z) - 31.5);
    }

    // ------------------------------------------------------------
    // Initial grain order parameters
    // ------------------------------------------------------------

    void initializeGrains()
    {
        /*
            We create two grains separated by a diffuse interface.

            Grain 1:
                n1 ~ 1
                n2 ~ 0

            Grain 2:
                n1 ~ 0
                n2 ~ 1

            At the GB:
                n1 and n2 both have intermediate values.
        */

        double interfaceWidth = 1.5;

        for (int z = 0; z < NZ; z++)
        {
            for (int y = 0; y < NY; y++)
            {
                for (int x = 0; x < NX; x++)
                {
                    int p = index(x, y, z);

                    double distance =
                        static_cast<double>(z) - 31.5;

                    n1[p] =
                        0.5 *
                        (
                            1.0 -
                            std::tanh(
                                distance / interfaceWidth
                            )
                        );

                    n2[p] =
                        0.5 *
                        (
                            1.0 +
                            std::tanh(
                                distance / interfaceWidth
                            )
                        );
                }
            }
        }
    }

    // ------------------------------------------------------------
    // Generate 33 random bubble centers on the GB.
    // ------------------------------------------------------------

    std::vector<std::pair<double, double>>
    generateBubbleCenters()
    {
        std::vector<std::pair<double, double>> centers;

        std::uniform_real_distribution<double> randomX(
            0.0,
            static_cast<double>(NX - 1)
        );

        std::uniform_real_distribution<double> randomY(
            0.0,
            static_cast<double>(NY - 1)
        );

        int attempts = 0;

        while (
            static_cast<int>(centers.size()) < numberOfBubbles
        )
        {
            attempts++;

            if (attempts > 1000000)
            {
                std::cerr
                    << "Could not place all bubbles.\n";

                break;
            }

            double x = randomX(generator);
            double y = randomY(generator);

            bool valid = true;

            for (auto center : centers)
            {
                double dx1 = x - center.first;
                double dy1 = y - center.second;

                // Periodic distance in x.
                if (std::abs(dx1) > NX / 2.0)
                    dx1 =
                        NX -
                        std::abs(dx1);

                // Periodic distance in y.
                if (std::abs(dy1) > NY / 2.0)
                    dy1 =
                        NY -
                        std::abs(dy1);

                double distance =
                    std::sqrt(
                        dx1 * dx1 +
                        dy1 * dy1
                    );

                if (
                    distance <
                    minimumBubbleSeparation
                )
                {
                    valid = false;
                    break;
                }
            }

            if (valid)
            {
                centers.push_back(
                    {x, y}
                );
            }
        }

        return centers;
    }

    // ------------------------------------------------------------
    // Initial gas field
    // ------------------------------------------------------------

    void initializeGas()
    {
        std::vector<
            std::pair<double, double>
        > centers = generateBubbleCenters();

        std::uniform_real_distribution<double> noise(
            -noiseAmplitude,
            noiseAmplitude
        );

        double bubbleWidth = 1.5;

        for (int z = 0; z < NZ; z++)
        {
            for (int y = 0; y < NY; y++)
            {
                for (int x = 0; x < NX; x++)
                {
                    int p = index(x, y, z);

                    double value = cInitial;

                    // Create each initial bubble.
                    for (auto center : centers)
                    {
                        double dx1 =
                            static_cast<double>(x)
                            - center.first;

                        double dy1 =
                            static_cast<double>(y)
                            - center.second;

                        // Periodic x/y distance.
                        if (std::abs(dx1) > NX / 2.0)
                        {
                            if (dx1 > 0)
                                dx1 -= NX;
                            else
                                dx1 += NX;
                        }

                        if (std::abs(dy1) > NY / 2.0)
                        {
                            if (dy1 > 0)
                                dy1 -= NY;
                            else
                                dy1 += NY;
                        }

                        double dz =
                            static_cast<double>(z)
                            - 31.5;

                        double radius =
                            std::sqrt(
                                dx1 * dx1 +
                                dy1 * dy1 +
                                dz * dz
                            );

                        /*
                            Smooth spherical bubble.

                            Inside:
                                c -> 1

                            Outside:
                                c -> 0.12
                        */
                        double bubble =
                            0.5 *
                            (
                                1.0 -
                                std::tanh(
                                    (radius - initialRadius)
                                    / bubbleWidth
                                )
                            );

                        value =
                            std::max(
                                value,
                                cInitial +
                                (cBubble - cInitial)
                                * bubble
                            );
                    }

                    /*
                        Add a very small random perturbation.

                        This happens ONLY at t = 0.

                        We do not keep adding noise because gas is
                        supposed to be conserved by the CH equation.
                    */
                    value += noise(generator);

                    value =
                        std::clamp(
                            value,
                            0.0,
                            1.0
                        );

                    c[p] = value;
                }
            }
        }

        std::cout
            << "Number of initial bubbles = "
            << centers.size()
            << "\n";
    }

    // ------------------------------------------------------------
    // Initialize complete simulation
    // ------------------------------------------------------------

    void initialize()
    {
        initializeGrains();
        initializeGas();
    }

    // ------------------------------------------------------------
    // Calculate laplacians
    // ------------------------------------------------------------

    void calculateLaplacians()
    {
        #pragma omp parallel for collapse(3)
        for (int z = 0; z < NZ; z++)
        {
            for (int y = 0; y < NY; y++)
            {
                for (int x = 0; x < NX; x++)
                {
                    int p = index(x, y, z);

                    lapC[p] =
                        laplacian(c, x, y, z);

                    lapN1[p] =
                        laplacian(n1, x, y, z);

                    lapN2[p] =
                        laplacian(n2, x, y, z);
                }
            }
        }
    }

    // ------------------------------------------------------------
    // Calculate chemical potential
    // ------------------------------------------------------------

    void calculateChemicalPotential()
    {
        #pragma omp parallel for
        for (int p = 0; p < NX * NY * NZ; p++)
        {
            double cValue = c[p];

            double grainValue =
                n1[p] * n1[p]
                +
                n2[p] * n2[p];

            mu[p] =
                A0 *
                (
                    cValue * cValue * cValue
                    -
                    cValue
                    +
                    2.0 *
                    bS *
                    cValue *
                    grainValue
                )
                -
                j *
                lapC[p];
        }
    }

    // ------------------------------------------------------------
    // Mobility at one grid point
    // ------------------------------------------------------------

    double mobility(
        int p) const
    {
        double cValue = c[p];

        /*
            The factor 16*c^2*(1-c^2) is zero in the two pure
            phases and becomes large around the diffuse bubble
            interface.

            We protect against tiny numerical overshoots.
        */

        double cForMobility =
            std::clamp(
                cValue,
                0.0,
                1.0
            );

        double surfaceTerm =
            Ms *
            16.0 *
            cForMobility *
            cForMobility *
            (
                1.0 -
                cForMobility *
                cForMobility
            );

        double gbTerm =
            Mgb *
            n1[p] * n1[p] *
            n2[p] * n2[p];

        double bulkTerm =
            Mb;

        return
            surfaceTerm +
            gbTerm +
            bulkTerm;
    }

    // ------------------------------------------------------------
    // Calculate variable-mobility CH fluxes.
    //
    // Flux:
    //
    //      J = -M grad(mu)
    //
    // CH:
    //
    //      dc/dt = -div(J)
    //
    // which is the same as:
    //
    //      dc/dt = div(M grad(mu))
    //
    // ------------------------------------------------------------

    void calculateFluxes()
    {
        #pragma omp parallel for collapse(3)
        for (int z = 0; z < NZ; z++)
        {
            for (int y = 0; y < NY; y++)
            {
                for (int x = 0; x < NX; x++)
                {
                    int p = index(x, y, z);

                    // ------------------------------------------------
                    // X face
                    // ------------------------------------------------

                    int xp = periodicX(x + 1);
                    int px = index(xp, y, z);

                    double Mx =
                        0.5 *
                        (
                            mobility(p)
                            +
                            mobility(px)
                        );

                    fluxX[p] =
                        Mx *
                        (
                            mu[px] -
                            mu[p]
                        ) / dx;

                    // ------------------------------------------------
                    // Y face
                    // ------------------------------------------------

                    int yp = periodicY(y + 1);
                    int py = index(x, yp, z);

                    double My =
                        0.5 *
                        (
                            mobility(p)
                            +
                            mobility(py)
                        );

                    fluxY[p] =
                        My *
                        (
                            mu[py] -
                            mu[p]
                        ) / dx;

                    // ------------------------------------------------
                    // Z face
                    //
                    // At the top/bottom external boundary we use
                    // zero normal flux.
                    // ------------------------------------------------

                    if (z == NZ - 1)
                    {
                        fluxZ[p] = 0.0;
                    }
                    else
                    {
                        int pz =
                            index(x, y, z + 1);

                        double Mz =
                            0.5 *
                            (
                                mobility(p)
                                +
                                mobility(pz)
                            );

                        fluxZ[p] =
                            Mz *
                            (
                                mu[pz] -
                                mu[p]
                            ) / dx;
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------
    // One CH + Allen-Cahn timestep
    // ------------------------------------------------------------

    void update()
    {
        calculateLaplacians();

        calculateChemicalPotential();

        calculateFluxes();

        #pragma omp parallel for collapse(3)
        for (int z = 0; z < NZ; z++)
        {
            for (int y = 0; y < NY; y++)
            {
                for (int x = 0; x < NX; x++)
                {
                    int p = index(x, y, z);

                    // ------------------------------------------------
                    // Cahn-Hilliard
                    //
                    // dc/dt = div(M grad(mu))
                    //
                    // We calculate the divergence of the flux.
                    // ------------------------------------------------

                    int xm =
                        periodicX(x - 1);

                    int ym =
                        periodicY(y - 1);

                    int pXm =
                        index(xm, y, z);

                    int pYm =
                        index(x, ym, z);

                    double divFlux =
                        (
                            fluxX[p]
                            -
                            fluxX[pXm]
                        ) / dx
                        +
                        (
                            fluxY[p]
                            -
                            fluxY[pYm]
                        ) / dx;

                    // z direction
                    double fluxZMinus = 0.0;

                    if (z > 0)
                    {
                        int pZm =
                            index(x, y, z - 1);

                        fluxZMinus =
                            fluxZ[pZm];
                    }

                    double divFluxZ =
                        (
                            fluxZ[p]
                            -
                            fluxZMinus
                        ) / dx;

                    double dc =
                        divFlux +
                        divFluxZ;

                    cNew[p] =
                        c[p] +
                        dt * dc;

                    // ------------------------------------------------
                    // Allen-Cahn for grain 1
                    // ------------------------------------------------

                    double dFdn1 =
                        A0 *
                        (
                            n1[p] * n1[p] * n1[p]
                            -
                            n1[p]
                            +
                            bGB *
                            n1[p] *
                            n2[p] *
                            n2[p]
                            +
                            2.0 *
                            bS *
                            c[p] *
                            c[p] *
                            n1[p]
                        )
                        -
                        jGB *
                        lapN1[p];

                    n1New[p] =
                        n1[p]
                        -
                        dt *
                        LGB *
                        dFdn1;

                    // ------------------------------------------------
                    // Allen-Cahn for grain 2
                    // ------------------------------------------------

                    double dFdn2 =
                        A0 *
                        (
                            n2[p] * n2[p] * n2[p]
                            -
                            n2[p]
                            +
                            bGB *
                            n2[p] *
                            n1[p] *
                            n1[p]
                            +
                            2.0 *
                            bS *
                            c[p] *
                            c[p] *
                            n2[p]
                        )
                        -
                        jGB *
                        lapN2[p];

                    n2New[p] =
                        n2[p]
                        -
                        dt *
                        LGB *
                        dFdn2;

                    /*
                        Numerical protection.

                        The physical phase fields are between 0 and 1.
                    */

                    cNew[p] =
                        std::clamp(
                            cNew[p],
                            0.0,
                            1.0
                        );

                    n1New[p] =
                        std::clamp(
                            n1New[p],
                            0.0,
                            1.0
                        );

                    n2New[p] =
                        std::clamp(
                            n2New[p],
                            0.0,
                            1.0
                        );
                }
            }
        }

        c.swap(cNew);
        n1.swap(n1New);
        n2.swap(n2New);
    }

    // ------------------------------------------------------------
    // Calculate projected bubble mask on the GB plane.
    //
    // A pixel (x,y) belongs to the projected bubble area if
    // any z position has c >= 0.5.
    // ------------------------------------------------------------

    std::vector<unsigned char>
    projectedBubbleMask() const
    {
        std::vector<unsigned char> mask(
            NX * NY,
            0
        );

        #pragma omp parallel for collapse(2)
        for (int y = 0; y < NY; y++)
        {
            for (int x = 0; x < NX; x++)
            {
                bool bubbleFound = false;

                for (int z = 0; z < NZ; z++)
                {
                    int p =
                        index(x, y, z);

                    if (c[p] >= 0.5)
                    {
                        bubbleFound = true;
                        break;
                    }
                }

                if (bubbleFound)
                {
                    mask[y * NX + x] = 1;
                }
            }
        }

        return mask;
    }

    // ------------------------------------------------------------
    // Find connected components in the projected GB plane.
    //
    // We use 8-neighbour connectivity because two bubble regions
    // touching diagonally should be treated as connected.
    // ------------------------------------------------------------

    std::vector<int>
    projectedBubbleAreas(
        const std::vector<unsigned char>& mask) const
    {
        std::vector<unsigned char> visited(
            NX * NY,
            0
        );

        std::vector<int> areas;

        std::vector<int> queueX(
            NX * NY
        );

        std::vector<int> queueY(
            NX * NY
        );

        int dxNeighbour[8] =
        {
            -1, 0, 1,
            -1, 1,
            -1, 0, 1
        };

        int dyNeighbour[8] =
        {
            -1, -1, -1,
            0, 0,
            1, 1, 1
        };

        for (int y = 0; y < NY; y++)
        {
            for (int x = 0; x < NX; x++)
            {
                int start =
                    y * NX + x;

                if (
                    mask[start] == 0 ||
                    visited[start] != 0
                )
                {
                    continue;
                }

                int head = 0;
                int tail = 0;

                queueX[tail] = x;
                queueY[tail] = y;
                tail++;

                visited[start] = 1;

                int area = 0;

                while (head < tail)
                {
                    int cx = queueX[head];
                    int cy = queueY[head];

                    head++;

                    area++;

                    for (int q = 0; q < 8; q++)
                    {
                        int nx =
                            cx +
                            dxNeighbour[q];

                        int ny =
                            cy +
                            dyNeighbour[q];

                        // Periodic x/y.
                        nx = periodicX(nx);
                        ny = periodicY(ny);

                        int np =
                            ny * NX + nx;

                        if (
                            mask[np] != 0 &&
                            visited[np] == 0
                        )
                        {
                            visited[np] = 1;

                            queueX[tail] = nx;
                            queueY[tail] = ny;

                            tail++;
                        }
                    }
                }

                areas.push_back(area);
            }
        }

        return areas;
    }

    // ------------------------------------------------------------
    // Calculate metrics
    // ------------------------------------------------------------

    void writeMetrics(
        std::ofstream& output,
        int step)
    {
        std::vector<unsigned char> mask =
            projectedBubbleMask();

        std::vector<int> areas =
            projectedBubbleAreas(mask);

        double totalArea = 0.0;
        double largestArea = 0.0;

        for (int area : areas)
        {
            totalArea +=
                static_cast<double>(area);

            if (
                static_cast<double>(area)
                >
                largestArea
            )
            {
                largestArea =
                    static_cast<double>(area);
            }
        }

        double gbArea =
            static_cast<double>(NX * NY);

        double coverage =
            totalArea /
            gbArea;

        double largestFraction = 0.0;

        if (totalArea > 0.0)
        {
            largestFraction =
                largestArea /
                totalArea;
        }

        // Dimensionless time.
        double time =
            step * dt;

        // Convert to physical seconds.
        //
        // t_physical =
        // t_dimensionless * dx^2 / Dphysical
        //
        double timeSeconds =
            time *
            dxPhysical_m *
            dxPhysical_m /
            DPhysical;

        output
            << std::setprecision(10)
            << step << ","
            << time << ","
            << timeSeconds << ","
            << areas.size() << ","
            << totalArea << ","
            << largestArea << ","
            << largestFraction << ","
            << coverage
            << "\n";
    }

    // ------------------------------------------------------------
    // Save the 3-D gas field.
    //
    // Python will read this file and create the green 3-D
    // morphology using the c = 0.5 isosurface.
    // ------------------------------------------------------------

    void saveMorphology(int step)
    {
        std::string filename =
            "morphology_step_" +
            std::to_string(step) +
            ".raw";

        std::ofstream output(
            filename,
            std::ios::binary
        );

        if (!output)
        {
            std::cerr
                << "Could not write "
                << filename
                << "\n";

            return;
        }

        for (double value : c)
        {
            float valueFloat =
                static_cast<float>(value);

            output.write(
                reinterpret_cast<char*>(
                    &valueFloat
                ),
                sizeof(float)
            );
        }

        output.close();
    }

    // ------------------------------------------------------------
    // Run simulation
    // ------------------------------------------------------------

    void run()
    {
        initialize();

        std::ofstream output(
            "bubble_metrics.csv"
        );

        if (!output)
        {
            std::cerr
                << "Could not open bubble_metrics.csv\n";

            return;
        }

        output
            << "step,"
            << "time_dimensionless,"
            << "time_seconds,"
            << "number_of_projected_bubbles,"
            << "total_projected_area,"
            << "largest_projected_area,"
            << "largest_area_fraction,"
            << "coverage\n";

        std::cout
            << "\n";
        std::cout
            << "3-D grain-boundary bubble simulation\n";
        std::cout
            << "Domain = "
            << NX << " x "
            << NY << " x "
            << NZ << "\n";

        std::cout
            << "Initial bubbles = "
            << numberOfBubbles
            << "\n";

        std::cout
            << "Initial radius = "
            << initialRadius
            << " grid cells\n";

        std::cout
            << "Contact angle = "
            << THETA
            << " degrees\n";

        std::cout
            << "bS = "
            << bS
            << "\n";

        std::cout
            << "bGB = "
            << bGB
            << "\n";

        std::cout
            << "Mb = "
            << Mb
            << "\n";

        std::cout
            << "Mgb = "
            << Mgb
            << "\n";

        std::cout
            << "Ms = "
            << Ms
            << "\n";

        std::cout
            << "dt = "
            << dt
            << "\n";

        std::cout
            << "Total steps = "
            << totalSteps
            << "\n\n";

        writeMetrics(
            output,
            0
        );

        saveMorphology(0);

        for (
            int step = 1;
            step <= totalSteps;
            step++
        )
        {
            update();

            if (
                step % outputEvery == 0
            )
            {
                writeMetrics(
                    output,
                    step
                );

                saveMorphology(
                    step
                );

                std::cout
                    << "step = "
                    << step
                    << " / "
                    << totalSteps
                    << "\n";
            }
        }

        output.close();

        std::cout
            << "\nSimulation finished.\n";

        std::cout
            << "Metrics: bubble_metrics.csv\n";

        std::cout
            << "Morphology files: morphology_step_*.raw\n";
    }
};

int main()
{
    PhaseField3D simulation;

    simulation.run();

    return 0;
}
