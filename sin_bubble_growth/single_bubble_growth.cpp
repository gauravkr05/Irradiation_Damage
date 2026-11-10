/*
 single_bubble_growth.cpp

 2-D phase-field demonstration of one gas bubble growing in a solid matrix.

 Model:
   F = integral [ A0*( c^4/4 - c^2/2
                     + phi^4/4 - phi^2/2
                     + bS*c^2*phi^2 )
                  + j/2*|grad(c)|^2
                  + jphi/2*|grad(phi)|^2 ] dV

 Gas chemical potential:
   mu = dF/dc
      = A0*(c^3 - c + 2*bS*c*phi^2) - j*laplacian(c)

 Cahn-Hilliard:
   dc/dt = M * laplacian(mu)

 Mobility calibration used in the paper:
   M = D / [ A0*(2*bS - 1) ]

 Grain/matrix order parameter (Allen-Cahn):
   dphi/dt = -Lphi * dF/dphi
           = -Lphi * [ A0*(phi^3 - phi + 2*bS*c^2*phi)
                        - jphi*laplacian(phi) ]

 Here:
   c   = gas phase field; c~1 inside bubble, c~0 in matrix
   phi = matrix order parameter; phi~0 in bubble, phi~1 in matrix

 IMPORTANT:
 - This is a 2-D cross-section, so the measured "projected area" is
   simply the area where c >= 0.5.
 - It is a compact validation/demo model, NOT the full 3-D bicrystal
   simulation from the paper.
 - Periodic boundary conditions are used.
 - The initial matrix concentration is c_i = 0.12, representing
   supersaturation.
 - No gas production occurs; total gas is conserved by CH.

 Compile:
   g++ -O3 -std=c++17 single_bubble_growth.cpp -o single_bubble_growth

 Run:
   ./single_bubble_growth

 Output:
   bubble_growth.csv
*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

class PhaseField2D {
public:
    // -----------------------------
    // Numerical/grid parameters
    // -----------------------------
    static constexpr int N = 128;       // 128 x 128 grid
    static constexpr double dx = 1.0;   // dimensionless grid spacing

    // We use a small explicit timestep because CH contains a
    // fourth-order spatial derivative.
    static constexpr double dt = 5.0e-4;
    static constexpr int nSteps = 100000;
    static constexpr int dumpEvery = 1000;

    // -----------------------------
    // Physical/model parameters
    // -----------------------------
    static constexpr double A0 = 4.0;
    static constexpr double bS = 1.0;

    // Gradient coefficient.
    // The exact paper uses its interface-width calibration.
    // j=2 is consistent with j/bS = 2 for bS=1.
    static constexpr double j = 2.0;

    // Grain/matrix interface relaxation rate.
    static constexpr double Lphi = 0.5;
    static constexpr double jphi = 2.0;

    // Dimensionless gas diffusivity.
    // Change this if you want a different diffusion timescale.
    static constexpr double D = 1.0;

    // From the paper:
    // M = D / [ A0 * (2*bS - 1) ]
    static constexpr double M =
        D / (A0 * (2.0 * bS - 1.0));

    // Initial gas concentration in the matrix.
    static constexpr double ci = 0.12;

    // Initial bubble concentration.
    static constexpr double cb = 1.0;

    // Initial bubble radius in grid units.
    // With dx_phys = 10 nm, R0=8 corresponds to about 80 nm.
    static constexpr double R0 = 8.0;

    // Smoothness of initial interface.
    static constexpr double interfaceWidth = 2.0;

    // Physical dimensionalization used only for the CSV.
    // Paper value: D_g = 3.4e-18 m^2/s and length scale ~10 nm.
    static constexpr double dxPhysical_m = 10.0e-9;
    static constexpr double DPhysical_m2s = 3.4e-18;

    std::vector<double> c;
    std::vector<double> phi;
    std::vector<double> cNew;
    std::vector<double> phiNew;
    std::vector<double> mu;
    std::vector<double> lapC;
    std::vector<double> lapPhi;
    std::vector<double> lapMu;

    PhaseField2D()
        : c(N * N, 0.0),
          phi(N * N, 0.0),
          cNew(N * N, 0.0),
          phiNew(N * N, 0.0),
          mu(N * N, 0.0),
          lapC(N * N, 0.0),
          lapPhi(N * N, 0.0),
          lapMu(N * N, 0.0) {}

    inline int id(int i, int j) const {
        return i * N + j;
    }

    inline int periodic(int x) const {
        if (x < 0) return x + N;
        if (x >= N) return x - N;
        return x;
    }

    inline double laplacian(const std::vector<double>& a,
                            int i, int k) const {
        const int im = periodic(i - 1);
        const int ip = periodic(i + 1);
        const int km = periodic(k - 1);
        const int kp = periodic(k + 1);

        return (a[id(im, k)] + a[id(ip, k)]
              + a[id(i, km)] + a[id(i, kp)]
              - 4.0 * a[id(i, k)]) / (dx * dx);
    }

    void initialize() {
        const double center = 0.5 * (N - 1);

        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < N; ++k) {
                const double x = (i - center) * dx;
                const double y = (k - center) * dx;
                const double r = std::sqrt(x * x + y * y);

                // Smooth bubble:
                // c -> 1 inside, c -> ci outside.
                const double s =
                    0.5 * (1.0 - std::tanh(
                        (r - R0) / interfaceWidth));

                c[id(i, k)] = ci + (cb - ci) * s;

                // phi -> 0 inside bubble and ->1 in matrix.
                phi[id(i, k)] =
                    0.5 * (1.0 + std::tanh(
                        (r - R0) / interfaceWidth));
            }
        }
    }

    void computeLaplacians() {
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < N; ++k) {
                const int p = id(i, k);
                lapC[p] = laplacian(c, i, k);
                lapPhi[p] = laplacian(phi, i, k);
            }
        }
    }

    void computeChemicalPotential() {
        for (int p = 0; p < N * N; ++p) {
            mu[p] =
                A0 * (c[p] * c[p] * c[p]
                    - c[p]
                    + 2.0 * bS * c[p] * phi[p] * phi[p])
                - j * lapC[p];
        }
    }

    void computeLapMu() {
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < N; ++k) {
                lapMu[id(i, k)] = laplacian(mu, i, k);
            }
        }
    }

    void update() {
        computeLaplacians();
        computeChemicalPotential();
        computeLapMu();

        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < N; ++k) {
                const int p = id(i, k);

                // Cahn-Hilliard update.
                cNew[p] = c[p] + dt * M * lapMu[p];

                // Allen-Cahn update for the matrix/bubble interface.
                const double dF_dphi =
                    A0 * (phi[p] * phi[p] * phi[p]
                        - phi[p]
                        + 2.0 * bS * c[p] * c[p] * phi[p])
                    - jphi * lapPhi[p];

                phiNew[p] = phi[p] - dt * Lphi * dF_dphi;

                // Keep numerical values in a physical range.
                // This is only a numerical safeguard.
                cNew[p] = std::clamp(cNew[p], 0.0, 1.2);
                phiNew[p] = std::clamp(phiNew[p], 0.0, 1.0);
            }
        }

        c.swap(cNew);
        phi.swap(phiNew);
    }

    double projectedArea() const {
        // In this 2-D cross-section, define bubble area as the region
        // with c >= 0.5.
        long long count = 0;

        for (double value : c) {
            if (value >= 0.5)
                ++count;
        }

        return static_cast<double>(count) * dx * dx;
    }

    double equivalentRadius() const {
        return std::sqrt(projectedArea() / M_PI);
    }

    double averageMatrixConcentration() const {
        double sum = 0.0;
        for (double value : c)
            sum += value;
        return sum / static_cast<double>(N * N);
    }

    void writeCSVHeader(std::ofstream& out) {
        out << "step,time_dimensionless,time_seconds,"
               "projected_area_grid2,projected_area_um2,"
               "equivalent_radius_grid,"
               "equivalent_radius_nm,"
               "average_c\n";
    }

    void writeCSVRow(std::ofstream& out, int step) {
        const double tDim = step * dt;

        // Since dx is the length scale:
        // t_physical = t_dimensionless * dxPhysical^2 / DPhysical.
        const double tSeconds =
            tDim * dxPhysical_m * dxPhysical_m
                / DPhysical_m2s;

        const double area = projectedArea();
        const double areaPhysicalUm2 =
            area * dxPhysical_m * dxPhysical_m * 1.0e12;

        const double r = equivalentRadius();
        const double rNm = r * dxPhysical_m * 1.0e9;

        out << std::fixed << std::setprecision(8)
            << step << ","
            << tDim << ","
            << tSeconds << ","
            << area << ","
            << areaPhysicalUm2 << ","
            << r << ","
            << rNm << ","
            << averageMatrixConcentration()
            << "\n";
    }

    void run() {
        initialize();

        std::ofstream out("bubble_growth.csv");
        if (!out) {
            std::cerr << "Could not open bubble_growth.csv\n";
            return;
        }

        writeCSVHeader(out);

        std::cout << "2-D single-bubble phase-field simulation\n";
        std::cout << "Grid: " << N << " x " << N << "\n";
        std::cout << "dt = " << dt << "\n";
        std::cout << "D  = " << D << "\n";
        std::cout << "M  = " << M << "\n";
        std::cout << "Initial matrix c = " << ci << "\n";
        std::cout << "Initial radius = " << R0 << " grid units\n\n";

        writeCSVRow(out, 0);

        for (int step = 1; step <= nSteps; ++step) {
            update();

            if (step % dumpEvery == 0) {
                writeCSVRow(out, step);

                std::cout
                    << "step = " << step
                    << " / " << nSteps
                    << ", area = " << projectedArea()
                    << ", R_eq = " << equivalentRadius()
                    << "\n";
            }
        }

        out.close();

        std::cout << "\nFinished.\n";
        std::cout << "Results written to bubble_growth.csv\n";
    }
};

int main() {
    PhaseField2D sim;
    sim.run();
    return 0;
}
