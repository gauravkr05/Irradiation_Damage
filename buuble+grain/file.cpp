// ============================================================
//  3D intergranular gas-bubble growth, coalescence and percolation
//  on a grain boundary.
//  Model of  P.C. Millett, M.R. Tonks, S.B. Biner,
//  J. Nucl. Mater. 425 (2012) 130-135.
//
//  Coupled Cahn-Hilliard (conserved gas concentration cg)
//        +  Allen-Cahn   (non-conserved grain order parameters phi_i)
//
//      dc/dt      =  div( Mg grad(mu) ),   mu = df/dc - div(kappa grad c)   (9)
//      dphi_i/dt  = -L ( df/dphi_i - div(kappa grad phi_i) )                (10)
//
//  Bulk free energy (Eqs. 2-5) :
//      f = Ao [   (c^4/4 - c^2/2)                       <- bubble well  (3)
//               + sum_i (phi_i^4/4 - phi_i^2/2)         <- grain wells  (4)
//               + sum_{i<j} bGB phi_i^2 phi_j^2         <- grain-grain  (5)
//               + sum_i      bS  c^2 phi_i^2 ]          <- grain-bubble (5)
//
//  Minima :  bubble  (c = 1, phi_i = 0)
//            grain p (phi_p = 1, others 0, c = 0)   [needs bS > 0.5]
//
//  kappa is field dependent (Eq. 7) so the interface WIDTH stays fixed while
//  the interface ENERGY is varied :  kGB/bGB = kS/bS = 2.
//  Since gamma ~ sqrt(kappa*Df) ~ beta, the contact angle (Eq. 8) is
//      cos(theta) = gGB / (2 gS) = bGB / (2 bS).
//
//  Mg = Dg / ( Ao (2 bS - 1) )   (Eq. 13), so Fick's law is recovered in the
//  solid whatever the surface energy is.  Verify with  --validate .
//
//  Geometry :  two grains, ONE flat GB at x = nx/2.
//              periodic in y,z ; zero-flux in x, so the box holds exactly one
//              grain boundary and the gas reservoir is nx/2 cells on each side.
//
// ------------------------------------------------------------
//  A NOTE ON TIME  (read this before comparing with the paper)
// ------------------------------------------------------------
//  The paper quotes l = 10 nm and t0 = l^2/Dg = 625 s, and runs for 8.8 h.
//  Those three numbers do not close:
//      8.8 h / 625 s          = 50.7 reduced time units
//      diffusion length       = sqrt(50.7) = 7.1 cells = 71 nm
//      Eq. 14 bubble growth   = 0.12*sqrt(50.7) = 0.85 cells
//  i.e. in "8.8 h" a bubble grows by less than one cell and nothing coalesces.
//  (Also (10 nm)^2 / 3.4e-18 m^2/s = 29 s, not 625 s.)
//
//  Coalescence needs the bubbles to grow by roughly half their spacing:
//      rho0 = 10 /um^2 -> spacing 32 cells -> t ~ 9.0e3 reduced
//      rho0 = 20 /um^2 -> spacing 22 cells -> t ~ 3.2e3 reduced
//  and the slab must drain, which costs (nx/2)^2 reduced.
//
//  So TIME IS SET IN REDUCED UNITS here, not hours, and the slab is made
//  thinner while the supersaturation is raised to keep the GAS BUDGET per
//  unit grain-boundary area equal to the paper's:
//
//      budget = c0 * nx  =  0.12 * 100  =  12   (cells^3 of gas per cell^2 GB)
//
//  A thinner slab reaches the same bubble SIZE in far less time; only the
//  porosity (bubble volume / box volume) is rescaled by 100/nx.  Use
//  --paper for the literal paper geometry if you are willing to wait.
//
// ------------------------------------------------------------
//  OUTPUT  ->  <tag>/metrics.dat , columns :
//      t  t_hours  porosity  n_bub  Vmax/Vtot  Omega_C  Omega_V  Abar  rho_b
//                            ^^^^^^^^^^^^^^^^
//  Vmax/Vtot = V_largest_bubble / V_total_bubble  is the curve in the sketch:
//  it sits near 1/N while the bubbles are isolated, then jumps to ~1 when they
//  merge into one percolated tunnel network.  Omega_V (vented GB coverage,
//  paper Fig. 6) is the paper's version of the same transition.
//
//  build :  g++ -O3 -march=native -fopenmp gb_bubble_percolation.cpp -o gbbub
//  run   :  ./gbbub --validate                 Fick's-law check, seconds
//           ./gbbub 40 10                      theta=40 deg, rho0=10 /um^2
//           ./gbbub 60 10 100 800 40 0.30      theta ny t_end nx c0
//           ./gbbub 80 10 100 800 40 0.30 --vtk
//           ./gbbub 40 10 --paper              literal paper geometry (slow)
// ============================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <sys/stat.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

// ---------------- Grid globals ----------------
int nx, ny, nz;
long ntot;

inline long ID(int i, int j, int k) { return ((long)i * ny + j) * nz + k; }
inline int  wrap(int a, int n)      { return (a + n) % n; }                      // periodic
inline int  clampi(int a, int n)    { return a < 0 ? 0 : (a >= n ? n - 1 : a); } // zero-flux

// ============================================================
//  connected-component labelling, 6-neighbour
//  3D : periodic in y,z , free in x
// ============================================================
int label3d(const vector<double>& c, double thr, vector<int>& lab, vector<long>& sizes)
{
    fill(lab.begin(), lab.end(), -1);
    sizes.clear();
    vector<long> st;
    int nlab = 0;

    for (long s = 0; s < ntot; s++)
    {
        if (c[s] <= thr || lab[s] != -1) continue;

        st.clear();
        st.push_back(s);
        lab[s] = nlab;
        long cnt = 0;

        while (!st.empty())
        {
            long q = st.back(); st.pop_back(); cnt++;

            int k = (int)(q % nz);
            int j = (int)((q / nz) % ny);
            int i = (int)(q / ((long)nz * ny));

            long nb[6]; int m = 0;
            if (i > 0)      nb[m++] = ID(i - 1, j, k);
            if (i < nx - 1) nb[m++] = ID(i + 1, j, k);
            nb[m++] = ID(i, wrap(j - 1, ny), k);
            nb[m++] = ID(i, wrap(j + 1, ny), k);
            nb[m++] = ID(i, j, wrap(k - 1, nz));
            nb[m++] = ID(i, j, wrap(k + 1, nz));

            for (int p = 0; p < m; p++)
                if (c[nb[p]] > thr && lab[nb[p]] == -1)
                {
                    lab[nb[p]] = nlab;
                    st.push_back(nb[p]);
                }
        }
        sizes.push_back(cnt);
        nlab++;
    }
    return nlab;
}

// 2D labelling of the GB-plane projection, periodic in both directions
int label2d(const vector<char>& m, vector<int>& lab)
{
    fill(lab.begin(), lab.end(), -1);
    vector<int> st;
    int nlab = 0;
    int n2 = ny * nz;

    for (int s = 0; s < n2; s++)
    {
        if (!m[s] || lab[s] != -1) continue;
        st.clear(); st.push_back(s); lab[s] = nlab;

        while (!st.empty())
        {
            int q = st.back(); st.pop_back();
            int k = q % nz, j = q / nz;

            int nb[4];
            nb[0] = wrap(j - 1, ny) * nz + k;
            nb[1] = wrap(j + 1, ny) * nz + k;
            nb[2] = j * nz + wrap(k - 1, nz);
            nb[3] = j * nz + wrap(k + 1, nz);

            for (int p = 0; p < 4; p++)
                if (m[nb[p]] && lab[nb[p]] == -1) { lab[nb[p]] = nlab; st.push_back(nb[p]); }
        }
        nlab++;
    }
    return nlab;
}

// ============================================================
//  Fick's-law validation : 1D slab, bubble at x=0, gas supersaturation
//  in the solid.  Eq. 14 :  r - ro = (ci - ceq)/(cb - ceq) * sqrt(Dg t)
//  Run for several bS ; all must collapse on the same sqrt curve.
// ============================================================
void validate()
{
    // 1D, ONE order parameter (single grain, no GB) + gas.
    //   bubble  : c = 1, phi = 0        solid : c = ci, phi = 1
    // Both fields are evolved with the same equations as the 3D code, so this
    // tests the real formulation, not a stand-in for it.
    const double Ao = 4.0, Dg = 1.0, dx = 1.0, Lphi = 1.0;
    const int    n  = 3000;
    const double ci = 0.1;                  // paper Eq. 14 uses ci = 0.1
    const double betaS_list[3] = {1.0, 1.3, 1.7};

    mkdir("validate", 0777);
    ofstream o("validate/fick.dat");
    o << "# 1D bubble growth, Eq. 14 check\n";
    o << "# betaS  t  r_minus_r0  analytic\n";

    printf("\n--- Fick's law validation (Eq. 13 / Eq. 14) ---\n");
    printf("Eq. 14 :  r - ro = (ci-ceq)/(cb-ceq) * sqrt(Dg t) = %.2f*sqrt(t)\n", ci);
    printf("%8s %10s %12s %12s %10s\n", "betaS", "Mg", "r-r0 (sim)", "analytic", "ratio");

    for (int b = 0; b < 3; b++)
    {
        double betaS = betaS_list[b];
        double kS    = 2.0 * betaS;
        double Mg    = Dg / (Ao * (2.0 * betaS - 1.0));

        double lapmax = 4.0 / (dx * dx);
        double dtc = 2.0 / (Mg * (kS * lapmax * lapmax + 3.0 * Ao * lapmax));
        double dta = 2.0 / (Lphi * (kS * lapmax + 4.0 * Ao));
        double dt  = 0.35 * min(dtc, dta);
        if (dt > 0.005) dt = 0.005;

        vector<double> c(n), cn(n), mu(n), ph(n), phn(n);
        double r0 = 20.0;
        for (int i = 0; i < n; i++)
        {
            double s = 0.5 * (1.0 - tanh((i + 0.5 - r0) / 1.5));   // 1 inside bubble
            c [i] = max(ci, s);
            ph[i] = 1.0 - s;
        }

        double tend = 400.0;
        long   ns   = (long)(tend / dt);
        long   outn = max(1L, ns / 40);

        for (long t = 1; t <= ns; t++)
        {
            for (int i = 0; i < n; i++)
            {
                int im = clampi(i - 1, n), ip = clampi(i + 1, n);
                double cc = c[i], pp = ph[i];
                double lapc = (c[ip] - 2.0 * cc + c[im]) / (dx * dx);
                mu[i] = Ao * (cc * cc * cc - cc + 2.0 * betaS * cc * pp * pp) - kS * lapc;
            }
            for (int i = 0; i < n; i++)
            {
                int im = clampi(i - 1, n), ip = clampi(i + 1, n);
                double cc = c[i], pp = ph[i];
                cn[i] = cc + dt * Mg * (mu[ip] - 2.0 * mu[i] + mu[im]) / (dx * dx);

                double lapp = (ph[ip] - 2.0 * pp + ph[im]) / (dx * dx);
                double dfdp = Ao * (pp * pp * pp - pp + 2.0 * betaS * pp * cc * cc);
                phn[i] = pp - dt * Lphi * (dfdp - kS * lapp);
            }
            c.swap(cn); ph.swap(phn);

            if (t % outn == 0)
            {
                double r = 0.0;
                for (int i = 1; i < n; i++)
                    if (c[i - 1] > 0.5 && c[i] <= 0.5)
                    { r = (i - 1) + (c[i - 1] - 0.5) / (c[i - 1] - c[i]); break; }
                double tt = t * dt;
                o << betaS << " " << tt << " " << (r - r0) << " " << ci * sqrt(tt) << "\n";
            }
        }

        double r = 0.0;
        for (int i = 1; i < n; i++)
            if (c[i - 1] > 0.5 && c[i] <= 0.5)
            { r = (i - 1) + (c[i - 1] - 0.5) / (c[i - 1] - c[i]); break; }

        double ana = ci * sqrt(tend);
        printf("%8.2f %10.4f %12.4f %12.4f %10.3f\n",
               betaS, Mg, r - r0, ana, (r - r0) / ana);
    }
    o.close();
    printf("wrote validate/fick.dat\n");
    printf("(the three bS must agree -> Eq. 13 is doing its job; the ratio to\n"
           " the analytic value is common to all three and reflects the diffuse\n"
           " interface, not a bS-dependent error.)\n\n");
}

// ============================================================
int main(int argc, char** argv)
{
    // ---------------- Input ----------------
    double theta_deg = 60.0;    // contact angle to target
    double rho0      = 10.0;    // initial bubble density, bubbles / um^2
    int    ny_in     = 100;     // GB plane is ny x nz cells
    double t_end     = 800.0;   // total time, REDUCED units
    int    nx_in     = 40;      // slab thickness in cells
    double c_sol     = 0.30;    // initial gas supersaturation in the solid
    bool   write_vtk = false;
    bool   paper     = false;

    {
        vector<string> pos;
        for (int a = 1; a < argc; a++)
        {
            string s = argv[a];
            if      (s == "--vtk")      write_vtk = true;
            else if (s == "--paper")    paper     = true;
            else if (s == "--validate") { validate(); return 0; }
            else if (!s.empty() && s[0] != '-') pos.push_back(s);
        }
        if (pos.size() > 0) theta_deg = atof(pos[0].c_str());
        if (pos.size() > 1) rho0      = atof(pos[1].c_str());
        if (pos.size() > 2) ny_in     = atoi(pos[2].c_str());
        if (pos.size() > 3) t_end     = atof(pos[3].c_str());
        if (pos.size() > 4) nx_in     = atoi(pos[4].c_str());
        if (pos.size() > 5) c_sol     = atof(pos[5].c_str());
    }

    if (paper) { nx_in = 100; ny_in = 300; c_sol = 0.12; t_end = 12000.0; }

    nx = nx_in;  ny = ny_in;  nz = ny_in;
    ntot = (long)nx * ny * nz;

    const int    np   = 2;        // number of grain orientations
    const double Ao   = 4.0;      // free energy scale        (paper: 4.0)
    const double Dg   = 1.0;      // gas diffusivity, reduced units
    const double Lphi = 1.0;      // order parameter mobility
    const double dx   = 1.0;      // 1 cell = 10 nm

    const double l0_nm   = 10.0;                // length scale
    const double t0_s    = 625.0;               // paper's quoted time scale
    const double cell_um = l0_nm * 1e-3 * dx;   // 0.01 um per cell

    // ---- beta's for the requested contact angle  (paper, Section 3) ----
    double betaGB, betaS;
    if      (fabs(theta_deg - 40.0) < 1e-9) { betaGB = 1.6; betaS = 1.0; }
    else if (fabs(theta_deg - 50.0) < 1e-9) { betaGB = 1.3; betaS = 1.0; }
    else if (fabs(theta_deg - 60.0) < 1e-9) { betaGB = 1.0; betaS = 1.0; }
    else if (fabs(theta_deg - 80.0) < 1e-9) { betaGB = 0.7; betaS = 1.7; }
    else {                                    // generic : cos t = bGB/(2 bS)
        betaS  = 1.0;
        betaGB = 2.0 * cos(theta_deg * M_PI / 180.0);
        if (betaGB < 0.6) betaGB = 0.6;
        if (betaGB > 1.8) betaGB = 1.8;
    }
    const double kGB  = 2.0 * betaGB;                      // kappa/beta = 2
    const double kS   = 2.0 * betaS;
    const double kDef = 0.5 * (kGB + kS);
    const double Mg   = Dg / (Ao * (2.0 * betaS - 1.0));   // Eq. 13

    // ---------------- Time step (explicit stability) ----------------
    double lapmax = 4.0 * 3.0 / (dx * dx);    // |max eigenvalue| of 3D laplacian
    double kmax   = max(kGB, kS);
    double dt_ch  = 2.0 / (Mg * (kmax * lapmax * lapmax + 3.0 * Ao * lapmax));
    double dt_ac  = 2.0 / (Lphi * (kmax * lapmax + 4.0 * Ao));
    double dt     = 0.35 * min(dt_ch, dt_ac);
    if (dt > 0.005) dt = 0.005;

    long nsteps     = (long)(t_end / dt + 0.5);
    long meas_every = max(1L, nsteps / 150);              // ~150 points

    // ---------------- Initial bubbles ----------------
    double R0      = 44.0 / l0_nm;            // 44 nm radius  -> 4.4 cells
    double min_sep = 160.0 / l0_nm;           // 160 nm centre-to-centre
    double gb_area = (ny * cell_um) * (nz * cell_um);       // um^2
    int    nbub0   = (int)(rho0 * gb_area + 0.5);
    if (nbub0 < 1) nbub0 = 1;

    // same layout for every theta so the curves are directly comparable
    mt19937 rng(12345 + (int)(rho0 * 10));
    uniform_real_distribution<double> U(0.0, 1.0);

    vector<double> by, bz;
    {
        double sep = min_sep;
        int guard = 0;
        while ((int)by.size() < nbub0)
        {
            double y = U(rng) * ny, z = U(rng) * nz;
            bool ok = true;
            for (size_t b = 0; b < by.size(); b++)
            {
                double dy = fabs(y - by[b]); if (dy > ny / 2.0) dy = ny - dy;
                double dz = fabs(z - bz[b]); if (dz > nz / 2.0) dz = nz - dz;
                if (dy * dy + dz * dz < sep * sep) { ok = false; break; }
            }
            if (ok) { by.push_back(y); bz.push_back(z); guard = 0; }
            else if (++guard > 20000) { sep *= 0.9; guard = 0; }   // relax if crowded
        }
    }

    // ---------------- Output folder ----------------
    ostringstream tagss;
    tagss << "run_theta" << (int)theta_deg << "_rho" << (int)rho0;
    const string tag = tagss.str();
    mkdir(tag.c_str(), 0777);

    double budget = c_sol * nx;   // cells^3 of gas per cell^2 of GB

    printf("-------------------------------------------------\n");
    printf("theta      = %.0f deg   (cos = %.3f = bGB/2bS)\n", theta_deg, betaGB / (2.0 * betaS));
    printf("bGB, bS    = %.2f, %.2f\n", betaGB, betaS);
    printf("kGB, kS    = %.2f, %.2f     Mg = %.4f   Ao = %.1f\n", kGB, kS, Mg, Ao);
    printf("grid       = %d x %d x %d   (%.2f x %.2f x %.2f um)\n",
           nx, ny, nz, nx * cell_um, ny * cell_um, nz * cell_um);
    printf("GB area    = %.2f um^2,  bubbles = %d  (rho0 = %.1f /um^2)\n",
           gb_area, (int)by.size(), by.size() / gb_area);
    printf("c0         = %.3f   gas budget c0*nx = %.1f  (paper: 12.0)\n", c_sol, budget);
    printf("dt         = %.5f   steps = %ld   t_end = %.0f reduced (%.1f h at t0=%.0f s)\n",
           dt, nsteps, t_end, t_end * t0_s / 3600.0, t0_s);
    printf("slab drain time (nx/2)^2 = %.0f reduced\n", (nx / 2.0) * (nx / 2.0));
    printf("output     -> %s/\n", tag.c_str());
    printf("-------------------------------------------------\n");
    fflush(stdout);

    // ---------------- Variable Decl ----------------
    vector<double> c(ntot, 0.0), cnew(ntot, 0.0), mu(ntot, 0.0), kap(ntot, 0.0);
    vector< vector<double> > phi(np, vector<double>(ntot, 0.0));
    vector< vector<double> > phinew(np, vector<double>(ntot, 0.0));
    vector<int>  lab(ntot, -1);
    vector<long> sizes;
    vector<char> proj(ny * nz, 0);
    vector<int>  lab2(ny * nz, -1);

    // ---------------- I.C. ----------------
    // flat GB at x = nx/2, diffuse over ~2 cells ; bubbles sitting on it
    double xmid = 0.5 * nx, wgb = 2.0;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
            {
                long s = ID(i, j, k);
                double p0 = 0.5 * (1.0 - tanh((i + 0.5 - xmid) / wgb));
                double p1 = 1.0 - p0;

                double bub = 0.0;
                for (size_t b = 0; b < by.size(); b++)
                {
                    double dy = fabs(j + 0.5 - by[b]); if (dy > ny / 2.0) dy = ny - dy;
                    double dz = fabs(k + 0.5 - bz[b]); if (dz > nz / 2.0) dz = nz - dz;
                    double dxx = (i + 0.5) - xmid;
                    double d = sqrt(dxx * dxx + dy * dy + dz * dz);
                    double v = 0.5 * (1.0 - tanh((d - R0) / 1.5));
                    if (v > bub) bub = v;
                }

                c[s]      = max(c_sol, bub);
                phi[0][s] = p0 * (1.0 - bub);
                phi[1][s] = p1 * (1.0 - bub);
            }

    // ---------------- Metrics file ----------------
    ofstream met((tag + "/metrics.dat").c_str());
    met << "# theta = " << theta_deg << "  rho0 = " << by.size() / gb_area
        << " /um^2  grid " << nx << "x" << ny << "x" << nz
        << "  c0 = " << c_sol << "  budget = " << budget << "\n";
    met << "# t_reduced t_hours porosity n_bub Vmax/Vtot Omega_C Omega_V "
           "avg_area_um2 rho_b_um2\n";

    // projected bubble maps at these fractions of t_end
    double save_f[6] = {0.0, 0.1, 0.25, 0.5, 0.75, 1.0};
    int nsave = 6, isave = 0;

    // ============================================================
    //  Time loop
    // ============================================================
    for (long t = 0; t <= nsteps; t++)
    {
        // ---- (1) kappa field  (Eq. 7) ----
        #pragma omp parallel for schedule(static)
        for (long s = 0; s < ntot; s++)
        {
            double wGB = 0.0, wS = 0.0;
            for (int a = 0; a < np; a++)
            {
                double pa2 = phi[a][s] * phi[a][s];
                wS += pa2 * c[s] * c[s];
                for (int b = 0; b < np; b++)
                    if (b != a) wGB += pa2 * phi[b][s] * phi[b][s];
            }
            double eps = 1e-6;
            kap[s] = (kGB * wGB + kS * wS + eps * kDef) / (wGB + wS + eps);
        }

        // ---- (2) chemical potential  (Eq. 11, variable kappa) ----
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < nx; i++)
        {
            int im = clampi(i - 1, nx), ip = clampi(i + 1, nx);
            for (int j = 0; j < ny; j++)
            {
                int jm = wrap(j - 1, ny), jp = wrap(j + 1, ny);
                for (int k = 0; k < nz; k++)
                {
                    int km = wrap(k - 1, nz), kp = wrap(k + 1, nz);
                    long s = ID(i, j, k);
                    double cc = c[s], kc = kap[s];

                    long n6[6] = { ID(im,j,k), ID(ip,j,k), ID(i,jm,k),
                                   ID(i,jp,k), ID(i,j,km), ID(i,j,kp) };
                    double div = 0.0;
                    for (int p = 0; p < 6; p++)
                        div += 0.5 * (kc + kap[n6[p]]) * (c[n6[p]] - cc);
                    div /= (dx * dx);

                    double sp2 = 0.0;
                    for (int a = 0; a < np; a++) sp2 += phi[a][s] * phi[a][s];

                    mu[s] = Ao * (cc * cc * cc - cc + 2.0 * betaS * cc * sp2) - div;
                }
            }
        }

        // ---- (3) update c (Cahn-Hilliard) and phi (Allen-Cahn) ----
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < nx; i++)
        {
            int im = clampi(i - 1, nx), ip = clampi(i + 1, nx);
            for (int j = 0; j < ny; j++)
            {
                int jm = wrap(j - 1, ny), jp = wrap(j + 1, ny);
                for (int k = 0; k < nz; k++)
                {
                    int km = wrap(k - 1, nz), kp = wrap(k + 1, nz);
                    long s = ID(i, j, k);
                    long n6[6] = { ID(im,j,k), ID(ip,j,k), ID(i,jm,k),
                                   ID(i,jp,k), ID(i,j,km), ID(i,j,kp) };

                    // Cahn-Hilliard : dc/dt = Mg lap(mu)
                    double lapmu = -6.0 * mu[s];
                    for (int p = 0; p < 6; p++) lapmu += mu[n6[p]];
                    lapmu /= (dx * dx);
                    cnew[s] = c[s] + dt * Mg * lapmu;

                    // Allen-Cahn for every order parameter
                    double kc = kap[s], cc = c[s];
                    for (int a = 0; a < np; a++)
                    {
                        double pa = phi[a][s];
                        double sumj = 0.0;
                        for (int b = 0; b < np; b++)
                            if (b != a) sumj += phi[b][s] * phi[b][s];

                        double dfdp = Ao * (pa * pa * pa - pa
                                            + 2.0 * betaGB * pa * sumj
                                            + 2.0 * betaS  * pa * cc * cc);

                        double div = 0.0;
                        for (int p = 0; p < 6; p++)
                            div += 0.5 * (kc + kap[n6[p]]) * (phi[a][n6[p]] - pa);
                        div /= (dx * dx);

                        phinew[a][s] = pa - dt * Lphi * (dfdp - div);
                    }
                }
            }
        }
        c.swap(cnew);
        for (int a = 0; a < np; a++) phi[a].swap(phinew[a]);

        // ============================================================
        //  Measurements
        // ============================================================
        if (t % meas_every == 0 || t == nsteps)
        {
            double tr = t * dt;
            double th = tr * t0_s / 3600.0;

            // --- 3D bubbles ---
            int nb = label3d(c, 0.5, lab, sizes);
            long vtot = 0, vmax = 0;
            for (size_t b = 0; b < sizes.size(); b++)
            { vtot += sizes[b]; if (sizes[b] > vmax) vmax = sizes[b]; }
            double porosity = (double)vtot / (double)ntot;
            double frac_max = vtot > 0 ? (double)vmax / (double)vtot : 0.0;

            // --- projection on the GB plane ---
            fill(proj.begin(), proj.end(), 0);
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++)
                        if (c[ID(i, j, k)] > 0.5) proj[j * nz + k] = 1;

            long aproj = 0;
            for (int s = 0; s < ny * nz; s++) aproj += proj[s];
            double OmegaC = (double)aproj / (double)(ny * nz);
            double avgA   = nb > 0 ? aproj * cell_um * cell_um / nb : 0.0;

            // --- vented coverage : clusters inside the 0.9-side-length circle
            //     that touch its circumference  (paper, Section 3) ---
            int n2 = label2d(proj, lab2);
            double yc = 0.5 * ny, zc = 0.5 * nz;
            double rad = 0.45 * min(ny, nz);
            vector<long> ain(n2, 0);
            vector<char> touch(n2, 0);
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                {
                    int s = j * nz + k;
                    if (lab2[s] < 0) continue;
                    double rr = sqrt((j + 0.5 - yc) * (j + 0.5 - yc)
                                   + (k + 0.5 - zc) * (k + 0.5 - zc));
                    if (rr <= rad)
                    {
                        ain[lab2[s]]++;
                        if (rr > rad - 1.5) touch[lab2[s]] = 1;
                    }
                }
            long atot_in = 0, atouch = 0;
            for (int b = 0; b < n2; b++) { atot_in += ain[b]; if (touch[b]) atouch += ain[b]; }
            double OmegaV = atot_in > 0 ? (double)atouch / (double)atot_in : 0.0;

            met << tr << " " << th << " " << porosity << " " << nb << " "
                << frac_max << " " << OmegaC << " " << OmegaV << " "
                << avgA << " " << nb / gb_area << "\n";
            met.flush();

            printf("t = %8.1f  poro %.4f  nbub %4d  Vmax/Vtot %.3f  OmC %.3f  OmV %.3f\n",
                   tr, porosity, nb, frac_max, OmegaC, OmegaV);
            fflush(stdout);

            // --- projected bubble map for the figures ---
            if (isave < nsave && tr >= save_f[isave] * t_end - 1e-9)
            {
                ostringstream fn;
                fn << tag << "/proj_" << isave << ".dat";
                ofstream o(fn.str().c_str());
                o << "# t = " << tr << "\n";
                for (int j = 0; j < ny; j++)
                {
                    for (int k = 0; k < nz; k++)
                    {
                        double m = 0.0;
                        for (int i = 0; i < nx; i++) m = max(m, c[ID(i, j, k)]);
                        o << m << (k == nz - 1 ? '\n' : ' ');
                    }
                }
                o.close();

                if (write_vtk)
                {
                    ostringstream vn;
                    vn << tag << "/cg_" << isave << ".vtk";
                    ofstream v(vn.str().c_str());
                    v << "# vtk DataFile Version 3.0\ncg\nASCII\n"
                      << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << nz << " " << ny << " " << nx
                      << "\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA " << ntot
                      << "\nSCALARS cg double 1\nLOOKUP_TABLE default\n";
                    for (long s = 0; s < ntot; s++) v << c[s] << "\n";
                    v.close();
                }
                isave++;
            }
        }
    }

    met.close();
    printf("done -> %s/metrics.dat\n", tag.c_str());
    return 0;
}