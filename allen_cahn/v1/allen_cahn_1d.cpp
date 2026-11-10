// ============================================================
//  1D Allen-Cahn equation  (non-conserved order parameter)
//
//      d(eta)/dt = -L * delta(F/Nv)/delta(eta)
//                = -L ( df/d(eta) - 2*kappa*d2(eta)/dx2 )
//
//  Bulk free energy :  f(eta) = A eta^2 (1 - eta)^2
//      df/d(eta) = A ( 2*eta - 6*eta^2 + 4*eta^3 )
//
//  Parameters from the sheet :  kappa = 1 , L = 1
//
//  I.C. :  background eta = 0.2 everywhere,
//          with a block of eta = 1.0 of width 50, centred in the domain.
//          Domain length = 512.
//
//  NOTE : Allen-Cahn is NOT conserved -- the integral of eta is free to
//         change. That is the whole difference from Cahn-Hilliard, so there
//         is no conservative/non-conservative choice to make here, and no
//         4th derivative: only a 3-point laplacian is needed.
// ============================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/stat.h>          
using namespace std;

int main()
{
    // ---------------- Input ----------------
    double L     = 1.0;        // kinetic coefficient
    double kappa = 1.0;        // gradient energy coefficient
    double A     = 1.0;        // free energy barrier height

    double dx = 1.0;
    double dt = 0.01;
    int    nx = 512;                               // domain 0 .. 512
    double total_time = 400.0;
    long   total_timesteps = (long)(total_time / dt + 0.5);
    int    n_saves = 5;
    int    save_times[5] = {0, 1, 5, 200, 400};

    double eta_background = 0.2;                   // starting value outside the block
    double eta_block      = 1.0;                   // value inside the block
    int    block_width    = 50;                    // the "gap" of 50
    int    block_lo = (nx - block_width) / 2;      // centred: 231
    int    block_hi = block_lo + block_width;      //          281

    const string tag = "ac";

    // ---------------- Stability check ----------------
    // gradient term behaves like diffusion with coefficient 2*kappa*L
    double r = 2.0 * kappa * L * dt / (dx * dx);
    cout << "Allen-Cahn : L = " << L << ", kappa = " << kappa << ", A = " << A
         << ", dt = " << dt << ", nx = " << nx << "\n";
    cout << "2*kappa*L*dt/dx^2 = " << r << "  (needs <= 0.5)\n";
    cout << "block : i = " << block_lo << " .. " << block_hi
         << "  (width " << block_width << ")\n";
    if (r > 0.5) { cout << "!! UNSTABLE : reduce dt\n"; return 1; }

    // ---------------- Variable Decl ----------------
    vector<double> eta    (nx, 0.0);
    vector<double> eta_new(nx, 0.0);

    // ---------------- I.C. ----------------
    for (int i = 0; i < nx; i++)
    {
        if (i >= block_lo && i < block_hi)  eta[i] = eta_block;
        else                                eta[i] = eta_background;
    }

    // ---------------- Time loop ----------------
    for (long t = 0; t <= total_timesteps; t++)
    {
        if (t > 0)
        {
            for (int i = 0; i < nx; i++)
            {
                // B.C. : periodicclass Solution {
public:
    int unit(int n){
        while(n!=0){
        if(n/10==0)return n%10;
        else n=n/10;
        }
        return 0;
    }
    vector<int> prevPermOpt1(vector<int>& arr) {
        int n=arr.size();
        vector<int>hash(10,-1);
        for(int i=0;i<n;i++){
            hash[unit(arr[i])]=i;
        }
        for(int i=0;i<n;i++){
            for(int j=9;j<=1;j--){
                int num=unit(arr[i]);
                if(num>j&&i<hash[j]){
                    swap(arr[i],arr[hash[j]]);
                    return arr;
                }
            }
        }
        return arr;
    }
};
                int im1 = (i - 1 + nx) % nx;
                int ip1 = (i + 1) % nx;

                double ei = eta[i];

                // df/d(eta)
                double d1f = A * (2.0*ei - 6.0*ei*ei + 4.0*ei*ei*ei);

                // central difference laplacian
                double lap = (eta[ip1] - 2.0*ei + eta[im1]) / (dx*dx);

                // D.E.
                eta_new[i] = ei + dt * ( -L * (d1f - 2.0*kappa*lap) );
            }
            eta.swap(eta_new);
        }

        // ---------------- Output ----------------
        double time = t * dt;
        for (int k = 0; k < n_saves; k++)
        {
            if (t == (long)(save_times[k] / dt + 0.5))
            {
                ostringstream folder;
                folder << tag << "_time_" << save_times[k];
                mkdir(folder.str().c_str(), 0777);

                ostringstream fname;
                fname << folder.str() << "/eta_time_" << save_times[k] << ".dat";
                ofstream out(fname.str().c_str());
                out << "x eta\n";
                for (int i = 0; i < nx; i++)
                    out << i*dx << " " << eta[i] << "\n";
                out.close();

                cout << "saved time = " << time << "\n";
            }
        }
    }
    return 0;
}
