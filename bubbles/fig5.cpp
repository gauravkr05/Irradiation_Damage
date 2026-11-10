// ============================================================================
//  FIG 5 : U-Mo bubble + recrystallization vs fission density
//  Morphology reconstruction of Jiang et al., Nucl. Eng. Tech. 54 (2022), Fig 5.
//  Colour code (as in the paper): white = grain boundary, red = bubble,
//  blue = grain interior.
//
//  (a) 2.0e21 : BEFORE recrystallization. Two coarse grains meet at a single
//      vertical GB down the middle; a few bubbles sit on that GB (GB is the
//      preferential, lower-energy site for bubbles to form and grow).
//  (b)-(d)    : recrystallization nucleates on that GB and grows OUTWARD as a
//      vertical band of fine grains that widens with fission density. The
//      abundant new GBs act as sinks / fast channels, so bubbles collect on the
//      boundaries between recrystallized grains. The coarse matrix on either
//      side stays blue.
//
//  Built directly at each fission density (see fig7.cpp header for why the
//  fully-coupled model is reconstructed rather than time-integrated). One fixed
//  fine polycrystal is revealed by a central band of half-width hw(F).
//
//  Output per panel: h_<tag>.csv, sg_<tag>.csv.  Plot: python fig5_plot.py
//  Build : g++ -O3 -o fig5 fig5.cpp     Run : ./fig5
// ============================================================================
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <cstdio>
#include <sys/stat.h>
using namespace std;

static const int N = 256;
static const int NN = N*N;
static const unsigned SEED = 2024;
static const int NFINE = 320;              // fine recrystallized grains (fewer = larger, clearer)
static const int GL = NFINE;               // left coarse grain label
static const int GR = NFINE+1;             // right coarse grain label

inline int wrapj(int j){ return (j%N+N)%N; }
inline int idw(int i,int j){ return wrapj(j)*N + i; }   // periodic in y only

int main(){
    const string outdir="fig5_out"; mkdir(outdir.c_str(),0777);
    mt19937 rng(SEED);
    uniform_real_distribution<double> U(0,1);

    vector<double> sx(NFINE),sy(NFINE);
    for(int g=0;g<NFINE;g++){ sx[g]=U(rng)*N; sy[g]=U(rng)*N; }

    struct Panel{ const char* tag; int hw; double cov; };
    vector<Panel> panels = {
        {"F20",   0, 0.004},  // (a) 2.0e21 : before recrystallization, single GB
        {"F35",  34, 0.030},  // (b) 3.5e21 : recrystallized band established
        {"F40",  52, 0.050},  // (c) 4.0e21 : band wider
        {"F45",  72, 0.070},  // (d) 4.5e21 : band widest
    };
    { ofstream m((outdir+"/meta.txt").c_str());
      m<<"5\nF20 2.0e21 f/cm^3\nF35 3.5e21 f/cm^3\nF40 4.0e21 f/cm^3\nF45 4.5e21 f/cm^3\n"; }

    const int cx = N/2;
    for(const auto& P : panels){
        // label field: fine grains inside |x-cx|<hw, two coarse grains outside
        vector<int> lab(NN,0);
        for(int j=0;j<N;j++) for(int i=0;i<N;i++){
            if(abs(i-cx) < P.hw){
                int best=0; double bd=1e300;
                for(int g=0;g<NFINE;g++){
                    double dx=fabs(i-sx[g]); if(dx>N/2.0) dx=N-dx;
                    double dy=fabs(j-sy[g]); if(dy>N/2.0) dy=N-dy;
                    double d=dx*dx+dy*dy; if(d<bd){bd=d;best=g;}
                }
                lab[j*N+i]=best;
            } else {
                lab[j*N+i] = (i<cx? GL : GR);
            }
        }
        // grain boundaries (non-periodic in x so no seam line, periodic in y)
        vector<double> gb(NN,0.0);
        for(int j=0;j<N;j++) for(int i=0;i<N;i++){
            int k=j*N+i, L=lab[k]; bool b=false;
            if(i+1<N && lab[k+1]!=L) b=true;
            if(i-1>=0 && lab[k-1]!=L) b=true;
            if(lab[idw(i,j+1)]!=L) b=true;
            if(lab[idw(i,j-1)]!=L) b=true;
            if(b) gb[k]=1.0;
        }
        // keep boundary lines thin: light 4-neighbour thickening only
        vector<double> gbf=gb;
        for(int j=0;j<N;j++) for(int i=0;i<N;i++){
            int k=j*N+i; double m=gb[k];
            if(i+1<N) m=max(m,0.5*gb[k+1]);
            if(i-1>=0) m=max(m,0.5*gb[k-1]);
            m=max(m,0.5*gb[idw(i,j+1)]);
            m=max(m,0.5*gb[idw(i,j-1)]);
            gbf[k]=m;
        }
        vector<double> sg(NN);
        for(int k=0;k<NN;k++) sg[k]=1.0-gbf[k];

        // bubbles on GB cells up to target coverage
        vector<double> h(NN,0.0);
        vector<int> gbcells;
        for(int k=0;k<NN;k++) if(gb[k]>0.5) gbcells.push_back(k);
        shuffle(gbcells.begin(),gbcells.end(),rng);
        long target=(long)(P.cov*NN), placed=0; size_t idxc=0; int RB=2;
        while(placed<target && idxc<gbcells.size()){
            int k=gbcells[idxc++]; int i=k%N, j=k/N;
            for(int dj=-RB;dj<=RB;dj++) for(int di=-RB;di<=RB;di++){
                if(di*di+dj*dj>RB*RB) continue;
                int ii=i+di; if(ii<0||ii>=N) continue;
                h[idw(ii,j+dj)]=1.0;
            }
            placed += (long)(M_PI*RB*RB);
        }
        char f1[256],f2[256];
        snprintf(f1,sizeof f1,"%s/h_%s.csv",outdir.c_str(),P.tag);
        snprintf(f2,sizeof f2,"%s/sg_%s.csv",outdir.c_str(),P.tag);
        ofstream o1(f1),o2(f2);
        for(int j=0;j<N;j++) for(int i=0;i<N;i++){
            int k=j*N+i; o1<<h[k]<<(i+1<N?' ':'\n'); o2<<sg[k]<<(i+1<N?' ':'\n');
        }
    }
    printf("done fig5_out/  ->  python fig5_plot.py\n");
    return 0;
}