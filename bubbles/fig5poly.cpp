// ============================================================================
//  FIG 5 (polycrystalline variant) : recrystallization along ALL original GBs
//  of a coarse polycrystal, widening with fission density (Fig 5 mechanism on a
//  polycrystal instead of a single central boundary).
//  Colour code: white = grain boundary, red = bubble, blue = grain interior.
//
//  A coarse polycrystal is laid down first (its boundaries are the original GBs).
//  Recrystallization nucleates on those boundaries and grows outward: every cell
//  within a distance hw(F) of an original GB becomes fine recrystallized grain,
//  so bands widen along the whole GB network as fission density rises and
//  eventually consume the coarse grains. Bubbles sit on the boundaries and
//  triple junctions of the recrystallized grains, increasing in number with F.
//
//  Built directly at each fission density (reconstruction, see fig5.cpp note).
//  Coarse structure, its GB distance transform, and the fine polycrystal are all
//  fixed; only hw(F) and the bubble count change across panels, so the bands grow
//  continuously.
//
//  Output per panel: h_<tag>.csv, sg_<tag>.csv.  Plot: python fig5poly_plot.py
//  Build : g++ -O3 -o fig5poly fig5poly.cpp     Run : ./fig5poly
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
static const unsigned SEED = 21;
static const int NCOARSE = 7;              // original coarse grains
static const int NFINE   = 500;            // recrystallized fine-grain pool
static const int OFFSET  = NFINE;          // coarse labels start here

inline int wrap(int a){ return (a%N+N)%N; }
inline int id(int i,int j){ return wrap(j)*N + wrap(i); }

int main(){
    const string outdir="fig5poly_out"; mkdir(outdir.c_str(),0777);
    mt19937 rng(SEED);
    uniform_real_distribution<double> U(0,1);

    // --- coarse polycrystal (original grains) ---
    vector<double> cx(NCOARSE),cy(NCOARSE);
    for(int g=0;g<NCOARSE;g++){ cx[g]=U(rng)*N; cy[g]=U(rng)*N; }
    vector<int> coarse(NN,0);
    for(int j=0;j<N;j++) for(int i=0;i<N;i++){
        int best=0; double bd=1e300;
        for(int g=0;g<NCOARSE;g++){
            double dx=fabs(i-cx[g]); if(dx>N/2.0) dx=N-dx;
            double dy=fabs(j-cy[g]); if(dy>N/2.0) dy=N-dy;
            double d=dx*dx+dy*dy; if(d<bd){bd=d;best=g;}
        }
        coarse[j*N+i]=best;
    }
    // --- distance to nearest original GB (chamfer transform) ---
    const double INF=1e9;
    vector<double> D(NN,INF);
    for(int j=0;j<N;j++) for(int i=0;i<N;i++){
        int k=j*N+i, L=coarse[k];
        if(coarse[id(i+1,j)]!=L||coarse[id(i-1,j)]!=L||coarse[id(i,j+1)]!=L||coarse[id(i,j-1)]!=L)
            D[k]=0.0;
    }
    for(int j=0;j<N;j++) for(int i=0;i<N;i++){ int k=j*N+i;
        if(i>0)   D[k]=min(D[k],D[k-1]+1.0);
        if(j>0)   D[k]=min(D[k],D[k-N]+1.0);
        if(i>0&&j>0) D[k]=min(D[k],D[k-N-1]+1.4142);
        if(i<N-1&&j>0) D[k]=min(D[k],D[k-N+1]+1.4142);
    }
    for(int j=N-1;j>=0;j--) for(int i=N-1;i>=0;i--){ int k=j*N+i;
        if(i<N-1) D[k]=min(D[k],D[k+1]+1.0);
        if(j<N-1) D[k]=min(D[k],D[k+N]+1.0);
        if(i<N-1&&j<N-1) D[k]=min(D[k],D[k+N+1]+1.4142);
        if(i>0&&j<N-1) D[k]=min(D[k],D[k+N-1]+1.4142);
    }
    // --- fixed fine polycrystal (revealed inside the bands) ---
    vector<double> fx(NFINE),fy(NFINE);
    for(int g=0;g<NFINE;g++){ fx[g]=U(rng)*N; fy[g]=U(rng)*N; }
    vector<int> fine(NN,0);
    for(int j=0;j<N;j++) for(int i=0;i<N;i++){
        int best=0; double bd=1e300;
        for(int g=0;g<NFINE;g++){
            double dx=fabs(i-fx[g]); if(dx>N/2.0) dx=N-dx;
            double dy=fabs(j-fy[g]); if(dy>N/2.0) dy=N-dy;
            double d=dx*dx+dy*dy; if(d<bd){bd=d;best=g;}
        }
        fine[j*N+i]=best;
    }

    struct Panel{ const char* tag; double hw; int nbub; };
    vector<Panel> panels = {
        {"F20",  3.0,   6},   // (a) 2.0e21 : thin bands on the original GBs
        {"F35", 22.0,  40},   // (b) 3.5e21 : bands established
        {"F40", 40.0,  90},   // (c) 4.0e21 : bands wide
        {"F45", 70.0, 150},   // (d) 4.5e21 : recrystallization nearly complete
    };
    { ofstream m((outdir+"/meta.txt").c_str());
      m<<"5\nF20 2.0e21 f/cm^3\nF35 3.5e21 f/cm^3\nF40 4.0e21 f/cm^3\nF45 4.5e21 f/cm^3\n"; }

    for(const auto& P : panels){
        // combined label: fine grain inside a band, else coarse grain
        vector<int> lab(NN);
        for(int k=0;k<NN;k++) lab[k] = (D[k] < P.hw)? fine[k] : coarse[k]+OFFSET;

        vector<double> gb(NN,0.0); vector<char> tj(NN,0);
        for(int j=0;j<N;j++) for(int i=0;i<N;i++){
            int k=j*N+i, L=lab[k];
            if(lab[id(i+1,j)]!=L||lab[id(i-1,j)]!=L||lab[id(i,j+1)]!=L||lab[id(i,j-1)]!=L)
                gb[k]=1.0;
            int a=lab[id(i+1,j)],b=lab[id(i-1,j)],c=lab[id(i,j+1)],d=lab[id(i,j-1)];
            int s[5]={L,a,b,c,d}; sort(s,s+5); int u=1; for(int t=1;t<5;t++) if(s[t]!=s[t-1]) u++;
            if(u>=3) tj[k]=1;
        }
        vector<double> gbf=gb;
        for(int j=0;j<N;j++) for(int i=0;i<N;i++){ int k=j*N+i; double m=gb[k];
            m=max(m,0.5*gb[id(i+1,j)]); m=max(m,0.5*gb[id(i-1,j)]);
            m=max(m,0.5*gb[id(i,j+1)]); m=max(m,0.5*gb[id(i,j-1)]); gbf[k]=m; }
        vector<double> sg(NN); for(int k=0;k<NN;k++) sg[k]=1.0-gbf[k];

        // bubbles by number, triple junctions first then GBs, spaced
        vector<double> h(NN,0.0);
        vector<int> tjc,gbc;
        for(int k=0;k<NN;k++){ if(tj[k]) tjc.push_back(k); else if(gb[k]>0.5) gbc.push_back(k); }
        shuffle(tjc.begin(),tjc.end(),rng); shuffle(gbc.begin(),gbc.end(),rng);
        vector<int> order=tjc; order.insert(order.end(),gbc.begin(),gbc.end());
        vector<pair<int,int>> ctr; int RB=2, placed=0;
        for(int k:order){ if(placed>=P.nbub) break; int i=k%N,j=k/N; bool ok=true;
            for(auto&c:ctr){int di=i-c.first,dj=j-c.second; if(di*di+dj*dj<(2*RB+2)*(2*RB+2)){ok=false;break;}}
            if(!ok) continue;
            for(int dj=-RB;dj<=RB;dj++) for(int di=-RB;di<=RB;di++)
                if(di*di+dj*dj<=RB*RB) h[id(i+di,j+dj)]=1.0;
            ctr.push_back({i,j}); placed++; }
        printf("%s : hw=%.0f  bubbles=%d\n",P.tag,P.hw,placed);

        char f1[256],f2[256];
        snprintf(f1,sizeof f1,"%s/h_%s.csv",outdir.c_str(),P.tag);
        snprintf(f2,sizeof f2,"%s/sg_%s.csv",outdir.c_str(),P.tag);
        ofstream o1(f1),o2(f2);
        for(int j=0;j<N;j++) for(int i=0;i<N;i++){
            int k=j*N+i; o1<<h[k]<<(i+1<N?' ':'\n'); o2<<sg[k]<<(i+1<N?' ':'\n');
        }
    }
    printf("done fig5poly_out/  ->  python fig5poly_plot.py\n");
    return 0;
}
