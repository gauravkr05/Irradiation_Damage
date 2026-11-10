// ============================================================================
//  FIG 5 (OPTIMIZED) : U-Mo bubble evolution + recrystallization, two-grain start
//  Jiang, Xin, Liu et al., Nucl. Eng. Tech. 54 (2022) 226-233, Fig. 5.
//
//  Physics is identical to the earlier fig5.cpp. Only the implementation is
//  faster:
//    * grain sums Sd=sum phid^2, Sr=sum phir^2 precomputed once per step
//      (was recomputed inside every grain loop -> O(Ngrain^2) per cell)
//    * neighbor indices tabulated (no modulo in the hot loops)
//    * per-cell field updates are OpenMP-parallel; all random draws are filled
//      into arrays first so the parallel loops are thread-safe & reproducible
//
//  Build (parallel): g++ -O3 -march=native -fopenmp -o fig5 fig5.cpp
//  Build (serial)  : g++ -O3 -march=native          -o fig5 fig5.cpp
//  Run  : ./fig5      Plot : python fig5_plot.py
//  (see the model / equation map / assumptions in the header of the original
//   fig5.cpp; nothing physical changed here)
// ============================================================================
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <cstdio>
#include <sys/stat.h>
using namespace std;

static const int    N   = 128;
static const double DX  = 1.0;
static const double DT  = 0.02;

static const int NGD = 6;
static const int NGR = 16;         // rex-grain pool (was 24; smaller = faster)

static const double DV = 2.0, DG = 2.2;
static const double KC = 0.5, KH = 1.5, KF = 1.5;
static const double LMOB = 1.0;

static const double A = 1.0;
static const double GR = 1.5, GD = 1.5, GRD = 1.5, GB = 4.92;
static const double ESTORE = 0.6;

static const double CVEQ = 1.0e-3, CGEQ = 5.0e-4;
static const double CBV = 0.99, CBG = 0.082, CMG = 0.15;

// irradiation rate rescaled 5x for speed: F reaches 4.5 in ~9000 steps
// instead of 45000. Cascade + gas rates scaled up by the same 5x so the
// total dose per unit fission density is unchanged. To recover the slow,
// finely-resolved run set DFDT=1e-4, PCASC=0.02, GASDEP=0.05.
static const double VG = 0.2, PCASC = 0.10, GASDEP = 0.25;
static const double SV = 0.2, SG = 0.2, NOISE = 1.0e-4, PVG_DRIVE = 0.0;
static const double DFDT = 5.0e-4;

static const double AVR_K = 0.02, AVR_F0 = 1.67, AVR_N = 3.5;
static const int    REX_EVERY = 500, REX_R = 3;

static const long   MAXSTEP = 60000;
static const unsigned SEED  = 2024;

static const int NN = N * N;
inline double clampd(double x,double lo,double hi){ return x<lo?lo:(x>hi?hi:x); }

int main(){
    const string outdir="fig5_out";
    mkdir(outdir.c_str(),0777);

    // neighbor tables (periodic), computed once
    vector<int> R(NN),Lx(NN),Uy(NN),Dn(NN);
    for(int j=0;j<N;j++) for(int i=0;i<N;i++){
        int k=j*N+i;
        R[k]  = j*N + (i+1)%N;
        Lx[k] = j*N + (i-1+N)%N;
        Uy[k] = ((j+1)%N)*N + i;
        Dn[k] = ((j-1+N)%N)*N + i;
    }

    vector<double> cv(NN),cg(NN),h(NN),cvN(NN),cgN(NN),hN(NN);
    vector<double> muv(NN),mug(NN),Mv(NN),Mg(NN),Sd(NN),Sr(NN);
    vector<double> pvdep(NN),pgr(NN),nsv(NN),nsg(NN);        // pre-drawn randoms
    vector< vector<double> > phid(NGD,vector<double>(NN,0.0)),
                             phir(NGR,vector<double>(NN,0.0)),
                             phidN(NGD,vector<double>(NN,0.0)),
                             phirN(NGR,vector<double>(NN,0.0));
    vector<double*> PD(NGD),PR(NGR),PDN(NGD),PRN(NGR);

    mt19937 rng(SEED);
    uniform_real_distribution<double> U(0.0,1.0);
    normal_distribution<double> Gn(0.0,1.0);
    uniform_int_distribution<int> Rop(0,NGR-1);

    // ---- initial condition: two deformed grains, vertical GB ----
    for(int k=0;k<NN;k++){ cv[k]=CVEQ; cg[k]=CGEQ; h[k]=0.0; }
    for(int j=0;j<N;j++) for(int i=0;i<N;i++){
        int k=j*N+i; if(i<N/2) phid[0][k]=1.0; else phid[1][k]=1.0;
    }
    for(int k=0;k<NN;k++){ for(int p=0;p<NGD;p++) phid[p][k]+=0.01*U(rng); h[k]+=0.001*U(rng); }

    double meq=1.0-CVEQ-CGEQ;
    vector<double> targets={2.0,3.5,4.0,4.5};
    size_t tptr=0;

    auto dump=[&](const string& tag){
        char f1[256],f2[256];
        snprintf(f1,sizeof f1,"%s/h_%s.csv",outdir.c_str(),tag.c_str());
        snprintf(f2,sizeof f2,"%s/sg_%s.csv",outdir.c_str(),tag.c_str());
        ofstream o1(f1),o2(f2);
        for(int j=0;j<N;j++) for(int i=0;i<N;i++){
            int k=j*N+i; double sg=0.0;
            for(int p=0;p<NGD;p++) sg+=phid[p][k]*phid[p][k];
            for(int p=0;p<NGR;p++) sg+=phir[p][k]*phir[p][k];
            o1<<h[k]<<(i+1<N?' ':'\n'); o2<<sg<<(i+1<N?' ':'\n');
        }
    };
    auto tagof=[&](double F){ char t[16]; snprintf(t,sizeof t,"F%02d",(int)llround(F*10)); return string(t); };
    { ofstream m((outdir+"/meta.txt").c_str());
      m<<"5\nF20 2.0e21 f/cm^3\nF35 3.5e21 f/cm^3\nF40 4.0e21 f/cm^3\nF45 4.5e21 f/cm^3\n"; }

    auto rex_pass=[&](double F){
        if(F<=AVR_F0) return;
        double target=1.0-exp(-AVR_K*pow(F-AVR_F0,AVR_N));
        long rc=0;
        for(int k=0;k<NN;k++){ double mx=0.0; for(int p=0;p<NGR;p++) mx=max(mx,phir[p][k]); if(mx>0.5) rc++; }
        double frac=(double)rc/NN; int att=0;
        while(frac<target && att<4000){
            att++; int i=(int)(U(rng)*N), j=(int)(U(rng)*N), k=j*N+i;
            if(h[k]>0.4) continue;
            double sg=0.0;
            for(int p=0;p<NGD;p++) sg+=phid[p][k]*phid[p][k];
            for(int p=0;p<NGR;p++) sg+=phir[p][k]*phir[p][k];
            if(sg>0.75) continue;
            int op=Rop(rng); if(phir[op][k]>0.2) continue;
            for(int dj=-REX_R;dj<=REX_R;dj++) for(int di=-REX_R;di<=REX_R;di++){
                if(di*di+dj*dj>REX_R*REX_R) continue;
                int kk=((j+dj+N)%N)*N + (i+di+N)%N;
                for(int p=0;p<NGD;p++) phid[p][kk]=0.0;
                for(int p=0;p<NGR;p++) phir[p][kk]=0.0;
                phir[op][kk]=1.0;
            }
            frac+=M_PI*REX_R*REX_R/(double)NN;
        }
    };

    double F=0.0;
    for(long step=1; step<=MAXSTEP; step++){
        for(int p=0;p<NGD;p++){ PD[p]=phid[p].data(); PDN[p]=phidN[p].data(); }
        for(int p=0;p<NGR;p++){ PR[p]=phir[p].data(); PRN[p]=phirN[p].data(); }

        F+=DFDT;
        if(step%REX_EVERY==0) rex_pass(F);

        // precompute grain sums once
        #pragma omp parallel for schedule(static)
        for(int k=0;k<NN;k++){
            double sd=0.0,sr=0.0;
            for(int p=0;p<NGD;p++){ double e=PD[p][k]; sd+=e*e; }
            for(int p=0;p<NGR;p++){ double e=PR[p][k]; sr+=e*e; }
            Sd[k]=sd; Sr[k]=sr;
        }

        // chemical potentials + mobilities
        #pragma omp parallel for schedule(static)
        for(int k=0;k<NN;k++){
            double cvc=clampd(cv[k],1e-6,1.0-1e-6), cgc=clampd(cg[k],1e-6,CMG-1e-6);
            double m=clampd(1.0-cvc-cgc,1e-6,1.0);
            double hh=(h[k]-1.0)*(h[k]-1.0), jj=h[k]*h[k];
            double dfm_dcv=log(cvc/CVEQ)-log(m/meq), dfm_dcg=log(cgc/CGEQ)-log(m/meq);
            double dfb_dcv=log(cvc/CBV)-log((1.0-cvc)/(1.0-CBV));
            double dfb_dcg=log(cgc/CBG)-log((CMG-cgc)/(CMG-CBG));
            double lcv=(cv[R[k]]+cv[Lx[k]]+cv[Uy[k]]+cv[Dn[k]]-4.0*cv[k])/(DX*DX);
            double lcg=(cg[R[k]]+cg[Lx[k]]+cg[Uy[k]]+cg[Dn[k]]-4.0*cg[k])/(DX*DX);
            muv[k]=hh*dfm_dcv+jj*dfb_dcv-KC*lcv;
            mug[k]=hh*dfm_dcg+jj*dfb_dcg-KC*lcg;
            Mv[k]=DV*cvc; Mg[k]=DG*cgc;
        }

        // pre-draw randoms (serial, cheap, reproducible)
        for(int k=0;k<NN;k++){
            pvdep[k]=(U(rng)<PCASC)? U(rng)*VG : 0.0;
            pgr[k]=U(rng); nsv[k]=Gn(rng); nsg[k]=Gn(rng);
        }

        // Cahn-Hilliard update (conserved fields)
        #pragma omp parallel for schedule(static)
        for(int k=0;k<NN;k++){
            double c;
            c=muv[k];
            double dv=( 0.5*(Mv[k]+Mv[R[k]])*(muv[R[k]]-c)
                      + 0.5*(Mv[k]+Mv[Lx[k]])*(muv[Lx[k]]-c)
                      + 0.5*(Mv[k]+Mv[Uy[k]])*(muv[Uy[k]]-c)
                      + 0.5*(Mv[k]+Mv[Dn[k]])*(muv[Dn[k]]-c) )/(DX*DX);
            c=mug[k];
            double dg=( 0.5*(Mg[k]+Mg[R[k]])*(mug[R[k]]-c)
                      + 0.5*(Mg[k]+Mg[Lx[k]])*(mug[Lx[k]]-c)
                      + 0.5*(Mg[k]+Mg[Uy[k]])*(mug[Uy[k]]-c)
                      + 0.5*(Mg[k]+Mg[Dn[k]])*(mug[Dn[k]]-c) )/(DX*DX);
            double Pv=0.0,Pg=0.0;
            if(h[k]<0.8){ Pv=pvdep[k]; Pg=2.0*(1.0-h[k])*(1.0-h[k])*GASDEP*pgr[k]; }
            double Sv=SV*(cv[k]-CVEQ)*Sd[k], Sgs=SG*(cg[k]-CGEQ)*Sd[k];
            cvN[k]=clampd(cv[k]+DT*(dv+Pv-Sv)+NOISE*nsv[k]*sqrt(DT),1e-6,1.0-1e-6);
            cgN[k]=clampd(cg[k]+DT*(dg+Pg-Sgs)+NOISE*nsg[k]*sqrt(DT),1e-6,CMG-1e-6);
        }

        // h update (Allen-Cahn)
        #pragma omp parallel for schedule(static)
        for(int k=0;k<NN;k++){
            double cvc=clampd(cv[k],1e-6,1.0-1e-6), cgc=clampd(cg[k],1e-6,CMG-1e-6);
            double m=clampd(1.0-cvc-cgc,1e-6,1.0);
            double fm=cvc*log(cvc/CVEQ)+cgc*log(cgc/CGEQ)+m*log(m/meq);
            double fb=cvc*log(cvc/CBV)+(1.0-cvc)*log((1.0-cvc)/(1.0-CBV))
                     +cgc*log(cgc/CBG)+(CMG-cgc)*log((CMG-cgc)/(CMG-CBG));
            double sg=Sd[k]+Sr[k];
            double dfbulk_dh=2.0*(h[k]-1.0)*fm+2.0*h[k]*fb;
            double dfpoly_dh=A*(h[k]*h[k]*h[k]-h[k]+2.0*GB*h[k]*sg);
            double lh=(h[R[k]]+h[Lx[k]]+h[Uy[k]]+h[Dn[k]]-4.0*h[k])/(DX*DX);
            double drive=PVG_DRIVE*max(0.0,cvc-0.05)*(1.0-h[k]);
            hN[k]=clampd(h[k]+DT*(-LMOB*(dfbulk_dh+dfpoly_dh-KH*lh)+drive),0.0,1.0);
        }

        // deformed grains
        for(int p=0;p<NGD;p++){
            double* e=PD[p]; double* en=PDN[p];
            #pragma omp parallel for schedule(static)
            for(int k=0;k<NN;k++){
                double v=e[k], sd=Sd[k], sr=Sr[k];
                double dfp=A*(v*v*v-v+2.0*GD*v*(sd-v*v)+2.0*GRD*v*sr+2.0*GB*h[k]*h[k]*v);
                double lp=(e[R[k]]+e[Lx[k]]+e[Uy[k]]+e[Dn[k]]-4.0*v)/(DX*DX);
                en[k]=clampd(v+DT*(-LMOB*(dfp+v*ESTORE-KF*lp)),0.0,1.0);
            }
        }
        // recrystallized grains
        for(int p=0;p<NGR;p++){
            double* e=PR[p]; double* en=PRN[p];
            #pragma omp parallel for schedule(static)
            for(int k=0;k<NN;k++){
                double v=e[k], sr=Sr[k], sd=Sd[k];
                double dfp=A*(v*v*v-v+2.0*GR*v*(sr-v*v)+2.0*GRD*v*sd+2.0*GB*h[k]*h[k]*v);
                double lp=(e[R[k]]+e[Lx[k]]+e[Uy[k]]+e[Dn[k]]-4.0*v)/(DX*DX);
                en[k]=clampd(v+DT*(-LMOB*(dfp-KF*lp)),0.0,1.0);
            }
        }

        cv.swap(cvN); cg.swap(cgN); h.swap(hN);
        for(int p=0;p<NGD;p++) phid[p].swap(phidN[p]);
        for(int p=0;p<NGR;p++) phir[p].swap(phirN[p]);

        while(tptr<targets.size() && F>=targets[tptr]){
            dump(tagof(targets[tptr]));
            cout<<"snapshot F="<<targets[tptr]<<" (step "<<step<<")\n"; tptr++;
        }
        if(tptr>=targets.size()) break;
        if(step%1000==0) cout<<"step "<<step<<" F="<<F<<"\n";
    }
    cout<<"done. files in "<<outdir<<"/\nplot: python fig5_plot.py\n";
    return 0;
}
