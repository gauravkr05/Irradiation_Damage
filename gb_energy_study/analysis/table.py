from interface_energies import *
import json
res={}
# bubble/grain: both lattice positions + continuum
EpA,_,_=pg_energy(dict(BASE,bGB=1.0),1.0,z0=0.0)
EpB,_,_=pg_energy(dict(BASE,bGB=1.0),1.0,z0=0.5)
Epc,wp,_=pg_energy(dict(BASE,bGB=1.0),0.05)
print(f"gamma_PG lattice node/bond {EpA:.4f}/{EpB:.4f}, continuum {Epc:.4f}, width {wp:.2f}")
gPG=0.5*(EpA+EpB)
res['gPG_lattice']=gPG; res['gPG_cont']=Epc
print(" bGB  gGB_mid  gGB_node  gGB_cont  width  S_core  int16  thetaY_lat  thetaY_cont")
for bGB in (1.0,1.2,1.5,2.0,2.5,3.0):
    P=dict(BASE,bGB=bGB)
    Em,w,(z,c,a,b)=gb_energy(P,1.0,L=60,z0=0.5)
    En,_,_=gb_energy(P,1.0,L=60,z0=0.0)
    Ec,_,(zc,cc,ac,bc)=gb_energy(P,0.05,L=60)
    S=(a**2+b**2).min(); I=np.sum(16*a**2*b**2)
    tl=theta_young(Em,gPG); tc=theta_young(Ec,Epc)
    print(f" {bGB:4.2f}  {Em:.4f}   {En:.4f}   {Ec:.4f}   {w:5.2f}  {S:.3f}  {I:6.2f}   {tl:6.1f}   {tc:6.1f}")
    res[str(bGB)]=dict(gGB_mid=Em,gGB_node=En,gGB_cont=Ec,width=w,S=S,int16=I,thetaY_lat=tl,thetaY_cont=tc)
json.dump(res,open('interface_table.json','w'),indent=1)
