from interface_energies import *
gPG_lat = 1.0634; gPG_con = 1.0959
print("bGB=1.0 (baseline): GB energy vs box length (no barrier -> GB spreads)")
for L in (10,20,40,80):
    E,w,_ = gb_energy(dict(BASE,bGB=1.0),1.0,L=L)
    print(f"  L={L:3d}  gamma_GB={E:.4f}  width={w:.2f}")
print()
print("Independent control: j=2, jGB=2, bS=1 fixed; vary bGB only")
print(" bGB  gGB_lat  w_lat  gGB_con  theta_Young_lat  theta_Young_con")
for bGB in (1.05,1.1,1.2,1.3,1.5,1.6,2.0,2.5,3.0,4.0,6.0,10.0):
    P=dict(BASE,bGB=bGB)
    El,wl,_=gb_energy(P,1.0,L=60); Ec,wc,_=gb_energy(P,0.05,L=60)
    print(f" {bGB:4.2f}  {El:.4f}  {wl:5.2f}  {Ec:.4f}  {theta_young(El,gPG_lat):6.1f}  {theta_young(Ec,gPG_con):6.1f}")
print()
print("Legacy setTheta coupling (jGB = 2 bGB, bS from table)")
for th,bS,bGB in ((40,1.0,1.6),(50,1.0,1.3),(60,1.0,1.0),(80,1.7,0.7)):
    P=dict(A0=4.0,bS=bS,bGB=bGB,j=2*bS,jGB=2*bGB)
    E,w,_=gb_energy(P,1.0,L=60); Ep,wp,_=pg_energy(P,1.0)
    print(f" label {th}: bS={bS} bGB={bGB} jGB={2*bGB}  gGB={E:.4f} (w={w:.1f})  gPG={Ep:.4f}  Young theta={theta_young(E,Ep):.1f}")
