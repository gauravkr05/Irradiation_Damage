from interface_energies import *
# GB energy with gas allowed to segregate (grand-canonical at mu = 0), plus GB-core state
print(" bGB  gGB(c=0)  gGB(c relaxed)  S_core=n1^2+n2^2  c_core   int16n1^2n2^2dz")
for bGB in (1.1,1.2,1.3,1.5,1.6,2.0,2.5,3.0):
    P=dict(BASE,bGB=bGB)
    E0,w,(z,c,a,b)=gb_energy(P,1.0,L=60)
    # now relax c too, starting from a small seed
    c=0.05*np.exp(-z**2/4); 
    E1,_=relax(c,a.copy(),b.copy(),P,1.0,{"c","a","b"})
    a2=a; i=np.argmin(abs(z))
    E1b,w2,(z2,c2,a2,b2)=None,None,(z,c,a,b)
    S=(a**2+b**2); 
    print(f" {bGB:4.2f}  {E0:.4f}    {E1:.4f}          {S.min():.3f}          {c.max():.3f}    {np.sum(16*a**2*b**2):.3f}")
