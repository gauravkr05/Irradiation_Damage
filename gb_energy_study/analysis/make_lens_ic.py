"""Write an initial condition (PFGBIC01 format read by initIn=) containing one
lens-shaped (2-D, x-z) bubble on the GB with prescribed initial dihedral angle.
Area is fixed to pi*R0^2 so all starts contain the same gas in the bubble."""
import numpy as np, struct, sys
NX, NY, NZ = 160, 1, 96
ZGB = NZ//2 - 0.5
def lens(theta_deg, R0=16.0, cmat=0.011, w=1.5, fn="ic.bin"):
    al = np.radians(theta_deg)/2
    rho = np.sqrt(np.pi*R0**2/(2*(al - np.sin(al)*np.cos(al))))
    d = -rho*np.cos(al)                       # cap-centre height relative to GB (top cap)
    x = np.arange(NX)[None,:] - NX/2; z = np.arange(NZ)[:,None] - ZGB
    # signed distance (approx) to lens = max over the two discs (intersection)
    r_top = np.sqrt(x**2 + (z - d)**2) - rho     # top cap disc centred below the GB
    r_bot = np.sqrt(x**2 + (z + d)**2) - rho
    sd = np.maximum(r_top, r_bot)
    b = 0.5*(1 - np.tanh(sd/w))
    g1 = 0.5*(1 - np.tanh(z/w)); g2 = 0.5*(1 + np.tanh(z/w))
    n1 = (g1*(1-b))*np.ones_like(b); n2 = (g2*(1-b))*np.ones_like(b)
    c = cmat + (1-cmat)*b
    with open(fn,'wb') as f:
        f.write(b'PFGBIC01'); f.write(struct.pack('3i',NX,NY,NZ))
        f.write(struct.pack('4d',R0,cmat,0.0,0.0)); f.write(struct.pack('i',1))
        f.write(struct.pack('2d',NX/2,0.0))
        for a in (c,n1,n2): f.write(np.ascontiguousarray(a.reshape(NZ,NY,NX),dtype=np.float64).tobytes())
    print(fn, f"theta_init={theta_deg} rho={rho:.2f} half-width={rho*np.sin(al):.2f} gas={c.sum():.3f}")
lens(179.9, fn="ic_circle.bin")
lens(80.0,  fn="ic_flat.bin")
