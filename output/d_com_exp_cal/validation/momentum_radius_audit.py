"""Independent momentum-space <r^2> using Gauss-Legendre spectral derivatives.
This is validation only. User observables are calculated by production Fortran.
"""
from pathlib import Path
import json,sys
import numpy as np
source=Path(sys.argv[1])
a=np.loadtxt(source,delimiter=",",skiprows=1)
p,w,s,d=a[:,:4].T
pmax=1000.0
x=2*p/pmax-1
# Barycentric weights of roots of P_N. The common scale cancels.
b=(-1.0)**np.arange(len(p))*np.sqrt((1-x*x)*w)
separation=p[:,None]-p[None,:]
np.fill_diagonal(separation,1)
D=b[None,:]/b[:,None]/separation
np.fill_diagonal(D,0)
np.fill_diagonal(D,-D.sum(axis=1))
d_s=D@s;d_d=D@d
r2=197.327**2*np.sum(w*(p*p*(d_s*d_s+d_d*d_d)+6*d*d))
result=dict(source=str(source.resolve()),method='Gauss-Legendre barycentric derivative; r2=hbarc2 integral[dp p2 psi_prime2 + dp L(L+1) psi2]',point_radius_fm=float(np.sqrt(r2)/2),relative_radius_squared_fm2=float(r2))
print(json.dumps(result,indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
