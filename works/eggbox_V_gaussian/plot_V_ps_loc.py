from __future__ import print_function
import sys
from scipy.special import erf
from math import exp, sqrt
import numpy as np

pspFile = open(sys.argv[1])

pspFile.readline()
line = pspFile.readline()

Nparams = len( line.split() )
if Nparams < 3:
    raise RuntimeError('Too few parameters')

symb = line.split()[0]
Zval = float( line.split()[1] )
rloc = float( line.split()[2] )

C = np.zeros(4)
line1 = line.split()[3:]
NCparams = len(line1)
for i in range(NCparams):
    C[i] = float(line1[i])

r0 = 1e-10
NptsPlot = 500
r = np.linspace(r0,8.0,NptsPlot)
Vloc = np.zeros(NptsPlot)
Vloc_short = np.zeros(NptsPlot)
Vloc_long  = np.zeros(NptsPlot)

for i in range(NptsPlot):
    rrloc = r[i]/rloc
    Vloc_long[i] = -Zval/r[i]*erf(rrloc/sqrt(2.0))
    Vloc_short[i] = exp(-0.5*rrloc**2)*( C[0] + C[1]*rrloc**2 + C[2]*rrloc**4 + C[3]*rloc**6 )

Vloc[:] = Vloc_long[:] + Vloc_short[:]

import matplotlib.pyplot as plt

plt.clf()
plt.plot( r, Vloc, marker='o', linewidth=2, label='Full')
plt.plot( r, Vloc_long, marker='o', linewidth=2, label='Long')
plt.plot( r, Vloc_short, marker='o', linewidth=2, label='Short')
plt.grid()
plt.legend()
filplot = sys.argv[1] + '.png'
plt.savefig(filplot, dpi=300)
