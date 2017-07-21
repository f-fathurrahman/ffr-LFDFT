from __future__ import print_function
import numpy as np
from ase.units import Bohr
from ase.units import Ry
import sys
import os

def read_etot(logfile):
    f = open(logfile, 'r')
    while True:
        line = f.readline()
        if not line:
            break
        if 'Total    =' in line:
            Etot = float( line.split()[2] )*Ry
    return Etot
    f.close()


start_pos = 0.0
dx = 16.0*Bohr/45/20
Nmoves = 20

f = open('TEMPLATE_INP', 'r')
lines = f.readlines()
f.close()

pos = start_pos
x = []
Etot = []

fdat = open('ETOT.dat', 'w')
fdat.write('# dx = %18.10f , Nmoves = %d\n' % (dx, Nmoves))

for i in range(Nmoves):
    #
    infile  = 'INPUT_' + str(i+1)
    outfile = 'LOG_' + str(i+1)
    #
    f = open(infile, 'w')
    f.writelines(lines)
    pos = pos + dx
    print('pos = ', pos)
    f.write('H  %18.10f %18.10f %18.10f\n' % (pos,0.0,0.0))
    f.close()
    #
    os.system('../../../ffr_LFDFT_gfortran.x ' + infile + ' > ' + outfile)
    #
    x.append(pos)
    Etot.append(read_etot(outfile))
    fdat.write('%18.10f %18.10f\n' % (x[-1], Etot[-1]))

fdat.close()

import matplotlib.pyplot as plt
plt.clf()
plt.plot( x, Etot, marker='o')
plt.savefig('dx_Etot.png', dpi=300)
