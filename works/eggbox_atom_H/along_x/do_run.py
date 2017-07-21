import numpy as np
import os
from ase.units import Bohr
from read_etot import read_etot

start_pos = np.array([0.0, 0.0, 0.0])
dx = 16.0*Bohr/45/40
Nmoves = 20

f = open('TEMPLATE_INP', 'r')
lines = f.readlines()
f.close()

pos = start_pos[:]
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
    pos[0] = pos[0] + i*dx
    f.write('H  %18.10f %18.10f %18.10f\n' % (pos[0],pos[1],pos[2]))
    f.close()
    #
    os.system('../../../ffr_LFDFT_gfortran.x ' + infile + ' | tee ' + outfile)
    #
    x.append(i*dx)
    Etot.append(read_etot(outfile))
    fdat.write('%18.10f %18.10f\n' % (x[-1], Etot[-1]))

fdat.close()

import matplotlib.pyplot as plt
plt.clf()
plt.plot( x, Etot, marker='o')
plt.savefig('dx_Etot.png', dpi=300)
