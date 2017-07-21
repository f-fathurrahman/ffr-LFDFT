import numpy as np
import os
from ase.units import Bohr

start_pos = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])*Bohr
dx = 0.005
Nmoves = 20

f = open('TEMPLATE_INP', 'r')
lines = f.readlines()
f.close()

pos = start_pos[:]
for i in range(Nmoves):
    #
    infile  = 'INPUT_' + str(i+1)
    outfile = 'LOG_' + str(i+1)
    #
    f = open(infile, 'w')
    f.writelines(lines)
    pos[:,0] = pos[:,0] + i*dx
    f.write('H  %18.10f %18.10f %18.10f\n' % (pos[0,0],pos[0,1],pos[0,2]))
    f.write('H  %18.10f %18.10f %18.10f\n' % (pos[1,0],pos[1,1],pos[1,2]))
    f.close()
    #
    os.system('../../ffr_LFDFT_gfortran.x ' + infile + ' | tee ' + outfile)


