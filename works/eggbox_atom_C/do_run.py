from __future__ import print_function
import numpy as np
import os
import subprocess

start_pos = np.array([0.0, 0.0, 0.0])
dx = 0.005
Nmoves = 20

f = open('TEMPLATE_INP', 'r')
lines = f.readlines()
f.close()

# executable path, don't forget to append space
EXEPATH = '../../src/ffr_LFDFT_gfortran.x '

fdat = open('Etot.51.dat', 'w')

pos = start_pos[:]
for i in range(Nmoves):
    #
    infile  = 'INPUT_' + str(i+1)
    outfile = 'LOG_' + str(i+1)
    #
    f = open(infile, 'w')
    f.writelines(lines)
    pos[:] = pos[:] + i*dx
    # write the atoms here
    f.write('C  %18.10f %18.10f %18.10f\n' % (pos[0],pos[1],pos[2]))
    f.close()
    #
    os.system(EXEPATH + infile + ' | tee ' + outfile)
    #
    str1 = subprocess.Popen('grep "! Total" ' + outfile, shell=True, stdout=subprocess.PIPE).stdout.read()
    Etot = float( str1.split()[3] )
    fdat.write('%18.10f %18.10f\n' % ( pos[0], Etot ))

fdat.close()



