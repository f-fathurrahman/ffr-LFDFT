from __future__ import print_function
from ase.units import Ry
import sys

def read_etot(logfile):
    f = open(logfile, 'r')
    while True:
        line = f.readline()
        if not line:
            break
        if 'Total    =' in line:
            Etot = float( line.split()[2] )*Ry
    #print('Total energy = %18.10f eV\n' % (Etot))
    return Etot
    f.close()
