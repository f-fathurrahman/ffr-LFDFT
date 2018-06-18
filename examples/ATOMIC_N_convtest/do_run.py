from __future__ import print_function
import numpy as np
import os
from ase.units import Bohr
import subprocess

def write_control(f):
    str_control = """&CONTROL
  pseudo_dir = '../../HGH'
  etot_conv_thr = 1.0d-6
/
"""
    f.write(str_control)


# FIXME: this will only work for nat=1, nsp=1
def write_system(f, nr1=45, nr2=45, nr3=45):
    str_system = """
&SYSTEM
  ibrav = 8
  nat = 1
  ntyp = 1
  A = 8.4668d0
  B = 8.4668d0
  C = 8.4668d0
"""
    f.write(str_system)
    f.write('  nr1 = %d\n' % nr1)
    f.write('  nr2 = %d\n' % nr2)
    f.write('  nr3 = %d\n' % nr3)
    f.write('/\n')


def write_electron(f):
    str_electron = """
&ELECTRONS
  KS_Solve = 'Emin_pcg'
  cg_beta = 'PR'
  electron_maxstep = 150
  diagonalization = 'LOBPCG'
/
"""
    f.write(str_electron)

# FIXME: this will only work for nat=1, nsp=1
def write_atomic_species(f, species, pseudo):
    f.write('\nATOMIC_SPECIES\n')
    f.write('%s  1.0  %s\n' % (species,pseudo))

# FIXME: this will only work for nat=1, nsp=1
#        fixed at 0,0,0
def write_atomic_positions(f, species):
    str_atomic_positions = """
ATOMIC_POSITIONS angstrom
"""
    f.write(str_atomic_positions)
    f.write('%s  0.0   0.0   0.0\n\n' % (species))


EXE_PATH = '../../src/ffr_LFDFT_gfortran.x '

#SPECIES = ['H', 'He']
SPECIES = ['H']

#N = [45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71]
N = [35, 37]


for s in SPECIES:
    fout = open('Etot_' + s + '.dat', 'w')
    for n in N:
        infile = 'INPUT_' + s + '_' + str(n)
        f = open(infile,'w')
        write_control(f)
        write_system(f, nr1=n, nr2=n, nr3=n)
        write_electron(f)
        write_atomic_species(f, s, s + '.hgh')
        write_atomic_positions(f, s)
        f.close()
        #
        outfile = 'OUTPUT_' + s + '_' + str(n)
        #
        os.system(EXE_PATH + infile + ' > ' + outfile)
        #
        line = subprocess.Popen('grep "! Total" ' + outfile, shell=True, stdout=subprocess.PIPE).stdout.read()
        Etot = float( line.split()[3] )
        #
        fout.write('%d  %18.10f\n' % (n, Etot))
    #
    fout.close()


















