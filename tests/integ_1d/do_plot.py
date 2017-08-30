import numpy as np
import matplotlib.pyplot as plt
from ase.units import Hartree

for N in range(15,75+1,4):
    filnam = 'res_' + str(N) + '.dat'
    filout = 'fort.res_' + str(N) + '.png'
    dat = np.loadtxt(filnam)
    dat[:,1] = dat[:,1]*Hartree
    Emin = np.max(dat[:,1])
    Emax = np.min(dat[:,1])
    Emid = (Emin + Emax)/2
    dat[:,1] = dat[:,1] - Emid
    Emin = Emin - Emid
    Emax = Emax - Emid
    plt.clf()
    plt.plot( dat[:,0], dat[:,1], marker='o' )
    plt.grid()
    plt.ylim(Emin,Emax)
    plt.savefig(filout,dpi=300)



