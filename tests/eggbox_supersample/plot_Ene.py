import numpy as np
import matplotlib.pyplot as plt

plt.clf()

dat = np.loadtxt('Ene_N_35_std.dat')
# middle energy is calculated here
Emin = np.min( dat[:,1] )
Emax = np.max( dat[:,1] )
Emid = 0.5*(Emin + Emax)
#
plt.plot( dat[:,0], dat[:,2]-Emid, marker='^', color='k', label='standard' )

for rcut in ['1.0', '1.1', '1.2', '1.3', '1.4', '1.5']:
  dat = np.loadtxt('Ene_N_35_rcut_' + rcut + '.dat')
  plt.plot( dat[:,0], dat[:,1]-Emid, marker='o', label='rcut='+rcut )

plt.grid()
plt.legend()
plt.savefig('fort.plot.pdf')

