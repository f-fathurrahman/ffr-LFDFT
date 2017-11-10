import numpy as np
import matplotlib.pyplot as plt


plt.clf()


for rcut in ['1.0', '1.1', '1.2', '1.3', '1.4', '1.5']:
  dat = np.loadtxt('Ene_N_35_rcut_' + rcut + '.dat')
  plt.plot( dat[:,0], dat[:,1], marker='o', label='rcut='+rcut )

dat = np.loadtxt('Ene_N_35_rcut_1.0.dat')
plt.plot( dat[:,0], dat[:,2], marker='^', color='k', label='standard' )

plt.grid()
plt.legend()
plt.savefig('fort.plot.png', dpi=300)
