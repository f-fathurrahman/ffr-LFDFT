import numpy as np
import matplotlib.pyplot as plt


plt.clf()

dat = np.loadtxt('Ene_N_35_long.dat')
Emin = np.min( dat[:,1] )
Emax = np.max( dat[:,1] )
Emid = 0.5*(Emin + Emax)
plt.plot( dat[:,0], dat[:,1]-Emid, marker='o', label='long-35' )

plt.grid()
plt.legend()
plt.savefig('fort.plot.long.pdf')
