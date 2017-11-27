import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
rc('font',**{'family':'serif'})
rc('text', usetex=True)

dat1 = np.loadtxt('A_10_alpha_1_oct.dat')
dat2 = np.loadtxt('A_10_alpha_1_LF.dat')

plt.clf()
plt.plot( dat1[:,0], dat1[:,1], marker='o', label='OCTOPUS' )
plt.plot( dat2[:,0], dat2[:,1], marker='o', label='LFDFT' )
plt.grid()
plt.legend(loc='lower left')
plt.xlabel('Grid spacing (bohr)')
plt.ylabel('Electronic energy (Ha)')
plt.xlim(0.25,0.6)
plt.ylim(-7.6,-7.4)
plt.savefig('A_10_alpha_1_v2.pdf')
