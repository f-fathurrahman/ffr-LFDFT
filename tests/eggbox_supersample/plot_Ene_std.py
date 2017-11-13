import numpy as np
import matplotlib.pyplot as plt
import sys

N = sys.argv[1]

plt.clf()

dat = np.loadtxt('Ene_N_' + N + '_std.dat')
Emin = np.min( dat[:,1] )
Emax = np.max( dat[:,1] )
Emid = 0.5*(Emin + Emax)
delta = Emax-Emin

plt.plot( dat[:,0], (dat[:,1]-Emid)*1000, marker='o', label='std-' + N )
plt.ylabel('Diff ene (mHa)')

print('std      Emid = %18.10f, delta = %18.10f mHa' % (Emid,delta*1000))

plt.grid()
plt.legend()
plt.savefig('fort.plot.std.' + N + '.pdf')


