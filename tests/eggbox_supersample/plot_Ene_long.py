import numpy as np
import matplotlib.pyplot as plt
import sys

N = sys.argv[1]

plt.clf()

dat = np.loadtxt('Ene_N_' + N + '_long.dat')
Emin = np.min( dat[:,1] )
Emax = np.max( dat[:,1] )
Emid = 0.5*(Emin + Emax)
delta = Emax-Emin

plt.plot( dat[:,0], dat[:,1]-Emid, marker='o', label='long-' + N )

print('std      Emid = %18.10f, delta = %18.10f mHa' % (Emid,delta*1000))

plt.grid()
plt.legend()
plt.savefig('fort.plot.long.' + N + '.pdf')


