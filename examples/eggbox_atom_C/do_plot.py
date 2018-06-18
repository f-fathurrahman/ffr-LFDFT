import numpy as np
import matplotlib.pyplot as plt
import sys

N = sys.argv[1]

plt.clf()

dat = np.loadtxt('Etot.' + N + '.dat')
Emin = np.min( dat[:,1] )
Emax = np.max( dat[:,1] )
Emid = 0.5*(Emin + Emax)
delta = Emax-Emin

plt.plot( dat[:,0], (dat[:,1]-Emid)*1000, marker='o', label='Etot' )
plt.xlabel('Position (bohr)')
plt.ylabel('Diff E (mHa)')

print('std      Emid = %18.10f, delta = %18.10f mHa' % (Emid,delta*1000))

plt.grid()
plt.legend()
plt.savefig('fort.Etot.' + N + '.pdf')


