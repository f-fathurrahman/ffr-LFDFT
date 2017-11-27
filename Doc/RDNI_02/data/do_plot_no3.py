import numpy as np
import matplotlib.pyplot as plt
import os

PRE = '../images/'
BASEFILE = 'A_10_alpha_3.pdf'

from matplotlib import rc
rc('font',**{'family':'serif', 'size':16})
rc('text', usetex=True)

dat1 = np.loadtxt('A_10_alpha_3_oct.dat')
dat2 = np.loadtxt('A_10_alpha_3_LF.dat')

FILEPLOT = PRE + BASEFILE

plt.clf()
plt.plot( dat1[:,0], dat1[:,1], marker='o', linewidth=2, label='OCTOPUS' )
plt.plot( dat2[:,0], dat2[:,1], marker='o', linewidth=2, label='LFDFT' )
plt.grid()
plt.legend(loc='lower left')
plt.xlabel('Grid spacing (bohr)')
plt.ylabel('Electronic energy (Ha)')
plt.savefig(FILEPLOT)

os.system('pdfcrop ' + FILEPLOT + ' ' + FILEPLOT)

dEne_OCT = np.abs( dat1[:,1] - dat1[-1,1] )
dEne_LF  = np.abs( dat2[:,1] - dat2[-1,1] )

idx_end = len(dEne_OCT)

CONV_FILEPLOT = PRE + 'CONV_' + BASEFILE

plt.clf()
plt.plot( dat1[0:idx_end-1,0], np.log10(dEne_OCT[0:idx_end-1]), marker='o', linewidth=2, label='OCTOPUS' )
plt.plot( dat2[0:idx_end-1,0], np.log10(dEne_LF[0:idx_end-1]), marker='o', linewidth=2, label='LFDFT' )
plt.grid()
plt.legend(loc='upper left')
plt.xlabel('Grid spacing (bohr)')
plt.ylabel('Convergence (log(Ha))')
plt.text(0.5, -2.5, "$\\alpha=3.0$", fontsize=24)
plt.yticks( [1.0, 0.0, -1.0, -2.0, -3.0] )
plt.savefig(CONV_FILEPLOT)

os.system('pdfcrop ' + CONV_FILEPLOT + ' ' + CONV_FILEPLOT)
