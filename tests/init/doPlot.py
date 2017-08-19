import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib import rc
rc('font',**{'family':'serif', 'size':16})
rc('text', usetex=True)

types = ['sinc', 'c', 'p']
NBASIS = 5

for t in types:
    dat1 = np.loadtxt('LF1d_' + t + '.dat')
    plt.clf()
    for ibf in range(1,NBASIS+1):
        plt.plot( dat1[:,0], dat1[:,ibf] )
    plt.xlim( dat1[-1,0], dat1[0,0] )
    plt.grid()
    FIPLOT = 'LF1d_' + t + '.pdf'
    plt.savefig(FIPLOT)
    #
    os.system('pdfcrop ' + FIPLOT + ' ' + FIPLOT)
