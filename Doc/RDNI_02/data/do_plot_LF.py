import numpy as np
import matplotlib.pyplot as plt

FILEPLOT = '../images/plotN_5.pdf'

from matplotlib import rc
rc('font',**{'family':'serif', 'size':16})
rc('text', usetex=True)

plt.clf()
for i in range(1,6):
    dat = np.loadtxt('LF5_i_' + str(i) + '.dat')
    plt.plot( dat[:,0], dat[:,1], label='$L_'+str(i) + '$')
plt.legend()
plt.grid()
plt.savefig(FILEPLOT)

import os
os.system('pdfcrop ' + FILEPLOT + ' ' + FILEPLOT)
