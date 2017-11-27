import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib import rc
rc('font',**{'family':'serif', 'size':16})
rc('text', usetex=True)

A = 10.0
center = 8.0
x = np.linspace(0.0, 16.0, 200)
y0 = -A*np.exp(-0.5* (x-center)**2)
y1 = -A*np.exp(-1.0* (x-center)**2)
y2 = -A*np.exp(-2.0* (x-center)**2)
y3 = -A*np.exp(-3.0* (x-center)**2)

FILEPLOT = '../images/V_gauss.pdf'
plt.clf()
plt.plot(x, y0, linewidth=2, label='$\\alpha=0.5$')
plt.plot(x, y1, linewidth=2, label='$\\alpha=1.0$')
plt.plot(x, y2, linewidth=2, label='$\\alpha=2.0$')
plt.plot(x, y3, linewidth=2, label='$\\alpha=3.0$')
plt.grid()
plt.legend()
plt.xlabel('Position (bohr)')
plt.ylabel('Potential')
plt.savefig(FILEPLOT)

os.system('pdfcrop ' + FILEPLOT + ' ' + FILEPLOT)
