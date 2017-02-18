#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

NBASIS = 4
plt.clf()

for i in range(NBASIS):
  bf = np.loadtxt('fort.' + str(101 + i) )
  lb1 = '$L_{' + str(i+1) + '}(x)$'
  plt.plot( bf[:,0], bf[:,1], linewidth=2, label=lb1 )

plt.legend()
plt.grid()
plt.savefig('LF_N_4.pdf')
plt.savefig('LF_N_4.png',dpi=300)
