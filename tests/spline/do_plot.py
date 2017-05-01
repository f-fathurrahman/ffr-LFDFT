import numpy as np
import matplotlib.pyplot as plt

dat_orig = np.loadtxt('fort.201')
dat_spli = np.loadtxt('fort.202')

plt.clf()
plt.plot( dat_orig[:,0], dat_orig[:,1], marker='o', label='orig' )
plt.plot( dat_spli[:,0], dat_spli[:,1], marker='x', label='spli' )
plt.savefig('compare.png', dpi=300)

