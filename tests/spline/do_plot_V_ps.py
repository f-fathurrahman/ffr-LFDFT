import numpy as np
import matplotlib.pyplot as plt

V_ps_x = np.loadtxt('fort.101')

plt.clf()
plt.plot( V_ps_x[:,0], V_ps_x[:,1], marker='o', label='orig' )
plt.savefig('V_ps_x.png', dpi=300)

