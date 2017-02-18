import numpy as np
import matplotlib.pyplot as plt

pot = np.loadtxt('fort.11')

def plot_states( istates ):
  psi_c = np.loadtxt('fort.' + str(100+istates))  # coarse, coeffs
  psi_d = np.loadtxt('fort.' + str(200+istates))  # dense, analytic formula
  plt.plot( psi_c[:,0], psi_c[:,1], marker='o', label='coef' )
  plt.plot( psi_d[:,0], psi_d[:,1], label='ana' )


def plot_states_dense( istates ):
  psi_d = np.loadtxt('fort.' + str(200+istates))
  plt.plot( psi_d[:,0], psi_d[:,1], linewidth=2, label='st-'+str(istates-1) )

"""
plt.clf()
plot_states_dense(1)
plot_states_dense(2)
plot_states_dense(3)
plot_states_dense(4)
plt.legend()
plt.grid()
plt.savefig('wfs.pdf')
"""

for i in range(1,5):
  plt.clf()
  plot_states(i)
  plt.grid()
  plt.savefig('wf' + str(i) + '.pdf')

