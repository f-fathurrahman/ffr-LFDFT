import numpy as np
import matplotlib.pyplot as plt

plt.clf()
for i in range(1,6):
    dat = np.loadtxt('fort.' + str(10+i))
    plt.plot( dat[:,0], dat[:,1], label='$\phi_'+str(i) + '$')
plt.legend()
plt.grid()
plt.savefig('plotN_5.pdf')
