import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

dat = np.loadtxt('fort.101')
N1 = int( sqrt( dat.shape[0] ) )
N2 = int( dat.shape[1] )
print(N1)
print(N2)

x = np.reshape( dat[:,0], (N1,N1) )
y = np.reshape( dat[:,1], (N1,N1) )
z = np.reshape( dat[:,2], (N1,N1) )

if N2 == 4:
    z2 = np.reshape( dat[:,3], (N1,N1) )

from mpl_toolkits.mplot3d import Axes3D

ax = plt.figure().gca(projection='3d')
# rstride=cstride=1 is important for Fortran based indexing data
#ax.plot_surface( x, y, z, linewidth=0, rstride=1, cstride=1, cmap='jet' )
if N2 == 4:
    ax.plot_surface( x, y, z2, linewidth=0, rstride=1, cstride=1, cmap='jet' )

#ax.set_xlim(6,10)
#ax.set_ylim(6,10)
plt.show()


