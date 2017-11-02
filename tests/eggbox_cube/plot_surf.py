import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

dat = np.loadtxt('fort.100')
N1 = int( sqrt( dat.shape[0] ) )
print(N1)

x = np.reshape( dat[:,0], (N1,N1) )
y = np.reshape( dat[:,1], (N1,N1) )
z = np.reshape( dat[:,2], (N1,N1) )

print(x)
print(y)
print(z)

from mpl_toolkits.mplot3d import Axes3D

ax = plt.figure().gca(projection='3d')
# rstride=cstride=1 is important for Fortran based indexing data
ax.plot_surface( x, y, z, linewidth=0, rstride=1, cstride=1 )
plt.show()


