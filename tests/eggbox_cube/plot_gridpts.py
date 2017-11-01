from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

dat1 = np.loadtxt('fort.100')
dat2 = np.loadtxt('fort.101')

N1 = int( round( dat1.shape[0]**(1.0/3.0) ) )
N2 = int( round( dat2.shape[0]**(1.0/3.0) ) )
print(N1)
print(N2)

# notice that the indexing works in reverse-order !!!!
# (iz,iy,ix,i)
grid1 = np.reshape( dat1, [N1,N1,N1,3] )
grid2 = np.reshape( dat2, [N2,N2,N2,3] )

zslice = 0
plt.clf()
#
for i in range(N1):
    for j in range(N1):
        x = grid1[zslice,j,i,0]
        y = grid1[zslice,j,i,1]
        plt.plot( x, y, marker='o', color='b' )

for i in range(N2):
    for j in range(N2):
        x = grid2[zslice,j,i,0]
        y = grid2[zslice,j,i,1]
        plt.plot( x, y, marker='x', color='b' )

plt.xlim(0,16)
plt.ylim(0,16)
plt.axis('equal')
plt.grid()
plt.savefig('fort.101.png', dpi=300)
plt.savefig('fort.101.pdf')
