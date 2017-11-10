from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

plt.clf()
fig = plt.gcf()
fig.set_size_inches(10.0, 10.0)

dat1 = np.loadtxt('fort.100')
dat2 = np.loadtxt('fort.101')
dat3 = np.loadtxt('fort.102')
dat4 = np.loadtxt('fort.103')

Npts1 = dat1.shape[0]
Npts2 = dat2.shape[0]
Npts3 = dat3.shape[0]
Npts4 = dat4.shape[0]

# coarse grids
for i in range(Npts1):
    plt.plot( dat1[i,0], dat1[i,1], marker='o', color='b' )

for i in range(Npts2):
    plt.plot( dat2[i,0], dat2[i,1], marker='o', color='r' )

# fine grids
for i in range(Npts3):
    plt.plot( dat3[i,0], dat3[i,1], marker='^', color='b' )

for i in range(Npts4):
    plt.plot( dat4[i,0], dat4[i,1], marker='^', color='r' )

# center
c = np.loadtxt('fort.99')
plt.plot( c[0], c[1], marker='*', color='k' )
c1 = plt.Circle( (c[0],c[1]), radius=c[2], color='k', alpha=0.5 )
plt.gca().add_artist( c1 )

plt.axis('equal')
plt.xlim(0.0,16)
plt.ylim(0.0,16)
plt.grid()
plt.savefig('fort.gridpts.pdf')
