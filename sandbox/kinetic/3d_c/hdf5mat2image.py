#!/usr/bin/python

import h5py
import matplotlib.pyplot as plt
import sys
import numpy as np

print sys.argv[1]
f1 = h5py.File( sys.argv[1] )

k0 = f1.keys()[0]
mat = np.abs( np.array(f1[k0]) )
#print mat.shape
Nrow = mat.shape[0]
Ncol = mat.shape[1]

# Convert to "Binary"
Nonzero = 0
for i in range(Nrow):
  for j in range(Ncol):
    if( abs(mat[i][j]) > 1.e-13 ):
      #print '%d %d is not zero:' % (i,j)
      Nonzero = Nonzero + 1
      mat[i][j] = 1.0

print 'Nonzero = ', Nonzero
print 'Percentage of Nonzero = ', float(Nonzero)/(Nrow*Ncol)*100.0
#plt.imshow(mat)
plt.matshow(-mat) #
plt.set_cmap('gray')
#plt.colorbar()
#plt.savefig(sys.argv[1] + '.png', dpi=300)
plt.savefig(sys.argv[1] + '.pdf')

