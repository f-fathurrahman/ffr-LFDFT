#!/usr/bin/python

import h5py
import matplotlib.pyplot as plt
import sys

print sys.argv[1]
f1 = h5py.File( sys.argv[1] )

k0 = f1.keys()[0]
mat = f1[k0]
#plt.imshow(mat)
plt.matshow(mat)
plt.colorbar()
#plt.savefig(sys.argv[1] + '.png', dpi=300)
plt.savefig(sys.argv[1] + '.pdf')

