"""csr_matrix((data, indices, indptr), [shape=(M, N)])
    is the standard CSR representation where the column indices for row i are stored
    in indices[indptr[i]:indptr[i+1]] and their corresponding values are stored in
    data[indptr[i]:indptr[i+1]]. If the shape parameter is not supplied, the matrix
    dimensions are inferred from the index arrays.
"""

from scipy.sparse import *
import numpy as np
import matplotlib.pyplot as plt

values = np.loadtxt('fort.100')
column = np.loadtxt('fort.200',dtype=np.int64)
idxRow = np.loadtxt('fort.300',dtype=np.int64)

# convert to zero-based indexing
column = column - 1
idxRow = idxRow - 1

N = 8**3
nabla2 = csr_matrix( (values, column, idxRow), shape=(N,N) )

#print nabla2.todense()

nabla2_D = nabla2.todense()
for i in range(0,N):
  nabla2_D[i,i] = 1.0
  for j in range(i+1,N):
    if( abs( nabla2_D[i,j] ) > 0.0 ):
      nabla2_D[i,j] = 1.0
    nabla2_D[j,i] = nabla2_D[i,j]

plt.clf()
#plt.spy(nabla2, precision=0.1, markersize=5)
#plt.spy(nabla2,precision=2.0)
plt.matshow(-nabla2_D)
plt.set_cmap('gray')
plt.savefig('nabla2_8_8.pdf')





