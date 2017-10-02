import numpy as np
import matplotlib.pyplot as plt
import sys

plt.style.use('dark_background')

def do_print(N):
    dat = np.loadtxt('fort.' + str(N) )
    plt.clf()
    Ncols = dat.shape[1]
    for ic in range(1,Ncols):
        plt.plot( dat[:,0], dat[:,ic], marker='o', label='col-'+str(ic) )
    plt.grid()
    plt.xlim(7.0,9.0)
    plt.legend()
    plt.savefig('fort.' + str(N) + '.png', dpi=300)

do_print(int(sys.argv[1]))

