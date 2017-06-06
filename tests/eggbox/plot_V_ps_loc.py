import numpy as np
import matplotlib.pyplot as plt
import sys

def do_print(N):
    Vpot = np.loadtxt('fort.' + str(N) )
    plt.clf()
    plt.plot( Vpot[:,0], Vpot[:,1], marker='o', label='V_short' )
    plt.plot( Vpot[:,0], Vpot[:,2], marker='o', label='V_short ss' )
    plt.grid()
    plt.xlim(3.0,13.0)
    plt.legend()
    plt.savefig('fort.' + str(N) + '.png', dpi=300)

do_print(45)

