import numpy as np
import matplotlib.pyplot as plt

def do_print(N):
    Vpot = np.loadtxt('fort.' + str(N) )
    delta = 0.5*( Vpot[1,0] - Vpot[0,0] )
    plt.clf()
    plt.plot( Vpot[:,0], Vpot[:,1], marker='o', label='V_ps' )
    #plt.plot( Vpot[:,0], Vpot[:,2], marker='o', label='V_Ha' )
    plt.grid()
    plt.savefig('Be_V_ps_loc_' + str(N) + '.png', dpi=300)

#do_print(25)
#do_print(35)
#do_print(45)
#do_print(55)
#do_print(65)
do_print(75)

