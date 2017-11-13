import numpy as np
import matplotlib.pyplot as plt

plt.clf()

dat = np.loadtxt('Ene_N_35_std.dat')
# middle energy is calculated here
Emin = np.min( dat[:,1] )
Emax = np.max( dat[:,1] )
Emid = 0.5*(Emin + Emax)
delta = Emax-Emin
#
plt.plot( dat[:,0], dat[:,1]-Emid, marker='^', color='k', label='standard' )

print('std      Emid = %18.10f, delta = %18.10f mHa' % (Emid,delta*1000))

rcutFull = ['1.0', '1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '1.7', '1.8', '1.9', '2.0']

#for rcut in ['1.2', '1.4', '1.6', '1.8', '2.0']:
for rcut in rcutFull:
    dat = np.loadtxt('Ene_N_35_rcut_' + rcut + '.dat')
    plt.plot( dat[:,0], dat[:,1]-Emid, marker='o', label='rcut='+rcut )
    Emin = np.min(dat[:,1])
    Emax = np.max(dat[:,1])
    print('rcut %s Emid = %18.10f, delta = %18.10f mHa' %  ( rcut, 0.5*(Emin+Emax), (Emax-Emin)*1000 ) )

plt.grid()
plt.legend()
plt.savefig('fort.plot.pdf')

