import matplotlib.pyplot as plt
import numpy as np

def plot_pseudo(atm, lmax, rc, rlocmax=5.0):

  Vlocal = np.loadtxt(atm + "_Vlocal.dat")
  plt.clf()
  plt.plot( Vlocal[:,0], Vlocal[:,1], marker="o")
  plt.grid()
  plt.title("Local potential " + atm)
  plt.xlim( 0, rlocmax)
  plt.savefig(atm + "_Vlocal.png", dpi=300)

  for l in range(lmax+1):
    proj = np.loadtxt(atm + "_proj_" + str(l) + ".dat")
    plt.clf()
    plt.plot( proj[:,0], proj[:,1], label="prj-1", marker="o" )
    plt.plot( proj[:,0], proj[:,2], label="prj-2", marker="o" )
    plt.plot( proj[:,0], proj[:,3], label="prj-3", marker="o" )
    plt.xlim( 0, rc[l] )
    plt.legend()
    plt.grid()
    plt.title( "Projector for " + atm + " l = " + str(l) )
    plt.savefig( atm + "_proj_" + str(l) + ".png", dpi=300 )


#plot_pseudo( "C", 1, [1.7784266991, 1.3988746415] )
#plot_pseudo( "N", 1, [1.4854009530, 1.6253251935] )
#plot_pseudo( "Ni", 2, [2.4007648869, 3.3395438289, 1.7258544251] )
#plot_pseudo( "Ge", 2, [2.7893555880, 3.4412602616, 4.6453554168] )
plot_pseudo( "H", -1, [] )

