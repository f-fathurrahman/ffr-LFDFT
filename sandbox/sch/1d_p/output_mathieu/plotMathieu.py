import numpy as np
import matplotlib.pyplot as plt

pot = np.loadtxt('V_mathieu')
Ene = np.loadtxt('ene_k_0.2')
evecr = np.loadtxt('evec_real.dat')
evecc = np.loadtxt('evec_cplx.dat')

Nstates = len(Ene)
xgrid = evecr[:,0]
Xmin = np.min(xgrid)
Xmax = np.max(xgrid)
Nbasis = evecr.shape[0]
rhoik = np.zeros( (Nbasis,Nstates) )
for ii in range(Nstates):
  rhoik[:,ii] = np.sqrt(evecr[:,ii+1]**2 + evecc[:,ii+1]**2)
#
plt.gcf().set_size_inches(5,10)
plt.clf()
plt.plot( pot[:,0], pot[:,1], linewidth=2, label='$V(r)$', linestyle='--' )
for ii in range(Nstates):
  plt.plot( [Xmin,Xmax], [Ene[ii],Ene[ii]], linestyle='--' )
  minrho = np.min( rhoik[:,ii] )
  print minrho
  plt.plot( xgrid, rhoik[:,ii]-minrho+Ene[ii], label='$|\psi_{'+str(ii+1)+'}|$' )
plt.legend(loc='best')
plt.savefig('plotMathieu.pdf')
plt.savefig('plotMathieu.png',dpi=300)


