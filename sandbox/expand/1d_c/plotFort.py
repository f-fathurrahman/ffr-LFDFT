#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from math import pi

# HARDCODED
XLIM = [-5.0,5.0]
YLIM = [-0.2,1.4]

def do_Npoints(Npts):
  LF_coefs = np.loadtxt('COARSE_N_' + str(Npts))
  
  dat1 = np.loadtxt('DENSE_N_' + str(Npts))
  xx = dat1[:,0]
  yy = dat1[:,1]
  g  = dat1[:,2]
  
  plt.clf()
  plt.plot( LF_coefs[:,0], LF_coefs[:,1], label='lbf-coef', marker='o', linestyle='_' )
  plt.plot( xx, yy, label='lbf-eval' )
  plt.plot( xx, g, label='ana' )
  plt.legend()
  plt.grid()
  plt.xlim(XLIM)
  plt.ylim(YLIM)
  plt.title('N = ' + str(Npts))
  plt.savefig('plot_' + str(Npts) + '.pdf')
  plt.savefig('plot_' + str(Npts) + '.png', dpi=300)


Npts = [10, 15, 20, 25, 30]  # HARDCODED

for n in Npts:
  do_Npoints(n)

