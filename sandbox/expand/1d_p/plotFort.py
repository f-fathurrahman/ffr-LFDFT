#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from math import pi

LF_coefs = np.loadtxt('fort.11')

xx = np.loadtxt('fort.12')
yy = np.loadtxt('fort.13')
g  = np.loadtxt('fort.14')

plt.clf()
plt.plot( LF_coefs[:,0], LF_coefs[:,1], label='lbf-coef', marker='o', linestyle='_' )
plt.plot( xx, yy, label='lbf-eval' )
plt.plot( xx, g, label='ana' )
plt.legend()
plt.grid()
plt.savefig('plot0.pdf')


