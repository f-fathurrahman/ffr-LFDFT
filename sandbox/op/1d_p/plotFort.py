#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

dat11 = np.loadtxt('fort.11')

xx = np.loadtxt('fort.12')
yy = np.loadtxt('fort.13')
g  = np.loadtxt('fort.14')
yy1 = np.loadtxt('fort.15')
dg = np.loadtxt('fort.16')
yy2 = np.loadtxt('fort.17')
d2g = np.loadtxt('fort.18')

plt.clf()
plt.plot( dat11[:,0], dat11[:,1], label='lbf-coef', marker='o', linestyle='_' )
plt.plot( xx, yy, label='lbf-eval' )
plt.plot( xx, g, label='ana' )
plt.legend()
plt.grid()
plt.savefig('plot0.pdf')

plt.clf()
plt.plot( xx, yy1, label='lbf-d1-eval', marker='o', linestyle='_')
plt.plot( xx, dg, label='ana-d1' )
plt.legend()
plt.grid()
plt.savefig('plot1.pdf')

plt.clf()
plt.plot( xx, yy2, label='lbf-d2-eval', marker='o', linestyle='_' )
plt.plot( xx, d2g, label='ana-d2' )
plt.legend()
plt.grid()
plt.savefig('plot2.pdf')

