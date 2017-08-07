from __future__ import print_function
import scipy.special as sp
import numpy as np
from math import *

EPS = 2.220446049250313e-16

def gen_grid( h, N ):
    # generate grid points
    A = -(N-1)/2.0*h
    xgrid = []
    for i in range(N):
        xgrid.append(A + i*h)
    return xgrid

def eval_bfs( x, ibf, xgrid ):
    #
    xx = pi/h*(x-xgrid[ibf])
    if abs(xx) < EPS:
      xx = EPS
    return sin(xx)/xx / sqrt(h)

h = 0.2
N = 5
xgrid = gen_grid(h, N)
print(xgrid)
print(eval_bfs(0.2, 3, xgrid))

def compute_F( t, x_bar, h ):
    if x_bar < 1e-30:
        return sqrt(h)*sp.erf( pi/(2*h*t) )
    else:
        z = np.complex( pi/(2*h*t), t*x_bar )
        w_iz = sp.erfcx( z )
        f = exp( -t**2 * x_bar**2 )
        f = f - np.real( np.exp(-t**2 * x_bar**2 - z*z)*w_iz )
        f = f*sqrt(h)
        return f

t = 1.0
ibf1 = 4
ibf2 = 0
xx = xgrid[ibf2]
x_bar = fabs( xx - xgrid[ibf1] )
print('x_bar = %18.10f' % x_bar)
print( 'F = %18.10f' % compute_F( t, x_bar, h ))

#def eval_F(x):
#    global t, ibf1, ibf2
#    return exp(-t**2*( x - xgrid[ibf1] )**2) * eval_bfs(x, ibf1, xgrid)

# f = lambda x: exp(-t**2*( x - xgrid[ibf1] )**2) * eval_bfs(x, ibf1, xgrid)

"""
import matplotlib.pyplot as plt
plt.clf()
x = np.linspace(-5.0, 5.0, 200)
y = np.zeros(len(x))
for i in range(len(x)):
    y[i] = f(x[i])
plt.plot( x, y, label='f1' )
for i in range(len(x)):
    y[i] = eval_F(x[i])
plt.plot( x, y, label='f2' )
plt.legend()
plt.grid()
plt.savefig('integrand.png', dpi=300)
"""

# TFINAL = -10000

import mpmath as mp
# diagonal: ibf1 = ibf2, should gives erf
# non-diagonal: ibf1 != ibf2
def calc_F_v1( t, ibf1, ibf2 ):
    f = lambda x: exp(-t**2*( x - xgrid[ibf1] )**2) * eval_bfs(x, ibf2, xgrid)
    #print('mpmath v1 = %18.10f' % mp.quad( eval_F, [-mp.inf,mp.inf] ) )
    print('mpmath v2 = %18.10f' % mp.quad( f, [-mp.inf,mp.inf] ) )

# Lee-Tuckerman (2008)
def calc_F_v2( t, ibf1, ibf2 ):
    # integrand
    beta = pi/(h*t)
    xbar = xgrid[ibf1] - xgrid[ibf2]
    f = lambda x : exp(-x**2) * sin(beta*(x + t*xbar))/(x + t*xbar) * sqrt(h)/pi
    print('mpmath v3 = %18.10f' % mp.quad( f, [-mp.inf,mp.inf] ) )

calc_F_v1( t, ibf1, ibf2 )
#calc_F_v1( t, ibf1, ibf1 )

calc_F_v2( t, ibf1, ibf2 )
