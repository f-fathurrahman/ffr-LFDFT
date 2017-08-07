from __future__ import print_function
from sympy import *

# numerical value of PI
PI = pi.n()

def gen_basis_set( h, N ):
    # generate grid points
    A = -(N-1)/2.0*h
    xgrid = []
    for i in range(N):
        xgrid.append(A + i*h)
    #
    x = symbols('x')
    bfs = []
    for i in range(N):
        bfs.append( sin(PI/h*(x-xgrid[i]))/sqrt(h)/(PI*(x-xgrid[i])/h) )
    #
    return xgrid, bfs

h = 0.2
N = 5
x = symbols('x')
xgrid, bfs = gen_basis_set(h, N)


# need to use scipy
import scipy.special as sp
import numpy as np
import math

def compute_F( t, x_bar, h ):
    if x_bar < 1e-30:
        return math.sqrt(h)*sp.erf( np.pi/(2*h*t) )
    else:
        z = np.complex( np.pi/(2*h*t), t*x_bar )
        w_iz = sp.erfcx( z )
        f = math.exp( -t**2 * x_bar**2 )
        f = f - np.real( np.exp(-t**2 * x_bar**2 - z*z)*w_iz )
        f = f*math.sqrt(h)
        return f

t = 10.0
xx = 0.1
ibf1 = 0
x_bar = math.fabs( xx - xgrid[ibf1] )
print('x_bar = %18.10f' % x_bar)
print( 'F = %18.10f' % compute_F( t, x_bar, h ))

"""
t = 0.1
f = exp(-t**2*(x-x_alpha)**2) * u
print(f)

print( Integral( f, (x,-oo,oo) ).evalf() )
"""
