from __future__ import print_function
from sympy import *

h = 0.2
x = symbols('x')

PI = pi.n()
# basis function
x_alpha = 0.2
u = sin(PI/h*(x-x_alpha))/sqrt(h)/(PI*(x-x_alpha)/h)
print(u)

t = 0.1
f = exp(-t**2*(x-x_alpha)**2) * u
print(f)

print( Integral( f, (x,-oo,oo) ).evalf() )
