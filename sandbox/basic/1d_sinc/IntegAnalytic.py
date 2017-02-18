from sympy import *

x = Symbol('x')

xi = -3/2
h = 1

Li = 1/sqrt(h) * sin( pi*(x-xi)/h ) / (pi*(x-xi)/h)

pprint( simplify(Li) )

nn = integrate( Li*Li, (x,-oo,oo) )
print nn
