from sympy import *

x, L = symbols('x L')

f = exp( cos( 2*pi/L * (x-L/2) ) )

df = diff( f, x )
d2f = diff( f, x, 2 )

#pprint( df )

pprint( simplify(f) )
pprint( simplify(d2f) )

