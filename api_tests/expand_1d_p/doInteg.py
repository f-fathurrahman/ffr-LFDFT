from sympy import *

x, L, a = symbols('x L a')

f = exp( a*cos(2*pi/L*(x-L/2)) )

df = diff( f, x )
d2f = diff( f, x, 2 )

pprint( simplify(df) )

pprint( simplify(d2f) )

