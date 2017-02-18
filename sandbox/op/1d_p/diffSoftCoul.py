from sympy import *

x, Q, A, L = symbols('x Q A L')

f = -Q**2 / sqrt(A**2 + (x-L/2)**2)

df = diff( f, x )
d2f = diff( f, x, 2 )

pprint( df )

pprint( d2f )

