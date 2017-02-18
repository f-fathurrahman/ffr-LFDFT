from sympy import *

x, Vm, L = symbols('x Vm L')

f = Vm + Vm*cos(2*pi*x/L)

df = diff( f, x )
d2f = diff( f, x, 2 )

pprint( df )

pprint( d2f )

