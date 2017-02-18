from sympy import *

x, mu, sigma = symbols('x mu sigma')

f = 1/(sigma*sqrt(2*pi))*exp( -1/2*((x-mu)/sigma)**2 )

df = diff( f, x )
d2f = diff( f, x, 2 )

pprint( df )

pprint( d2f )

