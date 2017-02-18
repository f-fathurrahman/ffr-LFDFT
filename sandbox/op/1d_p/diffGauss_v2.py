from sympy import *

x, y, z, mu, sigma = symbols('x y z mu sigma')

r = sqrt( x**2 + y**2 + z**2 )
f = 1/(sqrt((2*pi*sigma**2))**3) * exp( -r**2/(2*sigma**2) )
# beware of the expression such as x**(3/2)
# probably I should use x**( Integer(3)/Integer(2) )


df = diff( f, x )
d2f = diff( f, x, 2 )

pprint( f )
print ""
pprint( df )
print ""
pprint( d2f )

