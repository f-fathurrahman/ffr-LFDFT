from sympy import *

x, y, z = symbols('x y z')

r = sqrt( x**2 + y**2 + z**2 )


rho = exp(-r**2)
phi = (2*pi)**(3/2) * erf( r/sqrt(2) )


term1 = diff( phi, x, 2 )
term2 = diff( phi, y, 2 )
term3 = diff( phi, z, 2 )

rho_ana = -4*pi*( term1 + term2 + term3 )

pprint( simplify(rho_ana) )

pprint( simplify(rho-rho_ana) )



