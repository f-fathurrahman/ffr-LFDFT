from sympy import *

#r = sqrt( x**2 + y**2 + z**2 )
r = symbols('r')

rho = exp(-r**2/2)
phi = (2*pi)**(3/2) * erf( r/sqrt(2) ) / r

# Laplacian operator (radial)
term1 = r**2 * diff( phi, r, 1 )
rho_ana = (1/r**2) * diff( term1, r, 1 )

rho_ana = -rho_ana/(4*pi)

pprint( simplify(rho_ana) )

pprint( simplify(rho-rho_ana) )



