from sympy import *

x, y, z = symbols('x y z')
Gx, Gy, Gz = symbols('Gx Gy Gz')

Gr = Gx*x + Gy*y + Gz*z

PW = exp(I*Gr)

dx_PW = diff( PW, x )
d2x_PW = diff( PW, x, 2 )

nabla_PW = diff( PW, x, 2 ) + diff( PW, y, 2 ) + diff( PW, z, 2 )

pprint( dx_PW)

pprint( d2x_PW )

pprint( nabla_PW )

