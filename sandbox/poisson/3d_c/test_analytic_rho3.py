from sympy import *

x, y, z = symbols('x y z')

fx = exp(-x**2)
fy = exp(-y**2)
fz = exp(-z**2)

phi = -4*pi*fx*fy*fz

#rho = diff(phi,x,2) + diff(phi,y,2) + diff(phi,z,2)
rho = diff(fx,x,2)*fy*fz + fx*diff(fy,y,2)*fz + fx*fy*diff(fz,z,2)

pprint( simplify(rho) )

pprint( simplify(phi) )


