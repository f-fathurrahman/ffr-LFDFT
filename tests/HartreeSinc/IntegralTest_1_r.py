from __future__ import print_function
import mpmath as mp
from math import *

def func( t, r ):
    return exp( -r**2 * t **2)

r = 10.0

f = lambda t : func(t, r)

res1 = mp.quad( f, (0,mp.inf) ) * 2.0 / sqrt(pi)
print('res1 = %f' % res1)
