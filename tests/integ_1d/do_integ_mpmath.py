import mpmath as mp
from math import pi, exp, cos

def do_integ( scal=0.5 ):
    L = 10.0
    c1 = L*scal
    f = lambda x: exp( 1.5*cos(2.0*pi/L *(x - c1) ) )
    print('mpmath = %18.10f' % mp.quad( f, [0,L] ) )

do_integ( scal=0.1 )
do_integ( scal=0.5 )
do_integ( scal=0.9 )

