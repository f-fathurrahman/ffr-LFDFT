# Name for icpc:
# 0000000000002810 t _Z10erfcx_y100d
# 00000000000026e0 T _ZN8Faddeeva5erfcxEd
# 0000000000000d30 T _ZN8Faddeeva5erfcxESt7complexIdEd

import SpecialFunctions

z = 1.0 + im*1.0
println("z = ", z)

println( ccall( (:_ZN8Faddeeva5erfcxESt7complexIdEd, "Faddeeva.so"),
    Complex128, (Complex128,), z) )

println( erfcx(z) )

#println( ccall( (:_ZN8Faddeeva5erfcxEd, "Faddeeva.so"),#
#    Complex128, (Complex128,), 1.0 + im*1.0) )
