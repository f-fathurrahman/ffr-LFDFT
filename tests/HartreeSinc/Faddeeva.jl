function Faddeeva_erfcx( z::Complex128 )
    return ccall( (:_ZN8Faddeeva5erfcxESt7complexIdEd, "Faddeeva.so"),
        Complex128, (Complex128,), z)
end
