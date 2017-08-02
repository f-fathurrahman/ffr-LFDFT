import SpecialFunctions

function compute_F( t::Float64, x_bar::Float64, scaling::Float64 )
    return_val = 0.0

    if x_bar < 1.e-30
        return_val= sqrt(scaling) * erf(pi/(2.0*scaling*t))
    else
        z = Complex128( pi/(2*scaling*t), t*x_bar )
        w_iz = erfcx(z)
        return_val = exp( -t*t*x_bar*x_bar )
        return_val = return_val - ( exp( -t*t*x_bar*x_bar - z*z ) * w_iz).re
        return_val = sqrt(scaling)*return_val
    end
    return return_val;
end

function test_compute_F()
     t = 0.1
     x_bar = 0.12
     scaling = 0.2

     F = compute_F(t, x_bar, scaling)

     @printf("F = %18.10f\n", F)
end

test_compute_F()

