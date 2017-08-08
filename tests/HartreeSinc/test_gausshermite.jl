using FastGaussQuadrature
using SpecialFunctions

const NGaussHermite = 300
const xw_gauss = gausshermite(NGaussHermite)

function gen_grid( h, N )
    # generate grid points
    A = -(N-1)/2.0*h
    xgrid = Array{Float64}(N)
    for i = 1:N
        xgrid[i] = A + (i-1)*h
    end
    return xgrid
end


function eval_bfs( x, ibf, xgrid )
    h = xgrid[2] - xgrid[1]
    xx = pi/h*(x-xgrid[ibf])
    if abs(xx) < eps()
        xx = eps()
    end
    return sin(xx) / xx / sqrt(h)
end

# Choi2016
function compute_F_v1( xgrid, ibf1, ibf2, t )
    h = xgrid[2] - xgrid[1]
    x_bar = abs( xgrid[ibf2] - xgrid[ibf1] )
    if x_bar < 1e-30
        return sqrt(h)*erf( pi/(2*h*t) )
    else
        z = pi/(2*h*t) + im*t*x_bar
        w_iz = erfcx( z )
        f = exp( -t^2 * x_bar^2 )
        f = f - real( exp(-t^2 * x_bar^2 - z*z)*w_iz )
        f = f*sqrt(h)
        return f
    end
end

function integrand( x, xbar, beta, t )
    return sin( beta*(x + t*xbar) ) / ( x + t*xbar )
end

# LeeTuckerman2008
function compute_F_v2( xgrid, ibf1, ibf2, t )
    h = xgrid[2] - xgrid[1]
    xbar = xgrid[ibf2] - xgrid[ibf1]
    beta = pi/(h*t)
    ss = 0.0
    for i = 1:NGaussHermite
        ss = ss + xw_gauss[2][i] * integrand(xw_gauss[1][i], xbar, beta, t)
    end
    return ss*sqrt(h)/pi
end

function test_main()
    h = 0.2
    N = 5
    xgrid = gen_grid(h, N)
    t = 0.1
    ibf1 = 2
    ibf2 = 1
    println( compute_F_v1(xgrid, ibf1, ibf2, t) )
    println( compute_F_v2(xgrid, ibf1, ibf2, t) )
end

test_main()
