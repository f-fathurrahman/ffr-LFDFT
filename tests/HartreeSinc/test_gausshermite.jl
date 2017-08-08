using FastGaussQuadrature
using SpecialFunctions

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
function compute_F_v2( xgrid, ibf1, ibf2, t; NGaussHermite=50)
    x, w = gausshermite(NGaussHermite)
    h = xgrid[2] - xgrid[1]
    xbar = xgrid[ibf2] - xgrid[ibf1]
    beta = pi/(h*t)
    ss = 0.0
    for i = 1:NGaussHermite
        ss = ss + w[i] * integrand(x[i], xbar, beta, t)
    end
    return ss*sqrt(h)/pi
end

function test_main()
    h = 0.2
    N = 5
    xgrid = gen_grid(h, N)
    t = 1.0
    ibf1 = 2
    ibf2 = 1
    println( compute_F_v1(xgrid, ibf1, ibf2, t) )
    println( compute_F_v2(xgrid, ibf1, ibf2, t) )
end

function do_convergence_test( t; Nmax=500 )
    h = 0.2
    N = 5
    xgrid = gen_grid(h, N)
    ibf1 = 2
    ibf2 = 1
    F_choi = compute_F_v1(xgrid, ibf1, ibf2, t)
    F_old = 0.0
    @printf("t = %f\n", t)
    for n = 10:10:Nmax
        F = compute_F_v2(xgrid, ibf1, ibf2, t, NGaussHermite=n)
        Δ = abs(F-F_old)
        @printf("%8d %20.12f %20.12f %20.10e\n", n, F_choi, F, Δ)
        if Δ < 1e-12
            @printf("Already converge for n = %d\n", n)
            break
        end
        F_old = F
    end
end

import PyPlot
const plt = PyPlot

function do_plot( xgrid, ibf1, ibf2, t )
    h = xgrid[2] - xgrid[1]
    xbar = xgrid[ibf2] - xgrid[ibf1]
    beta = pi/(h*t)
    NptsPlot = 1000
    x = linspace(-1.0,1.0,1000)
    y = Array{Float64}(NptsPlot)
    for i = 1:NptsPlot
        y[i] = exp(-x[i]^2)*integrand( x[i], xbar, beta, t )
    end
    plt.clf()
    plt.plot(x, y, linewidth=2)
    plt.savefig("Integrand1.png", dpi=300)
end

function test_plot()
    h = 0.2
    N = 5
    xgrid = gen_grid(h, N)
    t = 1e-2
    ibf1 = 2
    ibf2 = 1
    do_plot( xgrid, ibf1, ibf2, t )
end

test_plot()

#test_main()

# Generally fast convergence for large t
#do_convergence_test(10000.0)
#do_convergence_test(1000.0)
#do_convergence_test(100.0)
#do_convergence_test(10.0)
#do_convergence_test(1.0)

#do_convergence_test(0.9)
#do_convergence_test(0.8)
#do_convergence_test(0.7)
#do_convergence_test(0.6)

#do_convergence_test(0.5)
#do_convergence_test(0.4)
#do_convergence_test(0.3)
#do_convergence_test(0.2)

#do_convergence_test(1e-1,Nmax=2000)
#do_convergence_test(1e-2,Nmax=2000)
#do_convergence_test(1e-3,Nmax=2000)
#do_convergence_test(1e-4,Nmax=2000)
#do_convergence_test(1e-5,Nmax=2000)
