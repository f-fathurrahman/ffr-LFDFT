# Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
# arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
# Legendre n-point quadrature formula.
function gauleg( x1::Float64, x2::Float64, n::Int64)

    EPS = 3.0E-11 # Relative precision

    # int m, j, i;
    # double z1, z, xm, xl, pp, p3, p2, p1; // High precision is a good idea for this routine.

    # Output
    x = Array{Float64}(n)
    w = Array{Float64}(n)

    m = round( Int64, (n + 1)/2 )      # The roots are symmetric in the interval, so
    xm = 0.5*(x2 + x1) # we only have to find half of them.
    xl = 0.5*(x2 - x1)

    z1 = 100.0  # some initial number to get the while loop enter the first time
    pp = 0.0    # make pp visible outside for loop

    for i = 1:m # Loop over the desired roots.
        z = cos( pi*(i-0.25)/(n+0.5) )
        # Starting with the above approximation to the ith root, we enter the main loop of
        # refinement by Newton’s method.
        while abs(z-z1) > EPS
            p1 = 1.0;
            p2 = 0.0;
            for j=1:n # Loop up the recurrence relation to get the
                p3 = p2       # Legendre polynomial evaluated at z.
                p2 = p1
                p1 = ( (2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
            end
            # p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
            # by a standard relation involving also p2, the polynomial of one lower order.
            pp = n*(z*p1-p2)/(z*z-1.0)
            z1 = z
            z  = z1 - p1/pp # Newton’s method.
        end # while
        #
        x[i] = xm - xl*z;                  # Scale the root to the desired interval,
        x[n-i+1] = xm + xl*z;              # and put in its symmetric counterpart.
        w[i] = 2.0*xl/((1.0-z*z)*pp*pp) # Compute the weight
        w[n-i+1] = w[i]                 # and its symmetric counterpart.
    end
  return x, w
end

function test_gauleg()
    x1 = 0.0
    x2 = 2.0
    n = 11
    x, w = gauleg(x1, x2, n)
    println("Test gauleg")
    for i = 1:n
        @printf("%4d %18.10f %18.10f\n", i, x[i], w[i])
    end
end

test_gauleg()
