function init_density( num_gaussian, position::Array{Float64}, coef::Array{Float64},
                       exponent::Array{Float64}, size, scaling, grid)

    weight = prod(scaling)

    total_size = prod(size)
    density = zeros(total_size)
    potential = zeros(total_size)

    for igauss = 1:num_gaussian
        exponent2 = exponent[igauss]^2
        for i_x = 1:size[1]
            for i_y = 1:size[2]
                for i_z = 1:size[3]
                    x = grid[1,i_x] - position[1,igauss]
                    y = grid[2,i_y] - position[2,igauss]
                    z = grid[3,i_z] - position[3,igauss]
                    r2 = x*x + y*y + z*z
                    # density
                    idx = (i_x-1)*size[2]*size[3] + (i_y-1)*size[3] + (i_z-1) + 1
                    density[idx]   += coef[igauss]*(exponent2/pi)^1.5 * exp(-exponent2*r2)*sqrt(weight)
                    # the potential (analytic) due to the Gaussian density
                    potential[idx] += coef[igauss]*erf(exponent[igauss]*sqrt(r2))/sqrt(r2)
                end
            end
        end
    end
    #
    return density, potential
end
