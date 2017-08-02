using SpecialFunctions

include("gauleg.jl")
include("t_sampling.jl")
include("gen_grid.jl")
include("init_density.jl")
include("construct_F.jl")

function main()

    # Input
    scaling = [0.2, 0.2, 0.2]
    size = [25,25,25]
    num_points1 = 5
    num_points2 = 10

    t_i = 0.0
    t_l = 2.0
    t_f = 10000.0

    num_gaussian = 2
    #
    positions = Array{Float64}(3,num_gaussian)
    coefs = Array{Float64}(num_gaussian)
    exponents = Array{Float64}(num_gaussian)

    # First Gaussian
    ipos = 1
    positions[1,ipos] = 0.1
    positions[2,ipos] = 0.0
    positions[3,ipos] = 0.0
    coefs[ipos]     = 1.0
    exponents[ipos] = 3.0

    # Second Gaussian
    ipos = 2
    positions[1,ipos] = 1.5
    positions[2,ipos] = 0.0
    positions[3,ipos] = 0.0
    coefs[ipos]     = 1.0
    exponents[ipos] = 4.0

    t_size = num_points1 + num_points2;
    total_size = prod(size)

    @printf("Calculated by Julia\n")
    @printf("======== input information ========\n")
    @printf("size: %d %d %d\n", size[1], size[2], size[3])
    @printf("total_size: %d\n", total_size)
    @printf("scaling: %f %f %f\n", scaling[1], scaling[2], scaling[3])
    @printf("num_points1: %d\n", num_points1)
    @printf("num_points2: %d\n", num_points2)
    @printf("num_gaussian: %d\n", num_gaussian)
    @printf("===================================\n")

    t_values, wt = t_sampling( num_points1, num_points2, t_i, t_l, t_f)

    for i = 1:t_size
        @printf("%18.10f %18.10f\n", t_values[i], wt[i])
    end

    grid = gen_grid(size, scaling)

    @printf("Grid points along x:\n")
    for j = 1:size[1]
        @printf("%18.10f\n", grid[1,j])
    end

    density, anal_potential =
        init_density( num_gaussian, positions, coefs, exponents, size, scaling, grid )

    density_norm = 0.0
    anal_energy = 0.0
    for j = 1:total_size
        density_norm += density[j]
        anal_energy  += 0.5*anal_potential[j]*density[j]
    end
    density_norm *= sqrt(prod(scaling))
    anal_energy  *= sqrt(prod(scaling))
    @printf("norm of density: %18.10f\n", density_norm)
    @printf("Analytic energy: %18.10f\n", anal_energy)

    return

    F_xs = construct_F_debug(1, t_size, t_values, size, scaling, grid)
    F_ys = construct_F(2, t_size, t_values, size, scaling, grid)
    F_zs = construct_F(3, t_size, t_values, size, scaling, grid)

    println("sum(F_xs) = ", sum(F_xs))
    println("sum(F_xs) = ", sum(F_ys))
    println("sum(F_xs) = ", sum(F_zs))

    println("Program ended normally")
end

main()
