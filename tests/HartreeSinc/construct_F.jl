include("compute_F.jl")

function construct_F(axis::Int64, t_size::Int64, t_values::Array{Float64},
        size, scaling, grid )
    F_values = Array{Float64,2}( t_size, size[axis]^2 )

    for i_t = 1:t_size
        for i = 1:size[axis]
            for j = 1:size[axis]
                idx = (i-1)*size[axis] + (j-1) + 1
                F_values[i_t,idx] = compute_F( t_values[i_t],
                                               abs( grid[axis,i]-grid[axis,j] ),
                                               scaling[axis] )
            end
        end
    end
    return F_values
end


function construct_F_debug(axis::Int64, t_size::Int64, t_values::Array{Float64},
        size, scaling, grid )
    F_values = Array{Float64,2}( t_size, size[axis]^2 )

    for i_t = 1:t_size
        for i = 1:size[axis]
            for j = 1:size[axis]
                idx = (i-1)*size[axis] + (j-1) + 1
                F_values[i_t,idx] = compute_F( t_values[i_t],
                                               abs( grid[axis,i]-grid[axis,j] ),
                                               scaling[axis] )
                @printf("%18.10f %18.10f %18.10f\n", t_values[i_t],
                                               abs( grid[axis,i]-grid[axis,j] ),
                                               scaling[axis])
                @printf("%8d %8d: %18.10f\n", i_t, idx, F_values[i_t,idx])
            end
        end
    end
    return F_values
end
