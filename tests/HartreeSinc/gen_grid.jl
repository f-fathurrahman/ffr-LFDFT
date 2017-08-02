function gen_grid( size::Array{Int64}, scaling::Array{Float64} )
  Nmax = maximum(size)
  grid = zeros(3,Nmax)
  for i = 1:3
      for j = 1:size[i]
          grid[i,j] = -scaling[i]*( 0.5*(size[i]-1) - j + 1 )
      end
  end
  return grid
end
