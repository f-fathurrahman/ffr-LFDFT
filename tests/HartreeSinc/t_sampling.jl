function t_sampling( num_points1::Int64, num_points2::Int64,
                     t_i::Float64, t_l::Float64, t_f::Float64 )
    # New method by Sundholm (JCP 132, 024102, 2010).
    # Within two-divided regions ([t_i,t_l], [t_l,t_f]),
    # num_points1, num_points2 quadrature points are made, respectively.

    t_values = Array{Float64}(num_points1 + num_points2)
    w_t = Array{Float64}(num_points1 + num_points2)

    # Linear coord region:  [t_i, t_l]
    x_leg, w_leg = gauleg(t_i, t_l, num_points1);

    for j=1:num_points1
        t_values[j] = x_leg[j]
        w_t[j]      = w_leg[j]*2.0/sqrt(pi)
    end

    # Logarithmic coord region: [t_l, t_f]
    x_leg, w_leg = gauleg( log(t_l), log(t_f), num_points2)

    # Return the log-coord-partitioned points back to linear t-space.
    s_p = 0.0
    w_p = 0.0
    for j=1:num_points2
        s_p = x_leg[j]
        w_p = w_leg[j]
        x_leg[j] = exp(s_p)
        w_leg[j] = w_p * exp(s_p)
    end

    for j=1:num_points2
        t_values[num_points1+j] = x_leg[j]
        w_t[num_points1+j] = w_leg[j]*2.0/sqrt(pi)
    end

    return t_values, w_t

end
