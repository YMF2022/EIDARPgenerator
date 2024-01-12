function generate_customer(n_c::Int64, opr_len, opr_width, operation_time, detour_factor, v_bus, tw, ts_coords, folder::String, location::Function)
    cus_array = Array{Float64}(undef, n_c, 7)
    opr_width_half = opr_width/2 - 1
    opr_len_half = opr_len/2 - 1
    cus_array[:,1:4], direct_dist = location(n_c, opr_len_half, opr_width_half, ts_coords)

    # maximum travel time
    direct_tt = direct_dist./v_bus
    cus_array[:,7] = direct_tt
    max_tt = detour_factor.*direct_dist./v_bus
    
    # Time window: need to have dep_ear, dep_late
    depot = params.depot[1]
    @show depot
    for c in 1:n_c
        xo, yo, xd, yd = cus_array[c,1:4]
        time_to_depot = sqrt((xo-depot[1])^2+(yo-depot[2])^2)/v_bus
        dep_ear = operation_time * rand()
        dep_late = dep_ear + tw
        arr_late = dep_late + max_tt[c]
        while (dep_ear + tw + max_tt[c] > operation_time) && (dep_late < time_to_depot)
            dep_ear = operation_time * rand()
            dep_late = dep_ear + tw
            arr_late = dep_late + max_tt[c]
        end
        arr_ear = dep_ear + direct_tt[c]
        cus_array[c,5:6] .= dep_ear, dep_late
    end

    open("$folder/customers.csv", "w") do f
        writedlm(f, ["x_o" "y_o" "x_d" "y_d" "ear_dep_time" "late_dep_time" "direct_ridetime"], ",")
        writedlm(f, cus_array, ",")
    end
    return cus_array
end

function random_spread(n_c, opr_len_half, opr_width_half, ts_coords)
    x_ori = rand(Uniform(-opr_len_half, opr_len_half), n_c)
    y_ori = rand(Uniform(-opr_width_half, opr_width_half), n_c)
    x_des = rand(Uniform(-opr_len_half, opr_len_half), n_c)
    y_des = rand(Uniform(-opr_width_half, opr_width_half), n_c)
    direct_dist = []
    for c in 1:n_c
        xo, xd, yo, yd = x_ori[c], x_des[c], y_ori[c], y_des[c]
        dist = sqrt((xo-xd)^2+(yo-yd)^2)
        while dist < 2.0
            xo, xd, yo, yd = rand(Uniform(-opr_len_half, opr_len_half), 4)
            dist = sqrt((xo-xd)^2+(yo-yd)^2)
        end
        x_ori[c], x_des[c], y_ori[c], y_des[c] = xo, xd, yo, yd
        push!(direct_dist, dist)
    end
    return hcat(x_ori, y_ori, x_des, y_des), direct_dist
end

function closeto_ts(n_c, opr_len_half, opr_width_half, ts_coords)
    r_max = 0.5
    n_ts = size(ts_coords)[1]
    cus_coords = Matrix{Float64}(undef, (n_c, 4))
    direct_dist = Vector{Float64}(undef, n_c)
    for c in 1:n_c
        θ = rand(0:360)
        r = r_max*rand()
        ts_o = rand(1:n_ts)
        ts_d = rand([i for i in 1:n_ts if i!=ts_o])
        xo = cosd(θ)*r + ts_coords[ts_o,1]
        yo = sind(θ)*r + ts_coords[ts_o,2]
        xd = cosd(θ)*r + ts_coords[ts_d,1]
        yd = sind(θ)*r + ts_coords[ts_d,2]
        dist = sqrt((xo-xd)^2+(yo-yd)^2)
        direct_dist[c] = dist
        cus_coords[c,:] .= xo, yo, xd, yd
    end
    return cus_coords, direct_dist
end

function farfrom_ts(n_c, opr_len_half, opr_width_half, ts_coords)
    cus_coords = Matrix{Float64}(undef, (n_c, 4))
    direct_dist = Vector{Float64}(undef, n_c)
    r_max = 3
    r_min = 1.5
    ts_list = [1, 3, 4, 6]
    for c in 1:n_c
        θ = rand(0:360)
        r = rand(Uniform(r_min,r_max))
        ts_o = rand(ts_list)
        ts_d = rand([i for i in ts_list if i!=ts_o])
        xo = cosd(θ)*r + ts_coords[ts_o,1]
        yo = sind(θ)*r + ts_coords[ts_o,2]
        xd = cosd(θ)*r + ts_coords[ts_d,1]
        yd = sind(θ)*r + ts_coords[ts_d,2]
        dist = sqrt((xo-xd)^2+(yo-yd)^2)
        while dist < 2.0
            θ = rand(0:360)
            r = rand(Uniform(r_min,r_max))
            ts_o = rand(ts_list)
            ts_d = rand([i for i in ts_list if i!=ts_o])
            xo = cosd(θ)*r + ts_coords[ts_o,1]
            yo = sind(θ)*r + ts_coords[ts_o,2]
            xd = cosd(θ)*r + ts_coords[ts_d,1]
            yd = sind(θ)*r + ts_coords[ts_d,2]
            dist = sqrt((xo-xd)^2+(yo-yd)^2)
        end
        cus_coords[c,:] .= xo, yo, xd, yd
        direct_dist[c] = dist
    end
    return cus_coords, direct_dist
end

function far_close_ts(n_c, opr_len_half, opr_width_half, ts_coords)   
    n_c_1 =  div(n_c,2)
    n_c_2 = n_c - n_c_1
    cus_coords1, dist1 = closeto_ts(n_c_1, opr_len_half, opr_width_half, ts_coords)
    cus_coords2, dist2 =farfrom_ts(n_c_2, opr_len_half, opr_width_half, ts_coords)
    cus_coords = vcat(cus_coords1, cus_coords2)
    direct_dist = vcat(dist1, dist2)
    return cus_coords, direct_dist
end