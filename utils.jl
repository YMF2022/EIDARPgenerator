using Plots
using Random
using JSON
using HDF5
using Distributions
using DelimitedFiles

function generate_timetable(ts_lines::Vector{Transitline}, start_t::Float64, end_t::Float64, folder::String; shift = 5.0)
    # sym = collect('A':'Z')
    # tt = 5 # need to be changed as dist/v_mass
    max_duration = 0
    ts_end = 1
    for (l, line) in enumerate(ts_lines)
        tt = line.dist_ts /(line.v_ts / 60)
        header = [ts_end:ts_end+line.n_ts-1; "Direction"]
        n_dep = Int(2*(end_t - start_t)*60 ÷ line.freq)
        timetable = Array{Float64}(undef, n_dep, line.n_ts+1)
        last_0 = 20.0 + 3*(l-1) # direction 0: right to left/down
        last_1 = last_0 + shift # direction 1: left to right /up
        for dep in 1:n_dep
            if isinteger(dep/2)
                timetable[dep, line.n_ts+1] = 0
                for ts in line.n_ts:-1:1
                    timetable[dep,ts] = last_1 + abs(ts-line.n_ts)*(tt+line.dt)
                end
                last_1 += line.freq
            else
                timetable[dep, line.n_ts+1] = 1
                for ts in 1:line.n_ts
                    timetable[dep,ts] = last_0 + (ts-1)*(tt+line.dt)
                end
                last_0 += line.freq
            end
            open("$folder/timetable_line$l.csv", "w") do f
                writedlm(f, reshape(header, 1, :), ",")
                writedlm(f, timetable, ",")
            end

            if maximum(timetable) > max_duration
                max_duration = maximum(timetable)
            end
        end
        ts_end = ts_end + line.n_ts
    end
    return max_duration
end

function generate_customer(n_c, opr_len, opr_width, operation_time, detour_factor, v_bus, folder)
    cus_array = Array{Float64}(undef, n_c, 8)
    opr_width_half = opr_width/2 - 1
    opr_len_half = opr_len/2 - 1
    x_ori = rand(Uniform(-opr_len_half, opr_len_half), n_c)
    y_ori = rand(Uniform(-opr_width_half, opr_width_half), n_c)
    x_des = rand(Uniform(-opr_len_half, opr_len_half), n_c)
    y_des = rand(Uniform(-opr_width_half, opr_width_half), n_c)
    for c in 1:n_c
        x1, x2, y1, y2 = x_ori[c], x_des[c], y_ori[c], y_des[c]
        while sqrt((x1-x2)^2+(y1-y2)^2) < 2.0
            x1, x2, y1, y2 = rand(Uniform(-opr_len_half, opr_len_half), 4)
        end
        x_ori[c], x_des[c], y_ori[c], y_des[c] = x1, x2, y1, y2
    end
    cus_array[:,1] = x_ori
    cus_array[:,2] = y_ori
    cus_array[:,3] = x_des
    cus_array[:,4] = y_des

    # maximum travel time
    direct_dist = sqrt.((x_ori .- x_des).^2 + (y_ori .- y_des).^2)
    max_duration = 2 .*direct_dist./v_bus # ear_dep + 2 * direct travel time = late_arr
    direct_tt = direct_dist./v_bus
    max_tt = detour_factor.*direct_dist./v_bus
    
    # Time window
    dep_time = Vector{Float64}(undef, n_c)
    arr_time = Vector{Float64}(undef, n_c)
    for c in 1:n_c
        dep = operation_time * rand()
        arr = dep + max_duration[c]
        while arr > operation_time
            dep = operation_time * rand()
            arr = dep + max_duration[c]
        end
        dep_time[c] = dep
        arr_time[c] = arr
    end
    cus_array[:,5] = dep_time
    cus_array[:,6] = dep_time .+ 15
    cus_array[:,7] = direct_tt
    cus_array[:,8] = max_tt

    open("$folder/customers.csv", "w") do f
        writedlm(f, ["x_o" "y_o" "x_d" "y_d" "ear_dep_time" "late_dep_time" "direct_ridetime" "max_ridetime"], ",")
        writedlm(f, cus_array, ",")
    end
    return cus_array
end

function generate_trainstop(ts_lines::Vector{Transitline}, max_opr_radius::Float64, folder::String)
    n_ts = sum(get_ts_number.(ts_lines))
    ts_coords = Array{Any}(undef, n_ts, 4) # create coordinates Array
    opr_len, opr_width = 0.0, 0.0 # initilize operational area length and width
    last = 1
    for (i,line) in enumerate(ts_lines) 
        trans = line.ts_transfer
        dist_ts = line.dist_ts
        n_ts_line = line.n_ts
        if !iseven(i) # Horizontal line
            ts_coords[last:n_ts_line+last-1,1] .= [x for x in -((n_ts_line-1)*dist_ts/2):dist_ts:((n_ts_line-1)*dist_ts/2)]
            ts_coords[last:n_ts_line+last-1,2] .= zeros(n_ts_line)
            opr_width = max(opr_width, 2 * max_opr_radius + (n_ts_line-1)*dist_ts)
        else # Vertical line
            ts_coords[last:n_ts_line+last-1,1] .= zeros(n_ts_line)
            ts_coords[last:n_ts_line+last-1,2] .= [x for x in -((n_ts_line-1)*dist_ts/2):dist_ts:((n_ts_line-1)*dist_ts/2)]
            opr_len = max(opr_len, 2 * max_opr_radius + (n_ts_line-1)*dist_ts)
        end
        ts_coords[last:n_ts_line+last-1,3] .= i
        ts_coords[last:n_ts_line+last-1,4] = [i == trans ? 1 : 0 for i in 1:n_ts_line]

        last = last + n_ts_line
    end

    open("$folder/trainStops.csv", "w") do f
        writedlm(f, ["x" "y" "line" "transfer"], ",")
        writedlm(f, ts_coords, ",")
    end

    println(ts_coords)
    return ts_coords, opr_len, opr_width
end

function generate_charger(ts_coords::Matrix, cgr_speed::Float64, folder::String)
    cgr_info = zeros(1,3)
    cgr_info[:,3] .= cgr_speed
    open("$folder/chargers.csv", "w") do f
        writedlm(f, ["x" "y" "charging_speed"], ",")
        writedlm(f, cgr_info, ",")
    end
    return cgr_info[:,1:2]
end

function generate_bus(n_c::Int64, buses::Vector{Bustype}, folder_name::String)
    # busdata = Matrix()
    n_types = length(bus_types)
    open("$folder_name/buses.csv", "w") do f
        writedlm(f, ["capacity" "num" "speed" "consumption" "maxBattery"], ',')
        for (i,bus) in enumerate(buses)
            n_bus = Int(n_c/n_types)
            writedlm(f, [bus.capacity n_bus bus.v_bus bus.β bus.maxbattery], ',')
        end
    end
end

function depot_other(params::Parameters, max_duration::Float64, folder_name::String)
    # depot
    depots = hcat(params.depot...)
    depots = hcat(1:size(depots)[1], depots)
    open("$folder_name/depots.csv.csv", "w") do f
        writedlm(f, ["ID" "x" "y"], ',')
        writedlm(f, depots, ',')
    end

    # Output other parameter
    open("$folder_name/other_parameters.csv", "w") do f
        writedlm(f, ["service_time" "max_wlk_dist" "wlk_speed" "dwel_time" "dummy_charger" "start_time" "duration"], ',')
        writedlm(f, [params.μ params.maxwalkdist params.v_walk 1.0 params.charger_dummies 0 max_duration+15], ',')
    end
end

# Create a folder name
function foldername(upperfolder::String, n_line::Int64, n_c::Int64, replace::Int64)
    # check if upper folder exists
    if !isdir(upperfolder)
        mkdir(upperfolder)
    end
    folder_name = upperfolder * "l$(n_line)c$n_c"
    if !isdir(folder_name)
        mkdir(folder_name)
        i = 0
    else
        if replace == 0
            target_word = "l$(n_line)c$n_c"
            items = readdir(upperfolder)
            matchfolders = filter(item -> isdir(joinpath(upperfolder, item)) && contains(lowercase(item), lowercase(target_word)), items) 
            i = parse(Int, matchfolders[end][end]) + 1
            folder_name = folder_name * "-$i"
            mkdir(folder_name)
        else
            i = 0
            @warn "The generated instance replace the existing one."
        end
    end
    return folder_name, i
end

function graph(ts_coords, c_array, n_c, opr_width, folder; flag_annotate = 1)
    image = plot(title="Scenario",legendfontsize=7, legend=:false)
    ylims!(-opr_width/2-1, opr_width/2+1)
    xlims!(-opr_width/2-1, opr_width/2+1)
    # plot customers
    scatter!(c_array[:,1],c_array[:,2], label="origin", markershape=:circle, markercolor=:black, markersize=4)
    scatter!(c_array[:,3],c_array[:,4], label="destination", markershape=:utriangle, markercolor=:phase, markersize=5, markerstrokewidth=0)
    for c in 1:n_c
        annotate!(c_array[c,1]+0.2, c_array[c,2]+0.3, text("$c",10,:black))
        annotate!(c_array[c,3]+0.2, c_array[c,4]+0.3, text("$c",10,:phase))
    end
    # plot transit stations and line
    scatter!(ts_coords[:,1], ts_coords[:,2], label="Transit stops", 
            markershape=:star5, markercolor=:salmon, markersize=8, markerstrokewidth=0)
    n_ts_each = size(ts_coords)[1] ÷ 2
    plot!(ts_coords[1:n_ts_each,1], ts_coords[1:n_ts_each,2], color=:salmon, linewidth=4)
    plot!(ts_coords[n_ts_each+1:end,1], ts_coords[n_ts_each+1:end,2], color=:salmon, linewidth=4)
    sym = collect('A':'Z')
    if flag_annotate == 1
        for ts in 1:size(ts_coords)[1]
            if ts == n_ts_each + ceil(n_ts_each/2)
                annotate!(ts_coords[ts,1]+0.7, ts_coords[ts,2]-0.4, text("($(sym[ts]))",10,:salmon))
            else
                annotate!(ts_coords[ts,1]+0.3, ts_coords[ts,2]-0.4, text("$(sym[ts])",10,:salmon))
            end 
        end
    else
        for ts in 1:size(ts_coords)[1]
            if ts == n_ts_each + ceil(n_ts_each/2)
                annotate!(ts_coords[ts,1]+0.7, ts_coords[ts,2]-0.4, text("($ts)",10,:salmon))
            else
                annotate!(ts_coords[ts,1]+0.3, ts_coords[ts,2]-0.4, text("$ts",10,:salmon))
            end
        end
    end
    display(image)
    savefig(image, "$folder/fig.png")
end

function save_network(folder::String, ts_tt::Array, ts_arcs::Dict)
    # export network OD Matrix
    h5open("$folder/OD-3D.h5", "w") do file
        write(file, "OD", ts_tt)
    end

    #export OD pairs
    open("$folder/OD-pairs.json", "w") do f
        json_string = JSON.json(ts_arcs)
        JSON.print(f, json_string)
    end
end

# retrive data
# jdata = JSON.parsefile("cross/c10/OD-pairs.json")
# jd = JSON.parse(jdata)