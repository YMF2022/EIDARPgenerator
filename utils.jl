using Plots
using Random
using Distributions
using DelimitedFiles
using Distances
using DataFrames
using CSV

function generate_timetable(ts_lines::Vector{T}, start_t::Float64, end_t::Float64, folder::String; shift = 5.0) where T<:Transitline
    max_duration = 0
    ts_end = 1
    for (l, line) in enumerate(ts_lines)
        tt = line.dist_ts /(line.v_ts / 60)
        header = [ts_end:ts_end+line.n_ts-1; "Direction"]
        n_dep = Int(2*(end_t - start_t)*60 ÷ line.freq)
        timetable = Matrix{Float64}(undef, n_dep, line.n_ts+1)
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

function generate_trainstop(ts_lines::Vector{Crossline}, max_opr_radius::Float64, folder::String)
    n_ts = sum(getfield.(ts_lines, :n_ts))
    ts_stops = Array{Any}(undef, n_ts, 4) # create coordinates Array
    opr_len, opr_width = 0.0, 0.0 # initilize operational area length and width
    last = 1
    for (i,line) in enumerate(ts_lines) 
        trans = line.ts_transfer
        dist_ts = line.dist_ts
        n_ts_line = line.n_ts
        if !iseven(i) # Horizontal line
            ts_stops[last:n_ts_line+last-1,1] .= [x for x in -((n_ts_line-1)*dist_ts/2):dist_ts:((n_ts_line-1)*dist_ts/2)]
            ts_stops[last:n_ts_line+last-1,2] .= zeros(n_ts_line)
            opr_width = max(opr_width, 2 * max_opr_radius + (n_ts_line-1)*dist_ts)
        else # Vertical line
            ts_stops[last:n_ts_line+last-1,1] .= zeros(n_ts_line)
            ts_stops[last:n_ts_line+last-1,2] .= [x for x in -((n_ts_line-1)*dist_ts/2):dist_ts:((n_ts_line-1)*dist_ts/2)]
            opr_len = max(opr_len, 2 * max_opr_radius + (n_ts_line-1)*dist_ts)
        end
        ts_stops[last:n_ts_line+last-1,3] .= i
        ts_stops[last:n_ts_line+last-1,4] = [i == trans ? 1 : 0 for i in 1:n_ts_line]

        last = last + n_ts_line
    end

    open("$folder/trainStops.csv", "w") do f
        writedlm(f, ["x" "y" "line" "transfer"], ",")
        writedlm(f, ts_stops, ",")
    end

    ts_stops = DataFrame(ts_stops, [:x, :y, :line, :transfer])

    @info("Operational area: $opr_len*$opr_width km")
    return ts_stops, opr_len, opr_width
end

function read_transit_network(networkshape::Symbol, folder::String, params::Parameters)
    ts_stops = CSV.read(folder * string(networkshape) * "-network.csv", DataFrame; header = true)
    ts_coords = Matrix(ts_stops[:,[:x,:y]])
    ts_dist_matrix = pairwise(Euclidean(), eachrow(ts_coords), eachrow(ts_coords))
    opr_len = maximum(ts_stops.x) - minimum(ts_stops.x) + 2*params.max_opr_radius
    opr_width = maximum(ts_stops.y) - minimum(ts_stops.y) + 2*params.max_opr_radius
    return ts_stops, opr_len, opr_width
end

function generate_charger(ts_stops::DataFrame, cgr_speed::Float64, folder::String)
    cgr_info = zeros(1,3)
    cgr_info[:,3] .= cgr_speed
    open("$folder/chargers.csv", "w") do f
        writedlm(f, ["x" "y" "charging_speed"], ",")
        writedlm(f, cgr_info, ",")
    end
    return cgr_info[:,1:2]
end

function generate_bus(n_c::Int64, buses::Vector{Bustype}, depots::Vector, folder_name::String)
    n_depots = length(depots)
    n_types = length(buses)
    open("$folder_name/buses.csv", "w") do f
        writedlm(f, ["ID" "type" "capacity" "speed" "consumption" "maxBattery" "depot"], ',')
        last = 0
        for type in 1:n_types
            bus = buses[type]
            n_bus = Int(ceil(n_c/n_types))
            for i in 1:n_bus
                depot = rand(1:n_depots)
                writedlm(f, Any[i+last type bus.capacity bus.v_bus bus.β bus.maxbattery depot], ',')
            end
            last += n_bus
        end
    end
end

function depot_other(params::Parameters, max_duration::Float64, folder_name::String)
    # depot
    depots = hcat(params.depot...)
    depots = hcat(1:size(depots)[1], depots)
    open("$folder_name/depots.csv", "w") do f
        writedlm(f, ["x" "y"], ',')
        writedlm(f, [depots[:,2] depots[:,3]], ',')
    end

    # Output other parameter
    open("$folder_name/other_parameters.csv", "w") do f
        writedlm(f, ["service_time" "max_wlk_dist" "wlk_speed" "dwel_time" "dummy_charger" "detour_factor" "max_wait_time" "start_time" "duration"], ',')
        writedlm(f, Any[params.μ params.max_walkdist params.v_walk 1.0 params.charger_dummies params.detour_factor params.max_waittime 0.0 max_duration+15], ',')
    end
end

# Create a folder name
function foldername(upperfolder::String, n_line::Int64, n_c::Int64, n_depot::Int64, n_bt::Int64, replace::Int64)
    # check if upper folder exists
    if !isdir(upperfolder)
        mkdir(upperfolder)
    end
    folder_name = upperfolder * "l$(n_line)-c$n_c-d$n_depot-bt$n_bt"
    if !isdir(folder_name)
        mkdir(folder_name)
        i = 0
    else
        if replace == 0
            target_word = "l$(n_line)-c$n_c-d$n_depot-bt$n_bt"
            items = readdir(upperfolder)
            matchfolders = filter(item -> isdir(joinpath(upperfolder, item)) && contains(lowercase(item), lowercase(target_word)), items) 
            i = length(matchfolders) + 1
            folder_name = folder_name * "_$i"
            mkdir(folder_name)
        else
            i = 0
            @warn "The generated instance replace an existing one. Set replace = 0 if you don't want to replace. "
        end
    end
    return folder_name, i
end

function graph(ts_stops, ts_lines, c_array, n_c, opr_len, opr_width, cgr_coords, folder; flag_annotate = 1)
    image = plot(title="Scenario",legendfontsize=7, legend=:true)
    ylims!(-opr_width/2-1, opr_width/2+1)
    xlims!(-opr_len/2-1, opr_len/2+1)
    # plot customers
    scatter!(c_array[:,1],c_array[:,2], label="origin", markershape=:circle, markercolor=:black, markersize=4)
    scatter!(c_array[:,3],c_array[:,4], label="destination", markershape=:utriangle, markercolor=:phase, markersize=5, markerstrokewidth=0)
    for c in 1:n_c
        annotate!(c_array[c,1]+0.2, c_array[c,2]+0.3, text("$c",10,:black))
        annotate!(c_array[c,3]+0.2, c_array[c,4]+0.3, text("$c",10,:phase))
    end

    # plot charging stations
    scatter!(cgr_coords[:,1], cgr_coords[:,2], label="Charger", 
    markershape=:utriangle, markercolor=:lightblue, markersize=4, markerstrokewidth=0)

    # plot transit stations and lines
    graph_ts(ts_stops, ts_lines; flag_annotate = flag_annotate)

    # display(image)
    savefig(image, "$folder/fig.png")
end

function graph_ts(ts_stops, ts_lines; flag_annotate = 1)
    image = plot!(title="Scenario",legendfontsize=7, legend=:true)

    # plot transit stations and lines
    colors = [:darkolivegreen, :navy, :firebrick4]
    n_line = length(ts_lines)
    for i in 1:n_line
        linecoords_x = ts_stops.x[ts_stops.line .== i]
        linecoords_y = ts_stops.y[ts_stops.line .== i]
        plot!(linecoords_x, linecoords_y, color=colors[i], linewidth=4, label=false)
        # plot!(linecoords_x, linecoords_y, color=:gray, linewidth=4, label=false)
        if ts_lines[i].shape == :circle 
            plot!([linecoords_x[1],linecoords_x[end]], [linecoords_y[1],linecoords_y[end]], color=colors[i], linewidth=4, label=false)
        end
    end
    scatter!(ts_stops.x, ts_stops.y, label="Transit stops", markershape=:star5, markercolor=:gray27, markersize=8, markerstrokewidth=0)
    sym = collect('A':'Z')

    if flag_annotate == 1 
        # name transit stops by letters
        plotted_coords = Vector{Tuple{Float64, Float64}}()
        for ts in 1:size(ts_stops)[1]
            c = colors[ts_stops.line[ts]]
            if (ts_stops.x[ts],ts_stops.y[ts]) ∈ plotted_coords
                annotate!(ts_stops.x[ts]+0.8, ts_stops.y[ts]-0.4, text("($(sym[ts]))",10, c)) 
            else
                annotate!(ts_stops.x[ts]+0.3, ts_stops.y[ts]-0.4, text("$(sym[ts])",10, c))
                push!(plotted_coords, (ts_stops.x[ts],ts_stops.y[ts]))
            end 
        end
    else 
        # name transit stops by numbers
        plotted_coords = Vector{Tuple{Float64, Float64}}()
        for ts in 1:size(ts_stops)[1]
            if (ts_stops.x[ts],ts_stops.y[ts]) ∈ plotted_coords
                annotate!(ts_stops.x[ts]+0.8, ts_stops.y[ts], text("($ts)",10,colors[ts_stops[ts,3]]))
            else
                annotate!(ts_stops.x[ts]+0.3, ts_stops.y[ts]-0.4, text("$ts",10,colors[ts_stops[ts,3]]))
                push!(plotted_coords, (ts_stops.x[ts],ts_stops.y[ts]))
            end
        end
    end

    display(image)
end