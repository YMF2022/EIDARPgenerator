include("EIDARPparams.jl")
include("customer.jl")
include("utils.jl")


"""
generate(n_c::Int64, params::Parameters; upperfolder = "cross/", replace = 1)

Generate synthetic data for a transportation system simulation.

# Arguments
- `n_c::Int64`: Number of customers.
- `params::Parameters`: A structure containing parameters for EIDARP.
- `upperfolder::String`: The parent folder for storing generated data. Default is "cross/".
- `replace::Int`: A flag indicating whether to replace existing instance (1) or not (0). Default is 1.

# Returns
- The function generates and saves synthetic data for a transportation system simulation. It doesn't return any value.

# Examples
```julia
julia> generate(10, my_params)
```
"""
function generate(n_c::Int64, ts_lines::Vector{Crossline}, params::Parameters, networkshape::Symbol; replace = 1, location = random_spread)
    upperfolder = string(networkshape) * "/"
    folder_name, seed = foldername(upperfolder, length(ts_lines), n_c, length(params.depot), length(params.buses), replace)
    Random.seed!(n_c * 10 + seed)
    v_bus = params.buses[1].v_bus/60

    ts_stops, opr_len, opr_width = generate_trainstop(ts_lines, params.max_opr_radius, folder_name)
    max_duration = generate_timetable(ts_lines, params.start_t, params.end_t, folder_name)
    max_duration = max(max_duration+20, opr_len/v_bus) # define the operational time for customers' timewindow generation
    cus_array = generate_customer(n_c, opr_len, opr_width, max_duration, params.detour_factor, v_bus, params.tw, ts_stops, folder_name, location)
    cgr_coords = generate_charger(ts_stops, params.α, folder_name)
    generate_bus(n_c, params.buses, params.depot, folder_name)
    depot_other(params, max_duration, folder_name)
    graph(ts_stops, ts_lines, cus_array, n_c, opr_len, opr_width, cgr_coords, folder_name, flag_annotate = 1)
end

"""
Generate demands for network like "crossing, spoon and fork"
"""
function generate(n_c::Int64, ts_lines::Vector{Userline}, params::Parameters, networkshape::Symbol; replace = 1, location = random_spread)
    upperfolder = string(networkshape) * "/"
    folder_name, seed = foldername(upperfolder, length(ts_lines), n_c, length(params.depot), length(params.buses), replace)
    Random.seed!(n_c * 10 + seed)
    v_bus = params.buses[1].v_bus/60
    ts_stops, opr_len, opr_width = read_transit_network(networkshape, upperfolder, params)
    # graph_ts(ts_stops, ts_lines, flag_annotate = 0)
    max_duration = generate_timetable(ts_lines, params.start_t, params.end_t, folder_name)
    max_duration = max(max_duration+20, opr_len/v_bus) # define the operational time for customers' timewindow generation
    cus_array = generate_customer(n_c, opr_len, opr_width, max_duration, params.detour_factor, v_bus, params.tw, ts_stops, folder_name, location)
    cgr_coords = generate_charger(ts_stops, params.α, folder_name)
    generate_bus(n_c, params.buses, params.depot, folder_name)
    depot_other(params, max_duration, folder_name)
    graph(ts_stops, ts_lines, cus_array, n_c, opr_len, opr_width, cgr_coords, folder_name, flag_annotate = 0)
end

function generate(n_c_list::Vector, ts_lines::Vector, params::Parameters, networkshape::Symbol; replace = 1, location = random_spread)
    for n_c in n_c_list
        generate(n_c, ts_lines, params, networkshape; replace = replace, location = location)
    end
end
 