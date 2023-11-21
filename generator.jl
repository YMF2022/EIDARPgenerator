include("EIDARPparams.jl")
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
generate(10, my_params)
"""
function generate(n_c::Int64, params::Parameters; upperfolder = "cross/", replace = 1)
    folder_name = foldername(upperfolder, length(params.ts_network), n_c, replace)
    Random.seed!(n_c)

    ts_coords, opr_len, opr_width = generate_trainstop(params.ts_network, params.max_opr_radius, folder_name)
    @info("Operational area: $opr_len*$opr_width km")
    max_duration = generate_timetable(params.ts_network, params.start_t, params.end_t, folder_name)
    cus_array = generate_customer(n_c, opr_len, opr_width, max_duration, params.detour_factor, params.bus_types[1].v_bus/60, folder_name)
    cgr_coords = generate_charger(ts_coords, params.α, folder_name)
    generate_bus(n_c, params.bus_types, folder_name)
    depot_other(params, max_duration, folder_name)
    graph(ts_coords, cus_array, n_c, opr_width, folder_name, flag = 0)
end