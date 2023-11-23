include("EIDARPparams.jl")
include("utils.jl")

function generate(n_c::Int64, params::Parameters; upperfolder = "cross/", replace = 1)
    folder_name, seed = foldername(upperfolder, length(params.ts_network), n_c, replace)
    Random.seed!(n_c * 10 + seed)

    ts_coords, opr_len, opr_width = generate_trainstop(params.ts_network, params.max_opr_radius, folder_name)
    @info("Operational area: $opr_len*$opr_width km")
    max_duration = generate_timetable(params.ts_network, params.start_t, params.end_t, folder_name)
    cus_array = generate_customer(n_c, opr_len, opr_width, max_duration, params.detour_factor, params.bus_types[1].v_bus/60, folder_name)
    cgr_coords = generate_charger(ts_coords, params.α, folder_name)
    generate_bus(n_c, params.bus_types, folder_name)
    depot_other(params, max_duration, folder_name)
    graph(ts_coords, cus_array, n_c, opr_width, folder_name, flag_annotate = 0)
end
|