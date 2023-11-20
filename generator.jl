include("EIDARPparams.jl")
include("utils.jl")

# define bus types
bus_types = [
    Bustype(
    40.0, # bus speed
    10, # bus capacity
    0.226 # energy consumption rate
    )
]

# define train lines
ts_network = [Transitline(3, 2, 5.0, 30.0, 60.0,1.0),
            Transitline(3, 2, 5.0, 30.0, 60.0,1.0)]

params = Parameters(
    ts_network,
    bus_types, 
    3.0,            # maximum operational radius around a station
    1.0,            # maximum walking distance for each customer
    10.0,           # start time of operational period
    11.0,           # end time of operational period
    [(0, 0)],       # location of depots
    0.085,          # average walking speed of each customer in km/h (https://en.wikipedia.org/wiki/Preferred_walking_speed)
    1.8,            # detour index for each customer
    3               # number of dummies at each charger
)

# Create a folder name
function foldername(upperfolder::String, n_c::Int64)
    i = 0
    folder_name = upperfolder * "c$n_c"
    if !isdir(folder_name)
        mkdir(folder_name)
    else
        target_word = "c$n_c"
        items = readdir("cross/")
        matchfolders = filter(item -> isdir(joinpath(upperfolder, item)) && contains(lowercase(item), lowercase(target_word)), items) 
        i = parse(Int, matchfolders[end][end]) + 1
        folder_name = folder_name * "-$i"
        mkdir(folder_name)
    end
end


function main(n_c::Int64, params::Parameters, upperfolder = "cross/")
    print("nothing")
end