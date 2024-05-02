include("generator.jl")

networkshape = :crossring

# define bus types
buses = [
    # type 1
    Bustype(
        25.0,   # bus speed in km/h
        15,     # bus capacity according to Sales-Lentz on-demand e-buse
        0.552,  # energy consumption kWh/km according to Sales-Lentz on-demand e-buse
        69,     # max battery capacity kWh according to Sales-Lentz on-demand e-buse
    ),

    # type 2
    Bustype(
        25.0,       # bus speed in km/h
        22,         # bus capacity
        0.552*1.5,  # energy consumption kWh/km
        69*1.5,     # max battery capacity kWh
        )
]

# define train lines
# ts_network = read_transit_network("crossring/")
ts_lines = [Transitline(:straight, 4, 5.0, 15.0, 60.0, 1.0),
            Transitline(:straight, 4, 5.0, 15.0, 60.0, 1.0)]

# define depots
depot_coord = [[-3.0, 3.0], [3.0, -3.0]]

# define chargers
chargers = [Charger(3.0, 3.0, 50.0),  
            Charger(-3.0, -3.0, 100.0)] # charging speed in kWh/h: Europe DC charging speed https://alternative-fuels-observatory.ec.europa.eu/general-information/recharging-systems 

# define all the parameters
params = Parameters(
    buses,                  # set of buses
    2.0,                    # maximum operational radius around a station
    1.0,                    # maximum walking distance for each customer
    10.0,                   # maximum waiting time at transit stations
    10.0,                   # start time of operational period
    15.0,                   # end time of operational period
    depot_coord,            # location of depots
    5.1,                    # average walking speed of each customer in km/h (https://en.wikipedia.org/wiki/Preferred_walking_speed)
    1.5,                    # detour index for each customer
    chargers,               # set of chargers
    3,                      # number of dummies at each charger
    0.5,                    # service time at each stop
    15.0                    # timewindow duration
)

annotate_offset = 0.4

generate(0, ts_lines, params, :scissor_deps; replace = 1)

# demand_list = collect(6:2:22)
# generate(demand_list, ts_lines, params, :scissor)
