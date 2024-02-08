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
ts_lines = [Transitline(:straight, 5, 2.5, 30.0, 60.0, 1.0),
            Transitline(:straight, 5, 2.5, 30.0, 60.0, 1.0),
            Transitline(:circle, 9, 2.5, 30.0, 60.0, 1.0)]

# define depot(s) coordinats
depot_coord = [
    [-2.5, 0],
    [2.5, 0]
]


# define all the parameters
params = Parameters(
    buses, 
    3.0,                    # maximum operational radius around a station
    1.0,                    # maximum walking distance for each customer
    10.0,                   # maximum waiting time at transit stations
    10.0,                   # start time of operational period
    11.0,                   # end time of operational period
    depot_coord,            # location of depots
    5.1,                    # average walking speed of each customer in km/h (https://en.wikipedia.org/wiki/Preferred_walking_speed)
    1.5,                    # detour index for each customer
    3,                      # number of dummies at each charger
    50,                     # charging speed in kWh/h: Europe DC charging speed https://alternative-fuels-observatory.ec.europa.eu/general-information/recharging-systems 
    0.5,                    # service time at each stop
    15.0                    # timewindow duration
)




generate(8, ts_lines, params, :crossring, replace = 1)

demand_list = collect(5:10)
generate(demand_list, ts_lines, params, :crossring)
