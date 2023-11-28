include("generator.jl")

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
ts_network = [Transitline(3, 2, 5.0, 30.0, 60.0, 1.0),
            Transitline(3, 2, 5.0, 30.0, 60.0, 1.0)]


# define all the parameters
params = Parameters(
    ts_network,
    buses, 
    3.0,                    # maximum operational radius around a station
    1.0,                    # maximum walking distance for each customer
    10.0,                   # maximum waiting time at transit stations
    10.0,                   # start time of operational period
    11.0,                   # end time of operational period
    [[0, 0], [5, 0]],       # location of depots
    5.1,                    # average walking speed of each customer in km/h (https://en.wikipedia.org/wiki/Preferred_walking_speed)
    1.5,                    # detour index for each customer
    3,                      # number of dummies at each charger
    0.83,                   # charging speed in kWh/min
    0.5,                    # service time at each stop
    15.0                    # timewindow duration
)

generate(4, params; upperfolder = "TY/", replace = 1)
