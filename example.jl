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
    50,                     # charging speed in kWh/h: Europe DC charging speed https://alternative-fuels-observatory.ec.europa.eu/general-information/recharging-systems 
    0.5,                    # service time at each stop
    15.0                    # timewindow duration
)


<<<<<<< HEAD
generate(1, params; upperfolder = "cross/", replace = 1, location = closeto_ts)


# cus = collect(22:30)
# generate(cus, params; upperfolder = "TY/", replace = 1, location = random_spread)
=======
generate(10, params; upperfolder = "cross/", replace = 1)
>>>>>>> main
