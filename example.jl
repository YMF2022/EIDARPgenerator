include("generator.jl")

# define bus types
bus_types = [
    Bustype(
    40.0,   # bus speed
    10,     # bus capacity
    0.226,  # kWh/km Volswagen ID.BUZZ https://en.wikipedia.org/wiki/Volkswagen_ID._Buzz # needs to be changed
    77,     # max battery capacity
    )
]

# define train lines
ts_network = [Transitline(3, 2, 5.0, 30.0, 60.0,1.0),
            Transitline(3, 2, 5.0, 30.0, 60.0,1.0)]


# define all the parameters
params = Parameters(
    ts_network,
    bus_types, 
    3.0,                    # maximum operational radius around a station
    1.0,                     # maximum walking distance for each customer
    10.0,                   # start time of operational period
    11.0,                   # end time of operational period
    [[0, 0], [5,0]],        # location of depots
    0.085,                  # average walking speed of each customer in km/h (https://en.wikipedia.org/wiki/Preferred_walking_speed)
    1.8,                    # detour index for each customer
    3,                      # number of dummies at each charger
    0.83,                   # charging speed in kWh/min
    0.5                     # service time at each stop
)

generate(10, params)