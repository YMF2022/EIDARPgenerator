# input parameters

struct Bustype
    v_bus::Float64                  # average speed of each buse in km/h
    capacity::Int64             # bus capacity
    β::Float64                      # bus energy consumption speed kWh/km
    maxbattery::Float64             # max battery capacity in kWh
end

struct Transitline
    n_ts::Int64                     # number of stations at this line
    ts_transfer::Int64              # the transfer id
    dist_ts::Float64                # average distance between every two stations
    freq::Float64                   # frequency of this train line, e.g. every 30 mins
    v_ts::Float64                   # average speed of this line in km/h
    dt::Float64                     # dwelling time at each station
end

get_ts_number(line::Transitline) = line.n_ts

struct Parameters
    ts_network::Vector{Transitline} # transit network info 
    bus_types::Vector{Bustype}      # types of bus
    max_opr_radius::Float64         # maximum operational radius around a station
    maxwalkdist::Float64            # maximum walking distance for each customer
    start_t::Float64                # start time of operational period
    end_t::Float64                  # end time of operational period
    depot::Vector           # location of depots
    v_walk::Float64                 # average walking speed of each customer in km/h
    detour_factor::Float64          # detour index for each customer
    charger_dummies::Int64          # number of dummies at each charger
    α::Float64                      # charging speed in kWh/min
    μ::Float64                      # service time at each stop
end



