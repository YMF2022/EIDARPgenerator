# EIDARPgenerator
![badge1](https://img.shields.io/badge/language-julia-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This is the repository to generate instances/scenarios for Electric Itegrated Dial-A-Ride Problem. 

## Usage
`generate()` is the function to generate the instances, but users need to define the following:
1. Buses
```julia
struct Bustype
    v_bus::Float64                  # average speed of each buse in km/h
    capacity::Int64                 # bus capacity
    β::Float64                      # bus energy consumption speed kWh/km
    maxbattery::Float64             # max battery capacity in kWh
end
```
2. Transit lines
```julia
struct Transitline
    n_ts::Int64                     # number of stations at this line
    ts_transfer::Int64              # the transfer id
    dist_ts::Float64                # average distance between every two stations
    freq::Float64                   # frequency of this train line, e.g. every 30 mins
    v_ts::Float64                   # average speed of this line in km/h
    dt::Float64                     # dwelling time at each station
end
```
3. And
```julia
struct Parameters
    buses::Vector{Bustype}          # set of bus
    max_opr_radius::Float64         # maximum operational radius around a station
    max_walkdist::Float64          # maximum walking distance for each customer
    max_waittime::Float64           # maximum wait time at transit station
    start_t::Float64                # start time of operational period
    end_t::Float64                  # end time of operational period
    depot::Vector                   # location of depots
    v_walk::Float64                 # average walking speed of each customer in km/h
    detour_factor::Float64          # detour index for each customer
    charger_dummies::Int64          # number of dummies at each charger
    α::Float64                      # charging speed in kWh/min
    μ::Float64                      # service time at each stop
    tw::Float64                     # timewindow duration for customers' pickup/drop-off
end
```


[example.jl](https://github.com/YMF2022/EIDARPgenerator/blob/main/example.jl) gives an example of how to use it. Users need to predefine all the parameters. For the usage of function `generate()`:
```julia
help?> generate
```

## Outputs
Output data are at folder [cross](https://github.com/YMF2022/EIDARPgenerator/tree/main/cross). Each instance has a dedicated folder. For example, folder *l2-c10-d2-bt2* stands for 2 transit lines, 10 customers, 2 depots, and 2 bus types.
