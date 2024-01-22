# input parameters

struct Bustype
    v_bus::Float64                  # average speed of each buse in km/h
    capacity::Int64                 # bus capacity
    β::Float64                      # bus energy consumption speed kWh/km
    maxbattery::Float64             # max battery capacity in kWh
end

const  LINE_SHAPES = Set([:straight, :circle])
abstract type Transitline end

struct Crossline <: Transitline
    shape::Symbol                   # shape of the line
    n_ts::Int64                     # number of stations at this line
    ts_transfer::Int64              # the transfer id
    dist_ts::Float64                # average distance between every two stations
    freq::Float64                   # frequency of this train line, e.g. every 30 mins
    v_ts::Float64                   # average speed of this line in km/h
    dt::Float64                     # dwelling time at each station

    function Crossline(shape, n_ts, ts_transfer, dist_ts, freq, v_ts, dt)
        if shape ∉ LINE_SHAPES
            error("The shape of line is only ':straight' or ':circle' ")
        else
            new(shape, n_ts, ts_transfer, dist_ts, freq, v_ts, dt)
        end
    end
end

struct Userline <: Transitline
    shape::Symbol                   # shape of the line
    n_ts::Int64                     # number of stations at this line
    dist_ts::Float64                # average distance between every two stations
    freq::Float64                   # frequency of this train line, e.g. every 30 mins
    v_ts::Float64                   # average speed of this line in km/h
    dt::Float64                     # dwelling time at each station

    function Userline(shape, n_ts, dist_ts, freq, v_ts, dt)
        if shape ∉ LINE_SHAPES
            error("The shape of line is only ':straight' or ':circle' ")
        else
            new(shape, n_ts, dist_ts, freq, v_ts, dt)
        end
    end
end

Transitline(shape::Symbol, n_ts::Int64, ts_transfer::Int64, dist_ts::Float64, freq::Float64, v_ts::Float64, dt::Float64) = Crossline(shape, n_ts, ts_transfer, dist_ts, freq, v_ts, dt)
Transitline(shape::Symbol, n_ts::Int64, dist_ts::Float64, freq::Float64, v_ts::Float64, dt::Float64) = Userline(shape, n_ts, dist_ts, freq, v_ts, dt)

struct Transitstop
    line::Int64
    x::Float64
    y::Float64
    transfer::Int64
end


mutable struct Parameters
    buses::Vector{Bustype}          # set of bus
    max_opr_radius::Float64         # maximum operational radius around a station
    max_walkdist::Float64           # maximum walking distance for each customer
    max_waittime::Float64           # maximum wait time at transit station
    start_t::Float64                # start time of operational period
    end_t::Float64                  # end time of operational period
    depot::Vector{Vector{Float32}}  # location of depots
    v_walk::Float64                 # average walking speed of each customer in km/h
    detour_factor::Float64          # detour index for each customer
    charger_dummies::Int64          # number of dummies at each charger
    α::Float64                      # charging speed in kWh/min
    μ::Float64                      # service time at each stop
    tw::Float64                     # timewindow duration for customers' pickup/drop-off
end



