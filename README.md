# EIDARPgenerator
This is the repository to generate instances/scenarios for Electric Itegrated Dial-A-Ride Problem

## How to use?
Run [generator.jl](https://github.com/YMF2022/EIDARPgenerator/blob/main/generator.jl) and the  function generate()

```julia
julia> include("generator.jl")
generate

julia> generate(20, params)
```

Description of the function
```julia
julia>? generate
```

## Outputs
Output data are at folder [cross](https://github.com/YMF2022/EIDARPgenerator/tree/main/cross). For example, *l2c10* stands for 2 lines and 10 customers.
