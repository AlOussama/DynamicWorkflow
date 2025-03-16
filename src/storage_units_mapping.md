# HydroEnergyReservoir
## name::String 
## available::Bool
## bus::ACBus
    active_power::Float64
    reactive_power::Float64
    rating::Float64
    prime_mover_type::PrimeMovers
    active_power_limits::MinMax
    reactive_power_limits::Union{Nothing, MinMax}
    ramp_limits::Union{Nothing, UpDown}
    time_limits::Union{Nothing, UpDown}
    base_power::Float64
    storage_capacity::Float64  # max_hours Maximum storage capacity in the reservoir (units can be p.u-hr or m^3)
    inflow::Float64  # 0 
    initial_storage::Float64  # t_state_of_charge [:,1]
    operation_cost::Union{HydroGenerationCost, StorageCost, MarketBidCost}
    storage_target::Float64
    conversion_factor::Float64
    status::Bool
    time_at_status::Float64
    services::Vector{Service}
    dynamic_injector::Union{Nothing, DynamicInjection}
    ext::Dict{String, Any}
    internal::InfrastructureSystemsInternal
end

## Pumped Storage 
https://github.com/NREL-Sienna/PowerSystems.jl/blob/476cc982292eaf3f5de7057da239354df49440fa/src/parsers/power_system_table_data.jl#L1334

https://github.com/NREL-Sienna/HydroPowerSimulations.jl/blob/53b12aec38324cd3f120e0efca1d31077755df65/docs/src/formulation.md#HydroDispatchPumpedStorage
