module Pypsa2PS

using PowerSystems
import NetCDF
using Dates
using TimeSeries
using DataFrames

const BASE_MVA = 100.0
const FREQ0 = 50.0
const N_inc_load = 1e6
const N_inc_gen = 2e6

# Utility: get variable from NetCDF with default
function get_nc_var(data, var::String="buses_i", default = nothing)
    return haskey(data.vars, var) ? NetCDF.readvar(data[var]) : default
end

# Open NetCDF file
function open_cdf(src_file::String)
    NetCDF.open(src_file)
end

# Create PowerSystems System
create_sys(base_power::Float64=BASE_MVA, freq::Float64=FREQ0) = PowerSystems.System(base_power, frequency=freq)

# General split bus and transformer function
function split_bus!(sys, orig_bus; suffix="", base_voltage=nothing, area=nothing, x_connection=0.01, transformer_base_MVA=2000, rating=1.0,N_inc=1e6)
    area = isnothing(area) ? PowerSystems.get_area(orig_bus) : area
    base_voltage = isnothing(base_voltage) ? PowerSystems.get_base_voltage(orig_bus) : base_voltage

    new_bus = PowerSystems.ACBus(
        number = N_inc + 1e3*PowerSystems.get_number(orig_bus),
        name = PowerSystems.get_name(orig_bus) * suffix,
        bustype = PowerSystems.ACBusTypes.PQ,
        angle = 0.0,
        magnitude = 1.0,
        area = area,
        voltage_limits = (min = 0.9, max = 1.1),
        base_voltage = base_voltage,
    )
    PowerSystems.add_component!(sys, new_bus)

    n_parallel_transformer = max(1, rating / transformer_base_MVA)
    x_tr = x_connection / n_parallel_transformer
    x_tr_pu = x_tr * BASE_MVA / transformer_base_MVA
    x_l = x_tr_pu
    rating_l = n_parallel_transformer * transformer_base_MVA / BASE_MVA

    var_t = PowerSystems.TapTransformer(
        tap = 1.0,
        name = PowerSystems.get_name(orig_bus) * suffix * "_connection",
        available = true,
        active_power_flow = 0,
        reactive_power_flow = 0,
        arc = PowerSystems.Arc(from = orig_bus, to = new_bus),
        x = x_l,
        r = x_l / 50,
        primary_shunt = 0.0,
        rating = rating_l,
    )
    PowerSystems.add_component!(sys, var_t)
    return new_bus, var_t
end

# Add buses
function add_nc_buses!(sys, data)
    names = get_nc_var(data, "buses_i")
    n = length(names)
    v_nom = get_nc_var(data, "buses_v_nom")
    control = get_nc_var(data, "buses_control", fill("PV", n))
    control[control .== "REF"] .= "PV"
    control[control .== "SLACK"] .= "PV"
    area = get_nc_var(data, "buses_country", fill("Not Defined", n))
    carrier = get_nc_var(data, "buses_carrier", fill("AC", n))

    for i in 1:n
        if carrier[i] == "AC"
            if !PowerSystems.has_component(PowerSystems.Area, sys, area[i])
                area_obj = PowerSystems.Area(name=area[i])
                PowerSystems.add_component!(sys, area_obj)
            else
                area_obj = PowerSystems.get_component(PowerSystems.Area, sys, area[i])
            end
            bustype = PowerSystems.ACBusTypes[control[i]][1]
            bus = PowerSystems.ACBus(
                number = i,
                name = names[i],
                bustype = bustype,
                angle = 0.0,
                magnitude = 1.0,
                area = area_obj,
                voltage_limits = (min=0.9, max=1.1),
                base_voltage = v_nom[i]
            )
            PowerSystems.add_component!(sys, bus)
        end
    end
    return sys
end

# Add lines
function add_nc_lines!(sys, data)
    names = get_nc_var(data, "lines_i")
    n = length(names)
    from = get_nc_var(data, "lines_bus0")
    to = get_nc_var(data, "lines_bus1")
    v_nom = get_nc_var(data, "lines_v_nom")
    z_base = (v_nom.^2)/BASE_MVA
    x = get_nc_var(data, "lines_x") ./z_base
    r = get_nc_var(data, "lines_r") ./z_base
    b = get_nc_var(data, "lines_b", zeros(n)) .*z_base
    g = get_nc_var(data, "lines_g", zeros(n)) .*z_base
    rating = get_nc_var(data, "lines_s_nom") ./BASE_MVA
    num_parallel = get_nc_var(data, "lines_num_parallel", ones(n))
    l_length = get_nc_var(data, "lines_length", ones(n))

    for i in 1:n
        line = PowerSystems.Line(
            name = names[i],
            available = true,
            arc = PowerSystems.Arc(from=PowerSystems.get_bus(sys, from[i]), to=PowerSystems.get_bus(sys, to[i])),
            r = r[i],
            x = x[i],
            b = (from=b[i]/2, to=b[i]/2),
            g = (from=g[i]/2, to=g[i]/2),
            rating = rating[i], #0.7 * factor for n-1 should not be included in the rating
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            angle_limits = (min=-π/2, max=π/2),
            ext = Dict("length" => l_length[i], "num_parallel" => num_parallel[i])
        )
        PowerSystems.add_component!(sys, line)
    end
    return sys
end

# Add loads, using split_bus!
function add_nc_loads!(sys, data, timestamps, snapshot, x_connection, power_factor, transformer_base_MVA; zip_loads=false)
    prefix = "loads_"
    tj = snapshot
    pq_factor = tan(acos(power_factor))
    names = get_nc_var(data, prefix * "i")
    n = length(names)
    available = get_nc_var(data, prefix * "active", fill(true, n))
    bus_name = get_nc_var(data, prefix * "bus")
    p_t = get_nc_var(data, prefix * "t_p_set")
    p_max = maximum(abs.(p_t), dims=2)
    base_powerV = BASE_MVA
    p_max = p_max ./ base_powerV
    q_max = p_max .* pq_factor
    pt_ts = p_t ./ (base_powerV .* p_max)
    q_t_pu = pt_ts .* pq_factor

    for i in 1:n
        load_timearray = TimeArray(timestamps, pt_ts[i, :])
        b_grid = PowerSystems.get_bus(sys, bus_name[i])
        if x_connection > 0
            b, _ = split_bus!(sys, b_grid;
                suffix="_load",
                base_voltage=PowerSystems.get_base_voltage(b_grid)/3,
                area=PowerSystems.get_area(b_grid),
                x_connection=x_connection,
                transformer_base_MVA=transformer_base_MVA,
                rating=p_max[i]*base_powerV,
                N_inc = N_inc_load;
            )
        else
            b = b_grid
        end
        var_t = PowerSystems.PowerLoad(
            name = names[i],
            available = available[i],
            bus = b,
            active_power = p_t[i, tj] / base_powerV,
            reactive_power = pq_factor * p_t[i, tj] / base_powerV,
            base_power = base_powerV,
            max_active_power = p_max[i],
            max_reactive_power = q_max[i]
        )
        p_ts = SingleTimeSeries(
            name = "active_power",
            data = load_timearray,
            scaling_factor_multiplier = get_max_active_power,
        )
        q_ts = SingleTimeSeries(
            name = "reactive_power",
            data = load_timearray,
            scaling_factor_multiplier = get_max_reactive_power,
        )
        PowerSystems.add_component!(sys, var_t)
        PowerSystems.add_time_series!(sys, var_t, p_ts)
        PowerSystems.add_time_series!(sys, var_t, q_ts)
    end
    return sys
end

# # Add generators (preprocessing + dispatch)
# function add_generators!(sys, data, timestamps, snapshot, config)
#     prefix = "generators_"
#     tj = snapshot
#     names = get_nc_var(data, prefix * "i")
#     n = length(names)
#     available = get_nc_var(data, prefix * "active", fill(true, n))
#     bus_name = get_nc_var(data, prefix * "bus")
#     carrier = get_nc_var(data, prefix * "carrier")
#     op_cost = get_nc_var(data, prefix * "marginal_cost")
#     renew_index = get_nc_var(data, prefix * "t_p_max_pu_i")
#     pmax_pu = get_nc_var(data, prefix * "t_p_max_pu")
#     p_t = get_nc_var(data, prefix * "t_p_set")
#     base_powerV = get_nc_var(data, prefix * "p_nom")
#     pt_ts = p_t ./ base_powerV

#     for i in 1:n
#         bi = PowerSystems.get_bus(sys, bus_name[i])
#         ci = carrier[i]
#         name_i = names[i]
#         config_j = config[config.pypsa_comp .== ci, :]
#         if isempty(config_j)
#             continue
#         end
#         comp = config_j.component[1]
#         fuel = config_j.fuel[1]
#         prime_mover = config_j.prime_mover_type[1]
#         pq_max = config_j.pq_max[1]
#         pq_nom = config_j.pq_nom[1]
#         p_timearray = TimeArray(timestamps, pt_ts[i, :])
#         q_timearray = TimeArray(timestamps, pq_nom .* pt_ts[i, :])
#         gen_st = abs(pt_ts[i, tj]) > 1e-8
#         base_i = base_powerV[i]
#         p_i = pt_ts[i, tj] * base_powerV[i] / base_i

#         if comp == "ThermalStandard"
#             add_thermal_generator!(
#                 sys, name_i, available[i], gen_st, bi, p_i, pq_nom, pq_max, base_i, op_cost[i], prime_mover, fuel, p_timearray, q_timearray
#             )
#         elseif comp == "RenewableDispatch"
#             ren_i = findfirst(isequal(name_i), renew_index)
#             add_renewable_generator!(
#                 sys, name_i, available[i], bi, pt_ts[i, tj], pq_nom, pq_max, base_powerV[i], op_cost[i], prime_mover, p_timearray, q_timearray, pmax_pu, ren_i, timestamps
#             )
#         elseif comp == "HydroDispatch"
#             ren_i = findfirst(isequal(name_i), renew_index)
#             add_hydro_dispatch_generator!(
#                 sys, name_i, available[i], bi, pt_ts[i, tj], pq_nom, pq_max, base_powerV[i], op_cost[i], prime_mover, p_timearray, q_timearray, pmax_pu, ren_i, timestamps
#             )
#         end
#     end
#     return sys
# end


function add_generators!(sys, data, timestamps, snapshot, config;
    x_connection=0.1, transformer_base_MVA=2000)
    prefix = "generators_"
    tj = snapshot
    names = get_nc_var(data, prefix * "i")
    n = length(names)
    available = get_nc_var(data, prefix * "active", fill(true, n))
    bus_name = get_nc_var(data, prefix * "bus")
    carrier = get_nc_var(data, prefix * "carrier")
    op_cost = get_nc_var(data, prefix * "marginal_cost")
    renew_index = get_nc_var(data, prefix * "t_p_max_pu_i")
    pmax_pu = get_nc_var(data, prefix * "t_p_max_pu")
    p_t = get_nc_var(data, prefix * "t_p_set")
    base_powerV = get_nc_var(data, prefix * "p_nom")
    pt_ts = p_t ./ base_powerV

    for i in 1:n
        b_grid = PowerSystems.get_bus(sys, bus_name[i])
        # Change original bus to PQ type
        b_grid.bustype = PowerSystems.ACBusTypes.PQ
        # Create new PV bus for generator connection
        b_gen, _ = split_bus!(sys, b_grid;
            suffix="_gen",
            base_voltage=PowerSystems.get_base_voltage(b_grid),
            area=PowerSystems.get_area(b_grid),
            x_connection=x_connection,
            transformer_base_MVA=transformer_base_MVA,
            rating=base_powerV[i]
        )
        # Set new bus type to PV
        b_gen.bustype = PowerSystems.ACBusTypes.PV

        ci = carrier[i]
        name_i = names[i]
        config_j = config[config.pypsa_comp .== ci, :]
        if isempty(config_j)
            continue
        end
        comp = config_j.component[1]
        fuel = config_j.fuel[1]
        prime_mover = config_j.prime_mover_type[1]
        pq_max = config_j.pq_max[1]
        pq_nom = config_j.pq_nom[1]
        p_timearray = TimeArray(timestamps, pt_ts[i, :])
        q_timearray = TimeArray(timestamps, pq_nom .* pt_ts[i, :])
        gen_st = abs(pt_ts[i, tj]) > 1e-8
        base_i = base_powerV[i]
        p_i = pt_ts[i, tj] * base_powerV[i] / base_i

        if comp == "ThermalStandard"
            gen = add_thermal_generator!(
                sys, name_i, available[i], gen_st, b_gen, p_i, pq_nom, pq_max, base_i, op_cost[i], prime_mover, fuel
            )
            p_ts = SingleTimeSeries(; name = "active_power", data = p_timearray, scaling_factor_multiplier = get_max_active_power)
            q_ts = SingleTimeSeries(; name = "reactive_power", data = q_timearray, scaling_factor_multiplier = get_max_active_power)
            PowerSystems.add_time_series!(sys, gen, p_ts)
            PowerSystems.add_time_series!(sys, gen, q_ts)
        elseif comp == "RenewableDispatch"
            ren_i = findfirst(isequal(name_i), renew_index)
            gen = add_renewable_generator!(
                sys, name_i, available[i], b_gen, pt_ts[i, tj], pq_nom, pq_max, base_powerV[i], op_cost[i], prime_mover
            )
            p_ts = SingleTimeSeries(; name = "active_power", data = p_timearray, scaling_factor_multiplier = get_max_active_power)
            q_ts = SingleTimeSeries(; name = "reactive_power", data = q_timearray, scaling_factor_multiplier = get_max_active_power)
            if !isnothing(ren_i)
                pmax_ta = TimeArray(timestamps, pmax_pu[ren_i, :])
                pmax_ts = SingleTimeSeries(; name = "max_active_power", data = pmax_ta, scaling_factor_multiplier = get_max_active_power)
                PowerSystems.add_time_series!(sys, gen, pmax_ts)
            end
            PowerSystems.add_time_series!(sys, gen, p_ts)
            PowerSystems.add_time_series!(sys, gen, q_ts)
        elseif comp == "HydroDispatch"
            ren_i = findfirst(isequal(name_i), renew_index)
            gen = add_hydro_dispatch_generator!(
                sys, name_i, available[i], b_gen, pt_ts[i, tj], pq_nom, pq_max, base_powerV[i], op_cost[i], prime_mover
            )
            p_ts = SingleTimeSeries(; name = "active_power", data = p_timearray, scaling_factor_multiplier = get_max_active_power)
            q_ts = SingleTimeSeries(; name = "reactive_power", data = q_timearray, scaling_factor_multiplier = get_max_active_power)
            if !isnothing(ren_i)
                pmax_ta = TimeArray(timestamps, pmax_pu[ren_i, :])
                pmax_ts = SingleTimeSeries(; name = "max_active_power", data = pmax_ta, scaling_factor_multiplier = get_max_active_power)
                PowerSystems.add_time_series!(sys, gen, pmax_ts)
            end
            PowerSystems.add_time_series!(sys, gen, p_ts)
            PowerSystems.add_time_series!(sys, gen, q_ts)
        end
    end
    return sys
end

function add_thermal_generator!(sys, name, available, status, bus, p, pq_nom, pq_max, base_power, op_cost, prime_mover, fuel, p_timearray, q_timearray)
    var_t = PowerSystems.ThermalStandard(
        name = name,
        available = available,
        status = status,
        bus = bus,
        active_power = max(p, 0.0),
        reactive_power = pq_nom * p,
        base_power = base_power,
        rating = 1.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = -pq_max, max = pq_max),
        ramp_limits = (up = 1.0, down = 1.0),
        operation_cost = ThermalGenerationCost(CostCurve(LinearCurve(op_cost)), 0, 0, 0),
        prime_mover_type = PrimeMovers[prime_mover][1],
        fuel = ThermalFuels[fuel][1],
    )
    p_ts = SingleTimeSeries(; name = "active_power", data = p_timearray, scaling_factor_multiplier = get_max_active_power)
    q_ts = SingleTimeSeries(; name = "reactive_power", data = q_timearray, scaling_factor_multiplier = get_max_active_power)
    PowerSystems.add_component!(sys, var_t)
    PowerSystems.add_time_series!(sys, var_t, p_ts)
    PowerSystems.add_time_series!(sys, var_t, q_ts)
end

function add_thermal_generator!(sys, name, available, status, bus, p, pq_nom, pq_max, base_power, op_cost, prime_mover, fuel)
    var_t = PowerSystems.ThermalStandard(
        name = name,
        available = available,
        status = status,
        bus = bus,
        active_power = max(p, 0.0),
        reactive_power = pq_nom * p,
        base_power = base_power,
        rating = 1.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = -pq_max, max = pq_max),
        ramp_limits = (up = 1.0, down = 1.0),
        operation_cost = ThermalGenerationCost(CostCurve(LinearCurve(op_cost)), 0, 0, 0),
        prime_mover_type = PrimeMovers[prime_mover][1],
        fuel = ThermalFuels[fuel][1],
    )
    PowerSystems.add_component!(sys, var_t)
    return var_t
end

function add_renewable_generator!(sys, name, available, bus, p, pq_nom, pq_max, base_power, op_cost, prime_mover)
    var_t = PowerSystems.RenewableDispatch(
        name = name,
        available = available,
        bus = bus,
        active_power = max(p, 0.0),
        reactive_power = pq_nom * p,
        base_power = base_power,
        rating = 1.0,
        power_factor = cos(atan(pq_nom)),
        reactive_power_limits = (min = -pq_max, max = pq_max),
        operation_cost = RenewableGenerationCost(CostCurve(LinearCurve(op_cost))),
        prime_mover_type = PrimeMovers[prime_mover][1],
    )
    PowerSystems.add_component!(sys, var_t)
    return var_t
end

function add_hydro_dispatch_generator!(sys, name, available, bus, p, pq_nom, pq_max, base_power, op_cost, prime_mover)
    var_t = PowerSystems.HydroDispatch(
        name = name,
        available = available,
        bus = bus,
        active_power = max(p, 0.0),
        reactive_power = pq_nom * p,
        base_power = base_power,
        rating = 1.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = -pq_max, max = pq_max),
        ramp_limits = (up = 1.0, down = 1.0),
        time_limits = (up = 0.0, down = 0.0),
        operation_cost = HydroGenerationCost(CostCurve(LinearCurve(op_cost)), 0),
        prime_mover_type = PrimeMovers[prime_mover][1],
    )
    PowerSystems.add_component!(sys, var_t)
    return var_t
end


# function add_renewable_generator!(sys, name, available, bus, p, pq_nom, pq_max, base_power, op_cost, prime_mover, p_timearray, q_timearray, pmax_pu, ren_i, timestamps)
#     var_t = PowerSystems.RenewableDispatch(
#         name = name,
#         available = available,
#         bus = bus,
#         active_power = max(p, 0.0),
#         reactive_power = pq_nom * p,
#         base_power = base_power,
#         rating = 1.0,
#         power_factor = cos(atan(pq_nom)),
#         reactive_power_limits = (min = -pq_max, max = pq_max),
#         operation_cost = RenewableGenerationCost(CostCurve(LinearCurve(op_cost))),
#         prime_mover_type = PrimeMovers[prime_mover][1],
#     )
#     PowerSystems.add_component!(sys, var_t)
#     p_ts = SingleTimeSeries(; name = "active_power", data = p_timearray, scaling_factor_multiplier = get_max_active_power)
#     q_ts = SingleTimeSeries(; name = "reactive_power", data = q_timearray, scaling_factor_multiplier = get_max_active_power)

#     # if !isnothing(ren_i)
#     #     pmax_ta = TimeArray(timestamps, pmax_pu[ren_i, :])
#     #     pmax_ts = SingleTimeSeries(; name = "max_active_power", data = pmax_ta, scaling_factor_multiplier = get_max_active_power)
#     #     PowerSystems.add_time_series!(sys, var_t, pmax_ts)
#     # end

#     # PowerSystems.add_time_series!(sys, var_t, p_ts)
#     # PowerSystems.add_time_series!(sys, var_t, q_ts)
# end

# function add_hydro_dispatch_generator!(sys, name, available, bus, p, pq_nom, pq_max, base_power, op_cost, prime_mover, p_timearray, q_timearray, pmax_pu, ren_i, timestamps)
#     var_t = PowerSystems.HydroDispatch(
#         name = name,
#         available = available,
#         bus = bus,
#         active_power = max(p, 0.0),
#         reactive_power = pq_nom * p,
#         base_power = base_power,
#         rating = 1.0,
#         active_power_limits = (min = 0.0, max = 1.0),
#         reactive_power_limits = (min = -pq_max, max = pq_max),
#         ramp_limits = (up = 1.0, down = 1.0),
#         time_limits = (up = 0.0, down = 0.0),
#         operation_cost = HydroGenerationCost(CostCurve(LinearCurve(op_cost)), 0),
#         prime_mover_type = PrimeMovers[prime_mover][1],
#     )
#     p_ts = SingleTimeSeries(; name = "active_power", data = p_timearray, scaling_factor_multiplier = get_max_active_power)
#     q_ts = SingleTimeSeries(; name = "reactive_power", data = q_timearray, scaling_factor_multiplier = get_max_active_power)
    
#     PowerSystems.add_component!(sys, var_t)
#     if !isnothing(ren_i)
#         pmax_ta = TimeArray(timestamps, pmax_pu[ren_i, :])
#         pmax_ts = SingleTimeSeries(; name = "max_active_power", data = pmax_ta, scaling_factor_multiplier = get_max_active_power)
#         PowerSystems.add_time_series!(sys, var_t, pmax_ts)
#     end

#     PowerSystems.add_time_series!(sys, var_t, p_ts)
#     PowerSystems.add_time_series!(sys, var_t, q_ts)
# end

# Add storages (preprocessing + dispatch)
function add_storages!(sys, data, timestamps, snapshot, config)
    prefix = "storage_units_"
    tj = snapshot
    names = get_nc_var(data, prefix * "i")
    n = length(names)
    available = get_nc_var(data, prefix * "active", fill(true, n))
    bus_name = get_nc_var(data, prefix * "bus")
    carrier = get_nc_var(data, prefix * "carrier")
    inflow_v = get_nc_var(data, prefix * "t_inflow")
    renew_index = get_nc_var(data, prefix * "t_inflow_i")
    base_powerV = get_nc_var(data, prefix * "p_nom")
    cap_v = get_nc_var(data, prefix * "max_hours")
    p_t = get_nc_var(data, prefix * "t_p_set")
    e_t = get_nc_var(data, prefix * "t_state_of_charge")
    eff_dispatch = get_nc_var(data, prefix * "efficiency_dispatch")
    eff_store = get_nc_var(data, prefix * "efficiency_store")
    pt_ts = p_t ./ base_powerV

    for i in 1:n
        bi = PowerSystems.get_bus(sys, bus_name[i])
        ci = carrier[i]
        name_i = names[i]
        config_j = config[config.pypsa_comp .== ci, :]
        if isempty(config_j)
            continue
        end
        comp = config_j.component[1]
        pq_max = config_j.pq_max[1]
        pq_nom = config_j.pq_nom[1]
        p_timearray = TimeArray(timestamps, pt_ts[i, :])
        q_timearray = TimeArray(timestamps, pq_nom .* pt_ts[i, :])

        if comp == "HydroEnergyReservoir"
            ren_i = findfirst(isequal(name_i), renew_index)
            add_hydro_energy_reservoir!(
                sys, name_i, available[i], bi, pt_ts[i, tj], pq_nom, pq_max, base_powerV[i], cap_v[i],
                config_j.prime_mover_type[1], inflow_v, ren_i, e_t, i, tj, eff_dispatch[i], p_timearray, q_timearray
            )
        elseif comp == "HydroPumpedStorage"
            print(p_timearray)
            add_hydro_pumped_storage!(
                sys, name_i, available[i], bi, pt_ts[i, tj], pq_nom, pq_max, base_powerV[i], cap_v[i],
                config_j.prime_mover_type[1], eff_dispatch[i], eff_store[i], e_t, i, tj, p_timearray, q_timearray
            )
        end
    end
    return sys
end

function add_hydro_energy_reservoir!(
    sys, name, available, bus, p, pq_nom, pq_max, base_power, cap, prime_mover,
    inflow_v, ren_i, e_t, i, tj, eff_dispatch, p_timearray, q_timearray
)
    var_t = PowerSystems.HydroEnergyReservoir(
        name = name,
        available = available,
        bus = bus,
        active_power = max(p, 0),
        reactive_power = pq_nom * p,
        base_power = base_power,
        rating = 1.0,
        storage_capacity = cap,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = -pq_max, max = pq_max),
        ramp_limits = (up = 1.0, down = 1.0),
        time_limits = (up = 0.0, down = 0.0),
        prime_mover_type = PrimeMovers[prime_mover][1],
        inflow = isnothing(ren_i) ? 0.0 : inflow_v[ren_i, tj],
        initial_storage = e_t[i, tj],
        storage_target = abs(e_t[i, tj] - p),
        conversion_factor = eff_dispatch,
        status = true,
    )
    p_ts = SingleTimeSeries(; name = "active_power", data = p_timearray, scaling_factor_multiplier = get_max_active_power)
    q_ts = SingleTimeSeries(; name = "reactive_power", data = q_timearray, scaling_factor_multiplier = get_max_active_power)
    PowerSystems.add_component!(sys, var_t)
    PowerSystems.add_time_series!(sys, var_t, p_ts)
    PowerSystems.add_time_series!(sys, var_t, q_ts)
end

function add_hydro_pumped_storage!(
    sys, name, available, bus, p, pq_nom, pq_max, base_power, cap, prime_mover,
    eff_dispatch, eff_store, e_t, i, tj, p_timearray, q_timearray
)
    # Split p_timearray into generation and pumping (load) time series
    gen_p = max.(p_timearray.values, 0.0)
    pump_p = abs.(min.(p_timearray.values, 0.0))

    # Generation time series (positive values only)
    gen_timearray = TimeArray(p_timearray.timestamp, gen_p)
    # Pumping time series (absolute value of negative values only)
    pump_timearray = TimeArray(p_timearray.timestamp, pump_p)

    # HydroDispatch generator (generation)
    var_t = PowerSystems.HydroDispatch(
        name = name,
        available = available,
        bus = bus,
        active_power = gen_p[tj],
        reactive_power = pq_nom * gen_p[tj],
        base_power = base_power,
        rating = 1.0,
        active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = -pq_max, max = pq_max),
        ramp_limits = (up = 1.0, down = 1.0),
        time_limits = (up = 0.0, down = 0.0),
        prime_mover_type = PrimeMovers.HY,
    )

    # PowerLoad (pumping)
    var_loadt = PowerSystems.PowerLoad(
        name = name * "_pump",
        available = true,
        bus = bus,
        active_power = pump_p[tj],
        reactive_power = pq_nom * pump_p[tj],
        base_power = base_power,
        max_active_power = 1,
        max_reactive_power = pq_max,
    )

    # Time series for generator
    p_ts = SingleTimeSeries(
        name = "active_power",
        data = gen_timearray,
        scaling_factor_multiplier = get_max_active_power,
    )
    q_ts = SingleTimeSeries(
        name = "reactive_power",
        data = TimeArray(gen_timearray.timestamp, pq_nom .* gen_p),
        scaling_factor_multiplier = get_max_active_power,
    )

    # Time series for load (pumping)
    p_pump_ts = SingleTimeSeries(
        name = "active_power",
        data = pump_timearray,
        scaling_factor_multiplier = get_max_active_power,
    )
    q_pump_ts = SingleTimeSeries(
        name = "reactive_power",
        data = TimeArray(pump_timearray.timestamp, pq_nom .* pump_p),
        scaling_factor_multiplier = get_max_active_power,
    )

    PowerSystems.add_component!(sys, var_t)
    PowerSystems.add_time_series!(sys, var_t, p_ts)
    PowerSystems.add_time_series!(sys, var_t, q_ts)

    PowerSystems.add_component!(sys, var_loadt)
    PowerSystems.add_time_series!(sys, var_loadt, p_pump_ts)
    PowerSystems.add_time_series!(sys, var_loadt, q_pump_ts)
end

"""
    convert_system(src_file::String; snapshot=1, x_connection=0.1, power_factor=0.999, transformer_base_MVA=2000, config=DataFrame())

Build a PowerSystems System from a pypsa NetCDF file, splitting buses for loads, generators, and storages.
"""
function convert_system(src_file::String; snapshot=1, x_connection=0.1, power_factor=0.999, transformer_base_MVA=2000, config=DataFrame())
    sys = create_sys(BASE_MVA, FREQ0)
    data = open_cdf(src_file)
    tv = get_nc_var(data, "snapshots")
    n_T = length(tv)
    timestamps = range(DateTime("2013-01-01T00:00:00"); step=Hour(1), length=n_T)
    sys = add_nc_buses!(sys, data)
    sys = add_nc_lines!(sys, data)
    sys = add_nc_loads!(sys, data, timestamps, snapshot, x_connection, power_factor, transformer_base_MVA)
    sys = add_generators!(sys, data, timestamps, snapshot, config)
    # sys = add_storages!(sys, data, timestamps, snapshot, config)
    return sys
end

end # module Pypsa2PS
