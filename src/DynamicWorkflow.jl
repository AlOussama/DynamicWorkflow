module DynamicWorkflow
using PowerSystems
import NetCDF
import JSON3
using DataFrames
using CSV
using Dates
using TimeSeries
# import NetCDF.readvar
using Base.Filesystem
const FS = Base.Filesystem
const PSY = PowerSystems

const BUS_PREFIX = "buses_"
const LOAD_PREFIX = "loads_"
const GEN_PREFIX = "generators_"
const LINE_PREFIX = "lines_"
const TFMR_PREFIX = "transformers_"
const STRG_PREFIX = "store_"
const STRG_UNIT_PREFIX = "storage_units_"
const BASE_MVA = 100.0
const FREQ0 = 50.0
const res = Dates.Hour(1)
const CONFIG_DF = CSV.read("mapping_config.csv",DataFrame,missingstring="NULL")
# snapshot = 230
# greet() = print("Hello World!")
function convert_system(src_file::String,snapshot=1)    
    sys = create_sys(BASE_MVA, FREQ0)
    data = open_cdf(src_file)
    # NetCDF.close(src_file)
    tv = get_nc_var(data, "snapshots")
    n_T = length(tv)
    timestamps = range(DateTime("2013-01-01T00:00:00"); step= res, length=n_T)
    sys = add_nc_buses!(sys,data)
    sys = add_nc_lines!(sys,data)
    sys = add_nc_loads!(sys, data, timestamps ,snapshot; zip_loads= false)
    sys = add_nc_generators!(sys, data,timestamps ,snapshot; config=CONFIG_DF)
    sys = add_nc_storage!(sys, data,timestamps ,snapshot; config=CONFIG_DF)
    return sys
end
function get_nc_var(data, var::String="buses_i", default = nothing)
    return haskey(data.vars, var) ? NetCDF.readvar(data[var]) : default
end
function create_sys(base_power::Float64 =1.0,freq::Float64 =50.0)
    sys = PSY.System(base_power, frequency = freq)
    return sys
end
function open_cdf(src_file::String = "data\\base_s_5_elec.nc")
    src_path = FS.joinpath(pwd(),src_file)
    println(src_path)
    data = NetCDF.open(src_path)

    return data
end

function add_nc_buses!(sys, data)
    #data = open_cdf(src_file)
    name_str = BUS_PREFIX .* "i"
    var_name = get_nc_var(data, name_str)
    n_vars = length(var_name)
    v_nom = get_nc_var(data, BUS_PREFIX*"v_nom")
    v_mag_pu_set = get_nc_var(data, BUS_PREFIX*"v_mag_pu_set", ones(n_vars))
    control = get_nc_var(data, BUS_PREFIX*"control", repeat(["PV"], n_vars))
    # control[control .== "PQ"] .= "PV"
    control[control .== "REF"] .= "PV"
    control[control .== "SLACK"] .= "PV"

    # control[10] = "SLACK"
    carrier = get_nc_var(data, BUS_PREFIX*"carrier", repeat(["AC"], n_vars))
    area = get_nc_var(data, BUS_PREFIX*"country",repeat(["Not Defined"], n_vars))

    for i=1: n_vars
        if carrier[i] =="AC"
            if !PSY.has_component(PSY.Area ,sys, area[i])
                area_c = PSY.Area(;name=area[i])
                PSY.add_component!(sys,area_c)
            else 
                area_c = PSY.get_component(PSY.Area,sys, area[i])
            end
            if i == 1
                # control[i]="REF"

            else 
                control[i] = "PV"  
            end            
            bus_t = PSY.ACBus(; 
            number = i,#parse(Int64,var_name[i])+1, 
            name = var_name[i],
            bustype = PSY.ACBusTypes[control[i]][1],
            angle = 0.0,
            magnitude = 1.0,
            area = area_c,
            voltage_limits = (min = 0.9, max = 1.1),
            base_voltage = v_nom[i],
            );
            PSY.add_component!(sys, bus_t)
        end
    end
    return sys
    
end


function add_nc_lines!(sys, data)
    name_str = LINE_PREFIX .* "i"
    var_name = get_nc_var(data, name_str)
    n_vars = length(var_name)
    available = get_nc_var(data, LINE_PREFIX*"active", repeat([true], n_vars))
    from = get_nc_var(data, LINE_PREFIX*"bus0")
    to = get_nc_var(data, LINE_PREFIX*"bus1")
    v_nom = get_nc_var(data, LINE_PREFIX*"v_nom")
    z_base = (v_nom .^2) ./ BASE_MVA
    x = get_nc_var(data, LINE_PREFIX*"x") ./z_base
    r = get_nc_var(data, LINE_PREFIX*"r") ./z_base
    b = get_nc_var(data, LINE_PREFIX*"b", repeat([0.0], n_vars)) .*z_base
    g = get_nc_var(data, LINE_PREFIX*"g", repeat([0.0], n_vars)) .*z_base
    
    rating = get_nc_var(data,LINE_PREFIX*"s_nom")./ BASE_MVA
    # available = get_nc_var(data, LINE_PREFIX*"active", av_vec)
    num_parallel = get_nc_var(data, LINE_PREFIX*"num_parallel")

    av_vec = [xi<1e6 for xi in x] .& [xi>1e-8 for xi in rating] .& [xi<1e6 for xi in r] .& [xi>0.1 for xi in num_parallel]

    l_length = get_nc_var(data,LINE_PREFIX*"length")

    for i=1: n_vars
        if av_vec[i]
            var_t = PSY.Line(; 
            name = var_name[i],
            available = av_vec[i], 
            arc = PSY.Arc(; from = PSY.get_bus(sys,from[i]), to = PSY.get_bus(sys,to[i])),
            r = r[i],
            x = x[i],
            b = (from = b[i]/2, to = b[i]/2),
            g = (from = g[i]/2, to = g[i]/2),
            rating = rating[i],
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            angle_limits = (min = -π/2, max = π/2 ),
            ext = Dict("length" => l_length[i], "num_parallel" => num_parallel[i]),
            );
            PSY.add_component!(sys,var_t)
        end
    end
    return sys
    
end


function add_nc_loads!(sys, data,timestamps ,snapshot; zip_loads)
    # zip_loads = false
    power_factor = 0.99
    zip_factors=(0.7,0.05,0.25)
    base_factor=1.0
    # snapshot=1
    prefix = LOAD_PREFIX
    tj = snapshot
    pq_factor = tan(acos(power_factor))
    name_str = prefix .* "i"
    var_name = get_nc_var(data, name_str)
    n_vars = length(var_name)
    available = get_nc_var(data, prefix*"active", repeat([true], n_vars))
    bus_name = get_nc_var(data, prefix*"bus")
    p_t = get_nc_var(data,prefix*"t_p_set")
    p_max=maximum(abs.(p_t),dims=2)
    base_powerV = BASE_MVA #
    # base_powerV = p_max .* base_factor

    p_max = p_max ./ base_powerV
    q_max = p_max .* pq_factor

    pt_ts = p_t ./ (base_powerV.*p_max)  
    q_t_pu = pt_ts .* pq_factor  

    for i=1: n_vars
        b = PSY.get_bus(sys,bus_name[i])
        load_timearray = TimeArray(timestamps, pt_ts[i,:]);
        if !zip_loads
            var_t = PSY.PowerLoad(
                name = var_name[i],
                available = available[i], 
                bus= b,
                active_power = p_t[i,tj]./base_powerV,
                reactive_power = pq_factor.*p_t[i,tj]./base_powerV,
                base_power = base_powerV,
                max_active_power = p_max[i],
                max_reactive_power = q_max[i]
                )
        else 
            
        end
        pmax_ts = SingleTimeSeries(;
           name = "max_active_power",
           data = load_timearray,
           scaling_factor_multiplier = get_max_active_power,
           );
        qmax_ts = SingleTimeSeries(;
            name = "max_reactive_power",
            data = load_timearray,
            scaling_factor_multiplier = get_max_reactive_power,
            );

        p_ts = SingleTimeSeries(;
            name = "active_power",
            data = load_timearray,
            scaling_factor_multiplier = get_max_active_power,
            );
        q_ts = SingleTimeSeries(;
             name = "reactive_power",
             data = load_timearray,
             scaling_factor_multiplier = get_max_reactive_power,
             );        

        PSY.add_component!(sys,var_t)
        PSY.add_time_series!(sys , var_t,pmax_ts)
        PSY.add_time_series!(sys , var_t, qmax_ts)
        PSY.add_time_series!(sys , var_t, p_ts)
        PSY.add_time_series!(sys , var_t, q_ts)
    end
    return sys
    
end
#ToDo add links, add stores add storage_units add generators, add loads

function add_nc_generators!(sys, data,timestamps ,snapshot; config=conf_df)
    # zip_loads = false
    # snapshot=1
    prefix = GEN_PREFIX
    tj = snapshot
    name_str = prefix .* "i"
    var_name = get_nc_var(data, name_str)
    # println(var_name)
    n_vars = length(var_name)
    available = get_nc_var(data, prefix*"active", repeat([true], n_vars))
    bus_name = get_nc_var(data, prefix*"bus")
    carrier = get_nc_var(data, prefix*"carrier")
    op_cost = get_nc_var(data, prefix*"marginal_cost")
    # cap_cost = get_nc_var(data, prefix*"capital_cost")
    # p_max = get_nc_var(data, prefix*"p_nom") 
    renew_index = get_nc_var(data, prefix*"t_p_max_pu_i")
    pmax_pu = get_nc_var(data, prefix*"t_p_max_pu")
    p_t = 1.0.*get_nc_var(data,prefix*"t_p_set") # matrix #TODO: factor wieder auf 1
    # p_max=maximum(abs.(p_t),dims=2)
    base_powerV = get_nc_var(data, prefix*"p_nom") #base_power for the generators 
    # base_powerV = p_max .* base_factor

    # p_max = p_max ./ base_powerV # should be 1 p.u
    
    pt_ts = p_t ./ (base_powerV)  
    # q_t_ts = pt_ts .* pq_nom  

    for i=1: n_vars
        bi = PSY.get_bus(sys,bus_name[i])
        ci = carrier[i]
        name_i = var_name[i] #deepcopy(var_name[i])
        # println(ci)
        config_j = config[config.pypsa_comp.==ci,:]
        # println(config_j)
        comp = config_j.component[1]
        fuel = config_j.fuel[1]
        prime_mover = config_j.prime_mover_type[1]
        pq_max = config_j.pq_max[1]
        pq_nom = config_j.pq_nom[1]
        p_timearray = TimeArray(timestamps, pt_ts[i,:]);
        q_timearray = TimeArray(timestamps, pq_nom.*pt_ts[i,:]);

        gen_st =true
        base_i = base_powerV[i];#minimum((1,maximum((0.1,1.3*pt_ts[i,tj]))))*base_powerV[i];
        p_i = pt_ts[i,tj]*base_powerV[i]/base_i;
        if abs(pt_ts[i,tj])<1e-8
            gen_st = false
        else
            gen_st = true
        end
        # cost_i = ThermalGenerationCost(CostCurve(LinearCurve(op_cost[i])),0,0,0)
        if comp=="ThermalStandard" 
            # println(name_i)
            var_t = PSY.ThermalStandard(
            name = var_name[i],
            available = available[i], 
            status = gen_st,
            bus= bi,
            active_power = maximum([p_i,0.0]),#maximum([pt_ts[i,tj],0.0]),
            reactive_power = pq_nom*p_i,#pt_ts[i,tj],
            base_power = base_i,#base_powerV[i],
            rating = sqrt(1.0+pq_max^2),
            active_power_limits = (min= 0.0,max =1.0),
            reactive_power_limits = (min= -pq_max,max =pq_max),
            ramp_limits = (up = 1.0, down = 1.0),
            operation_cost = ThermalGenerationCost(CostCurve(LinearCurve(op_cost[i])),0,0,0),
            prime_mover_type = PrimeMovers[prime_mover][1],
            fuel = ThermalFuels[fuel][1],
            );
            p_ts = SingleTimeSeries(;
                    name = "active_power",
                    data = p_timearray,
                    scaling_factor_multiplier = get_max_active_power,
                    );
            q_ts = SingleTimeSeries(;
                    name = "reactive_power",
                    data = q_timearray,
                    scaling_factor_multiplier = get_max_active_power, # should be scaled with active power, the series alrady consider the pq_scaling factor
                    );
            PSY.add_component!(sys,var_t)
            PSY.add_time_series!(sys , var_t, p_ts)
            PSY.add_time_series!(sys , var_t, q_ts)
    
        elseif ((comp == "RenewableDispatch") || (comp == "HydroDispatch"))
            # println(i)
            # println(name_i)
            # println(renew_index)

            ren_i = findfirst(isequal(name_i),renew_index)
            # var_t = nothing
            if comp=="RenewableDispatch"
                # var_t =     PSY.PowerLoad(
                #     name = var_name[i]*"_Renew",
                #     available = true, 
                #     bus= bi,
                #     active_power = -pt_ts[i,tj],
                #     reactive_power = pq_nom*-pt_ts[i,tj],
                #     base_power = base_powerV[i],
                #     max_active_power = 1,
                #     max_reactive_power = pq_max,
                #     );
                            var_t = PSY.RenewableDispatch(
                            name = var_name[i],
                            available = available[i], 
                            bus= bi,
                            active_power = maximum([pt_ts[i,tj],0]),
                            reactive_power = pq_nom*pt_ts[i,tj],
                            base_power = base_powerV[i],
                            rating = sqrt(1.0+pq_max^2),
                            power_factor = cos(atan(pq_nom)),
                            #active_power_limits = (min= 0.0,max =1.0), RenewableDispatch generators don't have active power limits, the have timeseries for this limits
                            reactive_power_limits = (min= -pq_max,max =pq_max),
                            operation_cost = RenewableGenerationCost(CostCurve(LinearCurve(op_cost[i]))),
                            prime_mover_type = PrimeMovers[prime_mover][1],
                            );
            elseif comp == "HydroDispatch"
                var_t = PSY.HydroDispatch(
                            name = var_name[i],
                            available = available[i], 
                            bus= bi,
                            active_power = maximum([pt_ts[i,tj],0]),
                            reactive_power = pq_nom*pt_ts[i,tj],
                            base_power = base_powerV[i],
                            rating = sqrt(1.0+pq_max^2),
                            # power_factor = cos(atan(pq_nom)), # HydroDispatch doesn't support power_factor
                            active_power_limits = (min= 0.0,max =1.0), 
                            reactive_power_limits = (min= -pq_max,max =pq_max),
                            ramp_limits = (up = 1.0, down = 1.0),
                            time_limits= (up= 0.0, down= 0.0),
                            operation_cost = HydroGenerationCost(CostCurve(LinearCurve(op_cost[i])),0),
                            prime_mover_type = PrimeMovers[prime_mover][1],
                            );
            else
                # var_t = PSY.RenewableDispatch(
                #     name = var_name[i],
                #     available = available[i], 
                #     bus= bi,
                #     active_power = pt_ts[i,tj],
                #     reactive_power = pq_nom*pt_ts[i,tj],
                #     base_power = base_powerV[i],
                #     rating = sqrt(1.0+pq_max^2),
                #     power_factor = cos(atan(pq_nom)),
                #     #active_power_limits = (min= 0.0,max =1.0), RenewableDispatch generators don't have active power limits, the have timeseries for this limits
                #     reactive_power_limits = (min= -pq_max,max =pq_max),
                #     operation_cost = RenewableGenerationCost(CostCurve(LinearCurve(op_cost[i]))),
                #     prime_mover_type = PrimeMovers[prime_mover][1],
                #     );                                            
            
                var_t =     PSY.PowerLoad(
                    name = var_name[i]*"_Renew",
                    available = true, 
                    bus= bi,
                    active_power = -pt_ts[i,tj],
                    reactive_power = pq_nom*-pt_ts[i,tj],
                    base_power = base_powerV[i],
                    max_active_power = 1,
                    max_reactive_power = pq_max,
                    );
            end 

            pmax_ta = TimeArray(timestamps, pmax_pu[ren_i,:]); 
            pmax_ts = SingleTimeSeries(;
                name = "max_active_power",
                data = pmax_ta,
                scaling_factor_multiplier = get_max_active_power,
                );
            qmax_ts = SingleTimeSeries(;
            name = "max_reactive_power",
            data = pmax_ta,
            scaling_factor_multiplier = get_max_reactive_power,
            );

            p_ts = SingleTimeSeries(;
            name = "active_power",
            data = p_timearray,
            scaling_factor_multiplier = get_max_active_power,
            );
            q_ts = SingleTimeSeries(;
            name = "reactive_power",
            data = q_timearray,
            scaling_factor_multiplier = get_max_active_power, # should be scaled with active power, the series alrady consider the pq_scaling factor
            );
        PSY.add_component!(sys,var_t)
        PSY.add_time_series!(sys , var_t,pmax_ts)
        PSY.add_time_series!(sys , var_t, qmax_ts)
        PSY.add_time_series!(sys , var_t, p_ts)
        PSY.add_time_series!(sys , var_t, q_ts)
        
        end

    end
    return sys
    
end



function add_nc_storage!(sys, data,timestamps ,snapshot; config=conf_df)
    # zip_loads = false
    snapshot=1
    prefix = STRG_UNIT_PREFIX
    tj = snapshot
    name_str = prefix .* "i"
    var_name = get_nc_var(data, name_str)
    # println(var_name)
    n_vars = length(var_name)
    available = get_nc_var(data, prefix*"active", repeat([true], n_vars))
    bus_name = get_nc_var(data, prefix*"bus")
    carrier = get_nc_var(data, prefix*"carrier")
    # op_cost = get_nc_var(data, prefix*"marginal_cost") # no information in some pypsa netwokrs
    # cap_cost = get_nc_var(data, prefix*"capital_cost")
    # p_max = get_nc_var(data, prefix*"p_nom") 
    renew_index = get_nc_var(data, prefix*"t_inflow_i")
    inflow_v = get_nc_var(data, prefix*"t_inflow") # must be handled like p_max_pu for generators
    base_powerV = get_nc_var(data, prefix*"p_nom") #base_power for the generators 
    cap_v = get_nc_var(data, prefix*"max_hours") #base_power for the generators
    # pmax_pu = get_nc_var(data, prefix*"t_p_max")./base_powerV 
    # pmin_pu = get_nc_var(data, prefix*"t_p_min")./base_powerV
    p_t = get_nc_var(data,prefix*"t_p_set") # matrix 
    e_t = get_nc_var(data, prefix*"t_state_of_charge") # state of charge time series in MWh
    # p_max=maximum(abs.(p_t),dims=2)
    eff_dispatch = get_nc_var(data, prefix*"efficiency_dispatch") # dispatch efficiency 
    eff_store = get_nc_var(data, prefix*"efficiency_store") # dispatch efficiency 

    # base_powerV = p_max .* base_factor

    # p_max = p_max ./ base_powerV # should be 1 p.u
    
    pt_ts = p_t ./ (base_powerV)  
    # q_t_ts = pt_ts .* pq_nom  

    for i=1: n_vars
        bi = PSY.get_bus(sys,bus_name[i])
        ci = carrier[i]
        name_i = var_name[i] #deepcopy(var_name[i])
        # println(ci)
        config_j = config[config.pypsa_comp.==ci,:]
        # println(config_j)
        comp = config_j.component[1]
        fuel = config_j.fuel[1]
        prime_mover = config_j.prime_mover_type[1]
        pq_max = config_j.pq_max[1]
        pq_nom = config_j.pq_nom[1]
        # pmax_ts = TimeArray(timestamps, pmax_pu[i,:]);
        # pmin_ts = TimeArray(timestamps, pmin_pu[i,:]); 
        p_timearray = TimeArray(timestamps, pt_ts[i,:]);
        q_timearray = TimeArray(timestamps, pq_nom.*pt_ts[i,:]);
        # cost_i = ThermalGenerationCost(CostCurve(LinearCurve(op_cost[i])),0,0,0)
        if comp == "HydroEnergyReservoir" 
            # println(name_i)
            ren_i = findfirst(isequal(name_i),renew_index)
            var_t = PSY.HydroEnergyReservoir(
                name = var_name[i],
                available = available[i], 
                bus= bi,
                active_power = maximum([pt_ts[i,tj],0]),
                reactive_power = pq_nom*pt_ts[i,tj],
                base_power = base_powerV[i],
                rating = sqrt(1.0+pq_max^2),
                storage_capacity= cap_v[i], # capacity in pu. hr 
                # power_factor = cos(atan(pq_nom)), # HydroDispatch doesn't support power_factor
                active_power_limits = (min= 0.0,max =1.0), 
                reactive_power_limits = (min= -pq_max,max =pq_max),
                ramp_limits = (up = 1.0, down = 1.0),
                time_limits= (up= 0.0, down= 0.0),
                # operation_cost = HydroGenerationCost(CostCurve(LinearCurve(op_cost[i])),0),
                prime_mover_type = PrimeMovers[prime_mover][1],
                inflow = inflow_v[ren_i,tj],
                initial_storage = e_t[i,tj],
                storage_target =  abs(e_t[i,tj]-pt_ts[i,tj]),  # should be adapted
                conversion_factor = eff_dispatch[i],
                status = true,
                );
            p_ts = SingleTimeSeries(;
                    name = "active_power",
                    data = p_timearray,
                    scaling_factor_multiplier = get_max_active_power,
                    );
            q_ts = SingleTimeSeries(;
                    name = "reactive_power",
                    data = q_timearray,
                    scaling_factor_multiplier = get_max_active_power, # should be scaled with active power, the series alrady consider the pq_scaling factor
                    );
            # pmax_timeseries = SingleTimeSeries(;
            #         name = "max_active_power",
            #         data = pmax_ts,
            #         scaling_factor_multiplier = get_max_active_power,
            #     );
            # qmax_timeseries = SingleTimeSeries(;
            #     name = "max_reactive_power",
            #     data = pmax_ts,
            #     scaling_factor_multiplier = get_max_reactive_power,
            # );
            # qmin_timeseries = SingleTimeSeries(;
            #     name = "min_reactive_power",
            #     data = pmax_ts,
            #     scaling_factor_multiplier = get_max_reactive_power,
            # );

            PSY.add_component!(sys,var_t)
            PSY.add_time_series!(sys , var_t, p_ts)
            PSY.add_time_series!(sys , var_t, q_ts)
            # PSY.add_time_series!(sys , var_t, pmax_timeseries)
            # PSY.add_time_series!(sys , var_t, qmax_timeseries)
            # PSY.add_time_series!(sys , var_t, qmin_timeseries)
        elseif comp == "HydroPumpedStorage"

            if pt_ts[i,tj] > 0.0000000000001
                hps_st = 1
            elseif pt_ts[i,tj] < -0.000000001
                hps_st = -1
            else
                hps_st = 0
            end
            # var_t = HydroPumpedStorage(;
            #     name = var_name[i],
            #     available = available[i], 
            #     bus= bi,
            #     active_power = maximum([pt_ts[i,tj],0]),
            #     reactive_power = pq_nom*pt_ts[i,tj],
            #     base_power = base_powerV[i],
            #     rating = sqrt(1.0+pq_max^2),
            #     prime_mover_type=PrimeMovers.HY,

            #     active_power_limits = (min= 0.0,max =1.0), 
            #     reactive_power_limits = (min= -pq_max,max =pq_max),
            #     ramp_limits = (up = 1.0, down = 1.0),
            #     time_limits= (up= 0.0, down= 0.0),
            #     operation_cost = HydroGenerationCost(CostCurve(LinearCurve(op_cost[i])),0),

            #     # pump attributes 
            #     rating_pump = sqrt(1.0+pq_max^2),
            #     active_power_limits_pump = (min=0.0, max=1.0),
            #     reactive_power_limits_pump = (min= -pq_max,max =pq_max),
            #     ramp_limits_pump = (up = 1.0, down = 1.0),
            #     time_limits_pump= (up = 0.0, down = 0.0),
            #     # general attributes
            #     storage_capacity=(up=cap_v[i], down=100*cap_v[i]), # the down storage is much higher than the up storage
            #     inflow=0.0,
            #     outflow=0.0,
            #     initial_storage=(up=e_t[i,tj]/base_powerV[i], down= 50*cap_v[i]),
            #     storage_target=(up=0.5*cap_v[i]*base_powerV[i], down= 50*cap_v[i]*base_powerV[i]),
            #     pump_efficiency=eff_store[i],
            #     conversion_factor = eff_dispatch[i],
            #     status=hps_st,
            #     # time_at_status=,
            #     # services=Device[],
            #     # dynamic_injector=nothing,
            #     # ext=Dict{String, Any}(),
            # );
            var_t = PSY.HydroDispatch(
                            name = var_name[i],
                            available = available[i], 
                            bus= bi,
                            active_power = maximum((pt_ts[i,tj],0.0)),
                            reactive_power = pq_nom*maximum((pt_ts[i,tj],0.0)),
                            base_power = base_powerV[i],
                            rating = sqrt(1.0+pq_max^2),
                            # power_factor = cos(atan(pq_nom)), # HydroDispatch doesn't support power_factor
                            active_power_limits = (min= 0.0,max =1.0), 
                            reactive_power_limits = (min= -pq_max,max =pq_max),
                            ramp_limits = (up = 1.0, down = 1.0),
                            time_limits= (up= 0.0, down= 0.0),
                            # operation_cost = HydroGenerationCost(CostCurve(LinearCurve(op_cost[i])),0),
                            prime_mover_type = PrimeMovers.HY,
                            );

                # var_loadt = PSY.PowerLoad(
                # name = var_name[i]*"_pump",
                # available = true, 
                # bus= bi,
                # active_power = maximum((-pt_ts[i,tj],0)),
                # reactive_power = pq_nom*maximum((-pt_ts[i,tj],0)),
                # base_power = base_powerV[i],
                # max_active_power = 1,
                # max_reactive_power = pq_max,
                # );

            # p_ts = SingleTimeSeries(;
            #         name = "active_power",
            #         data = p_timearray,
            #         scaling_factor_multiplier = get_max_active_power,
            #         );
            # q_ts = SingleTimeSeries(;
            #         name = "reactive_power",
            #         data = q_timearray,
            #         scaling_factor_multiplier = get_max_active_power, # should be scaled with active power, the series alrady consider the pq_scaling factor
            #         );
            # pmax_timeseries = SingleTimeSeries(;
            #         name = "max_active_power",
            #         data = pmax_ts,
            #         scaling_factor_multiplier = get_max_active_power,
            #     );
            # qmax_timeseries = SingleTimeSeries(;
            #     name = "max_reactive_power",
            #     data = pmax_ts,
            #     scaling_factor_multiplier = get_max_reactive_power,
            # );
            # qmin_timeseries = SingleTimeSeries(;
            #     name = "min_reactive_power",
            #     data = pmax_ts,
            #     scaling_factor_multiplier = get_max_reactive_power,
            # );

            PSY.add_component!(sys,var_t)
            # PSY.add_component!(sys,var_loadt)
            # PSY.add_time_series!(sys , var_t, p_ts)
            # PSY.add_time_series!(sys , var_t, q_ts)
            # PSY.add_time_series!(sys , var_t, pmax_timeseries)
            # PSY.add_time_series!(sys , var_t, qmax_timeseries)
            # PSY.add_time_series!(sys , var_t, qmin_timeseries)
        else 

        end

    end
    return sys
    
end



# function format_bus(src_file::AbstractString, out_path::AbstractString)
#     NetCDF.open(src_file) do data

#         name = BUS_PREFIX .* get_nc_var(data, "buses_i")
#         n_vars = length(name)
#         v_nom = get_nc_var(data, "buses_v_nom")
#         v_mag_pu_set = get_nc_var(data, "buses_v_mag_pu_set", ones(n_vars))
#         control = get_nc_var(data, "buses_control", repeat(["PV"], n_vars))
#         buses = DataFrames.DataFrame(
#             :id => 1:n_vars,
#             :name => name,
#             :v_nom => v_nom,
#             :v_mag_pu_set => v_mag_pu_set,
#             :control => control,
#             :angle => zeros(n_vars),
#         )
#         CSV.write(joinpath(out_path, "bus.csv"), buses)
#     end
# end

end # module DynamicWorkflow
