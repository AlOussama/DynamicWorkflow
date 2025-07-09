
sys_H = PowerSystems.System(100, frequency = 50.)


b = PowerSystems.ACBus(1,"1","REF", 0,1, (0.9, 1.1), 380)

add_component!(sys_H,b)

# gen_H = HydroPumpedStorage(;
#                 name = "PHS 1",
#                 available = true, 
#                 bus= b,
#                 active_power = maximum([50,0]),
#                 reactive_power = 0,
#                 base_power = 100,
#                 rating = 1,
#                 prime_mover_type=PrimeMovers.HY,
#                 active_power_limits = (min= 0.0,max =1.0), 
#                 reactive_power_limits = (min= -0.5,max =0.5),
#                 ramp_limits = (up = 1.0, down = 1.0),
#                 time_limits= (up= 0.0, down= 0.0),
#                 operation_cost = HydroGenerationCost(CostCurve(LinearCurve(1)),0),

#                 # pump attributes 
#                 rating_pump = 1,
#                 active_power_limits_pump = (min=0.0, max=1.0),
#                 reactive_power_limits_pump = (min= -0.5,max =0.5),
#                 ramp_limits_pump = (up = 1.0, down = 1.0),
#                 time_limits_pump= (up = 0.0, down = 0.0),
#                 # general attributes
#                 storage_capacity=(up=100, down=100), # the down storage is much higher than the up storage
#                 inflow=0.0,
#                 outflow=0.0,
#                 initial_storage=(up=50, down= 50),
#                 storage_target=(up=0.7, down= 0.3),
#                 pump_efficiency=0.9,
#                 conversion_factor = 0.9,
#                 status=1,
#                 # time_at_status=,
#                 # services=Device[],
#                 # dynamic_injector=nothing,
#                 # ext=Dict{String, Any}(),
#             );

gen_S =  EnergyReservoirStorage(;
        name="gen_St",
        available=true,
        bus=b,
        prime_mover_type=PowerSystems.PrimeMovers.BA,
        storage_technology_type=PowerSystems.StorageTech.OTHER_CHEM,
        storage_capacity=10.0,
        storage_level_limits=(min=0.0, max=1.0),
        initial_storage_capacity_level=0.5,
        rating=1.0,
        active_power=0.5,
        input_active_power_limits=(min=0.0, max=1.0),
        output_active_power_limits=(min=0.0, max=1.0),
        efficiency=(in=0.9, out=0.9),
        reactive_power=0.0,
        reactive_power_limits=(min=-1, max=1),
        base_power=100.0,

    )

add_component!(sys_H,gen_S)    


            load = PowerLoad(;
           name = "load 1",
           available = true,
           bus = b,
           active_power = 1.0, # Per-unitized by device base_power
           reactive_power = 0.0, # Per-unitized by device base_power
           base_power = 100.0, # MVA
           max_active_power = 1.0, # 10 MW per-unitized by device base_power
           max_reactive_power = 0.0,
       );

add_component!(sys_H,load)
    #    nf_source = Source(;
    #        name = "InfBus", #name
    #        available = true, #availability
    #        active_power = 0.0,
    #        reactive_power = 0.0,
    #        bus = b, #bus
    #        R_th = 0.0, #Rth
    #        X_th = 5e-6, #Xth
    #    );



       machine_oneDoneQ = OneDOneQMachine(;
           R = 0.0,
           Xd = 1.3125,
           Xq = 1.2578,
           Xd_p = 0.1813,
           Xq_p = 0.25,
           Td0_p = 5.89,
           Tq0_p = 0.6,
       )

       shaft_no_damping = SingleMass(;
           H = 3.01, #(M = 6.02 -> H = M/2)
           D = 1.5,
       )

    #    avr = AVRFixed(;Vf = 1.0, V_ref= 1.0)
       avr = AVRSimple(;Kv = 3)

       gov = TGFixed(;efficiency=1.0)
       pss =  PSSFixed(;V_pss=0.1)
    # pss =  PSSSimple(;K_ω = 0.01, K_p =0.1)


       dyn_gen= DynamicGenerator(name= "gen_St", ω_ref =1.0, machine=machine_oneDoneQ, shaft = shaft_no_damping, avr= avr , prime_mover = gov, pss=pss, base_power=100.0 )

       add_component!(sys_H,dyn_gen,gen_S)

       pert =  LoadChange(0.1,load, :P_ref, 0.0)

       time_span=(0.0,10.0);

       sim = Simulation!(ResidualModel, sys_H, pwd(), time_span, pert)


       PSD.execute!(sim, IDA(), dtmax = 0.01, saveat = 0.01, enable_progress_bar = true)
       results = read_results(sim)
       u1 = get_voltage_magnitude_series(results,1)
       δ_1 = get_state_series(results,("gen_St",:δ))
       ω_1 = get_state_series(results,("gen_St",:ω))
       plot(u1, xlabel = "time", ylabel = "voltage magnitude [pu]", label = "Storage")
   

       plot(δ_1, xlabel = "time", ylabel = "voltage angle [rad]", label = "Storage")
       plot(ω_1, xlabel = "time", ylabel = "frequency [rad/s]", label = "Storage")
   