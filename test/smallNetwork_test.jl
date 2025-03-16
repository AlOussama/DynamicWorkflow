using PowerSystems
using PowerSimulationsDynamics
using Dates
using TimeSeries
using Sundials
using Plots
## 2 Bus System for test components
sys1 = System(100.0;frequency=50.)
bus1 = ACBus(;
           number = 1,
           name = "bus1",
           bustype = ACBusTypes.REF,
           angle = 0.0,
           magnitude = 1.0,
           voltage_limits = (min = 0.9, max = 1.05),
           base_voltage = 230.0,
       );
bus2 = ACBus(;
            number = 2,
            name = "bus2",
            bustype = ACBusTypes.PV,
            angle = 0.0,
            magnitude = 1.0,
            voltage_limits = (min = 0.9, max = 1.05),
            base_voltage = 230.0,
            );

# add_component!(sys1, bus1)
# add_component!(sys1, bus2)



# load = PowerLoad(;
#            name = "load1",
#            available = true,
#            bus = bus1,
#            active_power = 0.0, # Per-unitized by device base_power
#            reactive_power = 0.0, # Per-unitized by device base_power
#            base_power = 10.0, # MVA
#            max_active_power = 1.0, # 10 MW per-unitized by device base_power
#            max_reactive_power = 0.0,
#        );

# add_component!(sys1, load)

# solar = RenewableDispatch(;
#            name = "solar1",
#            available = true,
#            bus = bus2,
#            active_power = 0.0, # Per-unitized by device base_power
#            reactive_power = 0.0, # Per-unitized by device base_power
#            rating = 1.0, # 5 MW per-unitized by device base_power
#            prime_mover_type = PrimeMovers.PVe,
#            reactive_power_limits = (min = 0.0, max = 0.05), # 0 MVAR to 0.25 MVAR per-unitized by device base_power
#            power_factor = 1.0,
#            operation_cost = RenewableGenerationCost(nothing),
#            base_power = 5.0, # MVA
#        );


# wind1 = RenewableDispatch(;
#            name = "wind1",
#            available = true,
#            bus = bus1,
#            active_power = 1.0, # Per-unitized by device base_power
#            reactive_power = 0.1, # Per-unitized by device base_power
#            rating = 1.0, # 10 MW per-unitized by device base_power
#            prime_mover_type = PrimeMovers.WT,
#            reactive_power_limits = (min = -0.5, max = 0.5), # per-unitized by device base_power
#            power_factor = 1.0,
#            operation_cost = RenewableGenerationCost(nothing),
#            base_power = 50.0, # MVA
#        );




load1 = PowerLoad(;
           name = "load1",
           available = true,
           bus = bus1,
           active_power = 1.0, # Per-unitized by device base_power
           reactive_power = 0.1, # Per-unitized by device base_power
           base_power = 20.0, # MVA
           max_active_power = 1.0, # 10 MW per-unitized by device base_power
           max_reactive_power = 0.1,
       );


load2 = PowerLoad(;
           name = "load2",
           available = true,
           bus = bus2,
           active_power = 1, # Per-unitized by device base_power
           reactive_power = 0.1, # Per-unitized by device base_power
           base_power = 30.0, # MVA
           max_active_power = 1.0, # 10 MW per-unitized by device base_power
           max_reactive_power = 0.1,
       );


# gen11 = HydroPumpedStorage(;
#         name = "gen11",
#         available = true,
#         bus = bus1,
#         active_power = -1.0,
#         reactive_power = 0.0,
#         rating = 1.0,
#         base_power = 50.0,
#         prime_mover_type = PrimeMovers.HY,
#         active_power_limits = (min = 0.0, max = 1),
#         reactive_power_limits = (min = 0.0, max = 1),
#         ramp_limits = (up = 0.1, down = 0.1),
#         time_limits = nothing,
#         operation_cost = HydroGenerationCost(CostCurve(LinearCurve(0.15)), 0.0),
#         rating_pump = 1,
#         active_power_limits_pump = (min = 0.0, max = 1),
#         reactive_power_limits_pump = nothing,
#         ramp_limits_pump = (up = 1, down = 1),
#         time_limits_pump = nothing,
#         storage_capacity = (up = 2, down = 2),    # 2 pu * hr (10 hrs of storage)
#         inflow = 0.0,                             # Simple system with no inflow
#         outflow = 0.0,                            # Simple system with no outflow
#         initial_storage = (up = 0.0, down = 0.0), # Device with no charge at the start
#         storage_target = (up = 0.0, down = 0.0),  # Parameter outadated and does not accept nothing.
#         conversion_factor = 1.0,
#         pump_efficiency = 0.8,
# );       
add_components!(sys1, [bus1,bus2, load1, load2]) # wind1

line = Line(;
           name = "line1",
           available = true,
           active_power_flow = 0.0,
           reactive_power_flow = 0.0,
           arc = Arc(; from = bus1, to = bus2),
           r = 0.00281, # Per-unit
           x = 0.0281, # Per-unit
           b = (from = 0.00356, to = 0.00356), # Per-unit
           rating = 2.0, # Line rating of 200 MVA / System base of 100 MVA
           angle_limits = (min = -0.7, max = 0.7),
       );
add_component!(sys1, line)
# wind_values = [6.0, 7, 7, 6, 7, 9, 9, 9, 8, 8, 7, 6, 5, 5, 5, 5, 5, 6, 6, 6, 7, 6, 7, 7];
# resolution = Dates.Hour(1);
# timestamps = range(DateTime("2025-01-01T00:00:00"); step = resolution, length = 24);
# wind_timearray = TimeArray(timestamps, wind_values);
# wind_time_series = SingleTimeSeries(;
#            name = "max_active_power",
#            data = wind_timearray,
#        );
# add_time_series!(sys1, wind1, wind_time_series)


# wind_time_series_p = SingleTimeSeries(;
#            name = "active_power",
#            data = wind_timearray,
#        );
# add_time_series!(sys1, wind1, wind_time_series_p)



# load_values = [0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.5, 0.5, 0.6, 0.6,
#            0.7, 0.8, 0.8, 0.8, 0.8, 0.8, 0.9, 0.8, 0.8, 0.8, 0.8, 0.8];
# load_timearray = TimeArray(timestamps, load_values);

# load_time_series_p = SingleTimeSeries(;
#            name = "active_power",
#            data = load_timearray,
#            scaling_factor_multiplier = get_max_active_power,
#        );
# load_time_series_q = SingleTimeSeries(;
#        name = "reactive_power",
#        data = load_timearray,
#        scaling_factor_multiplier = get_max_reactive_power,
#    );


# add_time_series!(sys1,load1,load_time_series_p)
# add_time_series!(sys1,load1,load_time_series_q)


# add_time_series!(sys1,load2,load_time_series_p)
# add_time_series!(sys1,load2,load_time_series_q)

# load_time_series_p_max = SingleTimeSeries(;
#            name = "max_active_power",
#            data = load_timearray,
#            scaling_factor_multiplier = get_max_active_power,
#        );
# load_time_series_q_max = SingleTimeSeries(;
#        name = "max_reactive_power",
#        data = load_timearray,
#        scaling_factor_multiplier = get_max_reactive_power,
#    );
# add_time_series!(sys1,load1,load_time_series_p_max)
# add_time_series!(sys1,load1,load_time_series_q_max)

gen11 =    HydroDispatch(;
       name = "gen11",
       available = true,
       bus = bus1,
       active_power = -0.21,
       reactive_power = 0.0,
       rating = 1.0,
       prime_mover_type = PrimeMovers.HY,
       active_power_limits = (min = 0.0, max = 1.0),
       reactive_power_limits = (min = -0.3, max = 0.3),
       ramp_limits = (up=1.0,down=1.0),
       time_limits = (up=1.0,down=1.0),
       base_power = 100.0,
   );
   add_component!(sys1,gen11)



   gen12 =    HydroDispatch(;
   name = "gen12",
   available = true,
   bus = bus2,
   active_power = 0.75,
   reactive_power = 0.0,
   rating = 1.0,
   prime_mover_type = PrimeMovers.HY,
   active_power_limits = (min = 0.0, max = 1.0),
   reactive_power_limits = (min = -0.3, max = 0.3),
   ramp_limits = (up=1.0,down=1.0),
   time_limits = (up=1.0,down=1.0),
   base_power = 100.0,
);
add_component!(sys1,gen12)

# gas = ThermalStandard(;
#             name = "gas1",
#             available = true,
#             status = true,
#             bus = bus1,
#             active_power = -0.5, # Per-unitized by device base_power
#             reactive_power = 0.0, # Per-unitized by device base_power
#             rating = 1.0, # 30 MW per-unitized by device base_power
#             active_power_limits = (min = 0.2, max = 1.0), # 6 MW to 30 MW per-unitized by device base_power
#             reactive_power_limits = (min=-0.3,max=0.3), # Per-unitized by device base_power
#             ramp_limits = (up = 0.2, down = 0.2), # 6 MW/min up or down, per-unitized by device base_power
#             operation_cost = ThermalGenerationCost(nothing),
#             base_power = 100, # MVA
#             time_limits = (up = 8.0, down = 8.0), # Hours
#             must_run = false,
#             prime_mover_type = PrimeMovers.CC,
#             fuel = ThermalFuels.NATURAL_GAS,
#    );
# add_components!(sys1, gas)
# gas_bus2 = ThermalStandard(;
#    name = "gas2",
#    available = true,
#    status = true,
#    bus = bus2,
#    active_power = 0.0, # Per-unitized by device base_power
#    reactive_power = 0.0, # Per-unitized by device base_power
#    rating = 1.0, # 30 MW per-unitized by device base_power
#    active_power_limits = (min = 0.2, max = 1.0), # 6 MW to 30 MW per-unitized by device base_power
#    reactive_power_limits = (min=-0.3,max=0.3), # Per-unitized by device base_power
#    ramp_limits = (up = 0.2, down = 0.2), # 6 MW/min up or down, per-unitized by device base_power
#    operation_cost = ThermalGenerationCost(nothing),
#    base_power = 100, # MVA
#    time_limits = (up = 8.0, down = 8.0), # Hours
#    must_run = false,
#    prime_mover_type = PrimeMovers.CC,
#    fuel = ThermalFuels.NATURAL_GAS,
# );
# add_components!(sys1, gas2)   

# gas3 = ThermalStandard(;
#    name = "gas3",
#    available = true,
#    status = true,
#    bus = bus2,
#    active_power = 0.0, # Per-unitized by device base_power
#    reactive_power = 0.0, # Per-unitized by device base_power
#    rating = 1.0, # 30 MW per-unitized by device base_power
#    active_power_limits = (min = 0.2, max = 1.0), # 6 MW to 30 MW per-unitized by device base_power
#    reactive_power_limits = (min=-0.3,max=0.3), # Per-unitized by device base_power
#    ramp_limits = (up = 0.2, down = 0.2), # 6 MW/min up or down, per-unitized by device base_power
#    operation_cost = ThermalGenerationCost(nothing),
#    base_power = 100, # MVA
#    time_limits = (up = 8.0, down = 8.0), # Hours
#    must_run = false,
#    prime_mover_type = PrimeMovers.CC,
#    fuel = ThermalFuels.NATURAL_GAS,
# );   
# add_components!(sys1, gas3)
## machines
SauerPai = SauerPaiMachine(;
        R = 0.0,
        Xd = 0.920,
        Xq = 0.130,
        Xd_p = 0.300,
        Xq_p = 0.228,
        Xd_pp = 0.220,
        Xq_pp = 0.290,
        Xl = 0.1,
        Td0_p = 5.2,
        Tq0_p = 0.85,
        Td0_pp = 0.029,
        Tq0_pp = 0.034,
    );

Mach2_benchmark = OneDOneQMachine(;
    R = 0.0,
    Xd = 1.3125,
    Xq = 1.2578,
    Xd_p = 0.1813,
    Xq_p = 0.25,
    Td0_p = 5.89,
    Tq0_p = 0.6,
)
KundurMachine = SimpleFullMachine(;
        R = 0.003, #Example 3.1 and 4.1 of Kundur
        R_f = 0.0006,
        R_1d = 0.0284, #RD in Machowski
        R_1q = 0.0062, #RQ on Machowski
        L_d = 1.81,
        L_q = 1.76,
        L_ad = 1.66, #k*M_f or k*M_D in Machowski
        L_aq = 1.61, #k*M_Q in Machowski
        L_f1d = 1.66, #L_fD in Machowski. Assumed to be equal to L_ad
        L_ff = 1.825,
        L_1d = 0.1713, #L_D in Machowski
        L_1q = 0.7525, #L_Q in Machowski
    )

## Shaft
BaseShaft = SingleMass(; H = 5.148, D = 2.0)

## PSS 
no_pss = PSSFixed(; V_pss = 0.0)

## Governer
degov_tg = DEGOV(;
        T1 = 0.0,
        T2 = 0.0,
        T3 = 0.0,
        K = 18.0,
        T4 = 12.0,
        T5 = 5.0,
        T6 = 0.2,
        Td = 0.0,
        P_ref = 0.0,
    );
typeI_tg = TGTypeI(;
    R = 0.02,
    Ts = 0.1,
    Tc = 0.45,
    T3 = 0.0,
    T4 = 0.0,
    T5 = 50.0,
    valve_position_limits = (min = 0.3, max = 1.2),
);

###
Basic = BaseMachine(; R = 0.0, Xd_p = 0.2995, eq_p = 1.05)

BaseShaft = SingleMass(; H = 5.148, D = 10)

fixed_avr = AVRFixed(; Vf = 1.05, V_ref = 1.0)

proportional_avr = AVRSimple(; Kv = 0.01) #  Kv = 5000.0)

sexs_avr = SEXS(; Ta_Tb = 0.1, Tb = 10.0, K = 100.0, Te = 0.1, V_lim = (-4.0, 5.0))

fixed_tg = TGFixed(; efficiency = 1.0)

no_pss = PSSFixed(; V_pss = 0.0)

oneDoneQ = OneDOneQMachine(;
    R = 0.0,
    Xd = 0.8979,
    Xq = 0.646,
    Xd_p = 0.2995,
    Xq_p = 0.04,
    Td0_p = 7.4,
    Tq0_p = 0.033,
)

##
Gen1AVR = DynamicGenerator(;
        name = "gen12",
        ω_ref = 1.0,
        machine = Basic,
        shaft = BaseShaft,
        avr = fixed_avr,
        prime_mover = fixed_tg,
        pss = no_pss,
    )
Gen12AVR = DynamicGenerator(;
    name = "gen11",
    ω_ref = 1.0,
    machine = Basic,
    shaft = BaseShaft,
    avr = fixed_avr,
    prime_mover = fixed_tg,
    pss = no_pss,
)

Gen13AVR = DynamicGenerator(;
    name = "gen13",
    ω_ref = 1.0,
    machine = Basic,
    shaft = BaseShaft,
    avr = fixed_avr,
    prime_mover = fixed_tg,
    pss = no_pss,
)

gen13 = HydroEnergyReservoir(;
name = "gen13",
available = true,
bus = bus2,
active_power = 0.05,
reactive_power = 0.0,
rating = 1,
prime_mover_type = PrimeMovers.HY,
active_power_limits = (min = 0.0, max = 1.0),
reactive_power_limits = (min = -1.0, max =1.0),
ramp_limits = (up = 1.0, down = 1.0),
time_limits = nothing,
base_power = 100.0,
storage_capacity = 1.0,
inflow = 0.2,
initial_storage = 0.5,
operation_cost = HydroGenerationCost(
    CostCurve(LinearCurve(15.0)), 0.0,),
    dynamic_injector = Gen13AVR,
);
add_component!(sys1,gen13)



#     Gen2AVR = DynamicGenerator(;
#     name = "gas1",
#     ω_ref = 1.0,
#     machine = oneDoneQ,
#     shaft = BaseShaft,
#     avr = proportional_avr,
#     prime_mover = fixed_tg,
#     pss = no_pss,
# )
# Gen2AVR2 = DynamicGenerator(;
#     name = "gas2",
#     ω_ref = 1.0,
#     machine = oneDoneQ,
#     shaft = BaseShaft,
#     avr = fixed_avr,
#     prime_mover = fixed_tg,
#     pss = no_pss,
# )

# Gen2AVR3 = DynamicGenerator(;
#     name = "gas3",
#     ω_ref = 1.0,
#     machine = oneDoneQ,
#     shaft = BaseShaft,
#     avr = proportional_avr,
#     prime_mover = fixed_tg,
#     pss = no_pss,
# )
add_component!(sys1,Gen1AVR,gen12 )

add_component!(sys1,Gen12AVR,gen11)

# add_component!(sys1,Gen2AVR,gas )
# add_component!(sys1,Gen2AVR2,gas_bus2)
# add_component!(sys1,Gen2AVR3,gas3)
pert = LoadTrip(1,load2)
sim=Simulation!(ResidualModel,sys1,pwd(),(0.0,60.0),pert)
sim_Mass=Simulation!(MassMatrixModel,sys1,pwd(),(0.0,60.0),pert)


small_sig = small_signal_analysis(sim)
summary_eigenvalues(small_sig)


small_sig_Mass = small_signal_analysis(sim_Mass)
summary_eigenvalues(small_sig_Mass)


# PowerSimulationsDynamics.execute!(sim, IDA(), dtmax = 0.02, saveat = 0.02, enable_progress_bar = false)
# results = read_results(sim)

# # angle_1 = get_state_series(results, ("gas1", :ω));
# # angle_2 = get_state_series(results, ("gas2", :ω));
# # angle_3 = get_state_series(results, ("gas3", :ω));
# angle_4 = get_state_series(results, ("gen11", :ω));
# p_1 = 
# plot([angle_4], xlabel = "time", ylabel = "rotor angle [rad]", label = "angle gas gen ".*("1","2","3"))