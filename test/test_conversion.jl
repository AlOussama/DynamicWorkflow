using Revise
using DynamicWorkflow
using PowerSystems
using Dates
using TimeSeries
using Logging
using PowerSimulationsDynamics
using PowerFlows
using Sundials
using Plots
# src_net= "data\\de_n5_solved2.nc"
src_net= "data\\de_n32_solved.nc"

# src_net= "data\\de120.nc"
# src_net= "data\\de40.nc"

# src_net= "data\\81BusGrid.nc"
# intialization success: 
# T =1, avr=proportional_avr ==> 1 
# T =1, avr=fixed_avr ==> 1 

T_vec = 1154# 3450 # 3457# 1# 2581 # 3457 #success#:4:1
N=length(T_vec)
conv = ones(N,1);

# configure_logging(; console_level = Logging.Error)

sys = DynamicWorkflow.convert_system(src_net,T_vec)
gens = collect(get_components(ThermalStandard,sys)); #convential generators
gen_nom_bus = zeros(length(sys.bus_numbers));
gen_act_bus = zeros(length(sys.bus_numbers));
gen_max_bus = zeros(length(sys.bus_numbers));

# [gen_nom_bus[get_number(get_bus(x))]+=get_rating(x) for x in gens];
l = length(sys.bus_numbers);
[gen_max_bus[get_number(get_bus(x))]+=get_max_active_power(x) for x in gens];

[gen_act_bus[get_number(get_bus(x))]+=get_active_power(x) for x in gens];

[gen_nom_bus[get_number(get_bus(x))]+=get_rating(x) for x in gens];

[set_bustype!(x,"PQ") for x in collect(get_components(ACBus,sys))];

ind_s = sortperm(gen_act_bus,rev=true);

pv_buses = collect(get_buses(sys,Set(ind_s[1:22])));
[set_bustype!(x,"PV") for x in pv_buses];

ind_s = sortperm(gen_max_bus,rev=true);
ref_buses = collect(get_buses(sys,Set(ind_s[1])));

[set_bustype!(x,"REF") for x in ref_buses];

# for i=1: 1
    
# sys = DynamicWorkflow.convert_system(src_net,T_vec)
# println("Power Flow for snapshot $i")
# conv_i  = PowerFlows.solve_powerflow(ACPowerFlow(check_reactive_power_limits=true),sys ;show_trace=false)

# conv[i] = !ismissing(conv_i)
# # if !ismissing(conv_i)
# # ac_sol[i] = conv_i
# # end

# end


Basic = BaseMachine(; R = 0.0, Xd_p = 0.2995, eq_p = 1.05)

BaseShaft = SingleMass(; H = 1, D = 2,);

fixed_avr = AVRFixed(; Vf = 1.05, V_ref = 1.05)

proportional_avr = AVRSimple(; Kv = 1) #  Kv = 5000.0)

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
function create_dyn_gen( str_name)
return DynamicGenerator(;
        name = str_name,
        ω_ref = 1.0,
        machine = oneDoneQ,
        shaft = BaseShaft,
        avr = proportional_avr,#,fixed_avr, #
        prime_mover = fixed_tg,
        pss = no_pss,
    )
end

thermal_gens= collect(get_components(ThermalStandard,sys));

for g in thermal_gens
    name =get_name(g);
    add_component!(sys, create_dyn_gen(name),g)
end

hydro_gens= collect(get_components(HydroDispatch,sys));

hydroRes_gens= collect(get_components(HydroEnergyReservoir,sys));
for g in hydro_gens
    name =get_name(g);
    add_component!(sys, create_dyn_gen(name),g)
end

for g in hydroRes_gens
    name =get_name(g);
    add_component!(sys, create_dyn_gen(name),g)
end

# phs_gens= collect(get_components(HydroPumpedStorage,sys))
# for g in phs_gens
#     name =get_name(g);
#     add_component!(sys, create_dyn_gen(name),g)
# end



ren_gens= collect(get_components(RenewableDispatch,sys));


converter = AverageConverter(600, 40) #S_rated goes in Watts
dc_source = FixedDCSource(600.0) #Not in the original data, guessed.
filt = LCLFilter(0.08, 0.003, 0.074, 0.2, 0.01)
pll = KauraPLL(500.0, 0.084, 4.69)
virtual_H = VirtualInertia(2.0, 400.0, 20.0, 2 * pi * 50.0)
Q_control = ReactivePowerDroop(0.2, 1000.0)
outer_control = OuterControl(virtual_H, Q_control)
vsc = VoltageModeControl(0.59, 736.0, 0.0, 0.0, 0.2, 1.27, 14.3, 0.0, 50.0, 0.2)

converter_high_power() = AverageConverter(;
           rated_voltage = 380,
           rated_current = 100.0,
       )
outer_control1() = OuterControl(
        VirtualInertia(; Ta = 2.0, kd = 400.0, kω = 20.0),
        ReactivePowerDroop(; kq = 0.2, ωf = 1000.0),
    )

inner_control1() = VoltageModeControl(;
    kpv = 0.59,     #Voltage controller proportional gain
    kiv = 10,#736.0,    #Voltage controller integral gain
    kffv = 0.0,     #Binary variable enabling voltage feed-forward in current controllers
    rv = 0.0,       #Virtual resistance in pu
    lv = 0.2,       #Virtual inductance in pu
    kpc = 1.27,     #Current controller proportional gain
    kic = 14.3,     #Current controller integral gain
    kffi = 0.0,     #Binary variable enabling the current feed-forward in output of current controllers
    ωad = 50.0,     #Active damping low pass filter cut-off frequency
    kad = 0.2,      #Active damping gain
)
dc_source_lv() = FixedDCSource(; voltage = 380);
pll1() = FixedFrequency();
# pll1() = KauraPLL(;
#            ω_lp = 500.0, #Cut-off frequency for LowPass filter of PLL filter.
#            kp_pll = 0.084,  #PLL proportional gain
#            ki_pll = 4.69,   #PLL integral gain
#        );

filt1() = LCLFilter(;
           lf = 0.08,
           rf = 0.003,
           cf = 0.074,
           lg = 0.2,
           rg = 0.01,
       )
inverter = DynamicInverter(
    "TestInverter",
    1.0,
    converter,
    outer_control,
    vsc,
    dc_source,
    pll,
    filt,
);


# #  DynamicInverter(;
# #            name = get_name(gen_103),
# #            ω_ref = 1.0, # frequency reference set-point
# #            converter = converter_high_power(),
# #            outer_control = outer_control(),
# #            inner_control = inner_control(),
# #            dc_source = dc_source_lv(),
# #            freq_estimator = pll(),
# #            filter = filt(),
# #        )

       for g in ren_gens
        # name =get_name(g);
        add_component!(sys,  DynamicInverter(;
        name = get_name(g),
        ω_ref = 1.0, # frequency reference set-point
        converter = converter_high_power(),
        outer_control = outer_control1(),
        inner_control = inner_control1(),
        dc_source = dc_source_lv(),
        freq_estimator = pll1(),
        filter = filt1(),
    ),g)
    end
    
loads = collect(get_components(PowerLoad,sys));


#TODO# ac_sol  = PowerFlows.solve_powerflow(ACPowerFlow(check_reactive_power_limits=true),sys;show_trace=true);
# bus_sol_ac = ac_sol["bus_results"];


# sort!(bus_sol_ac,:P_net)
# pq_buses = get_buses(sys,Set(bus_sol_ac.bus_number[1:10]));
# [set_bustype!(b,"PQ") for b in pq_buses];

# bus_sol_ac = sort!(bus_sol_ac,:P_gen)
# pq_buses = get_buses(sys,Set(bus_sol_ac.bus_number[1:25]));
# for b in pq_buses 
#     if (get_bustype(b)==ACBusTypes.PV)
#         set_bustype!(b,"PQ") 
#     end
# end


#TODO# ac_sol  = PowerFlows.solve_powerflow(ACPowerFlow(check_reactive_power_limits=true),sys;show_trace=true);
# bus_sol_ac = ac_sol["bus_results"];
# ordered_buses = get_buses(sys,Set(bus_sol_ac.bus_number));
# vm = bus_sol_ac.Vm

# [set_magnitude!(ordered_buses[i],vm[i]) for i in 1:length(vm)];



# bus_sol_ac = sort!(bus_sol_ac,:P_gen,rev=true);
# pv_buses = get_buses(sys,Set(bus_sol_ac.bus_number[1:25]));

# for b in pv_buses 
#     if (get_bustype(b)==ACBusTypes.PV)
#         set_magnitude!(b,1.05) 
#     end
# # end

# [set_magnitude!(ordered_buses[i],vm[i]) for i in 1:length(vm)];

# ac_sol  = PowerFlows.solve_powerflow(ACPowerFlow(check_reactive_power_limits=true),sys;show_trace=true);
# bus_sol_ac = ac_sol["bus_results"];
# ordered_buses = get_buses(sys,Set(bus_sol_ac.bus_number));
# vm = bus_sol_ac.Vm
# [set_magnitude!(ordered_buses[i],vm[i]) for i in 1:length(vm)];
# [set_magnitude!(b,1.05) for b in pv_buses];


pert = LoadTrip(1,loads[2])
pert1 = LoadChange(1.0,loads[4],:P_ref,0.5);
sim=Simulation!(ResidualModel,sys,pwd(),(0.0,2.0),pert1)
# sim_Mass=Simulation!(MassMatrixModel,sys,pwd(),(0.0,20.0),pert)


small_sig = small_signal_analysis(sim)
summary_eigenvalues(small_sig)


small_sig_Mass = small_signal_analysis(sim_Mass)
summary_eigenvalues(small_sig_Mass)

execute!(sim, IDA(), dtmax = 0.02, saveat = 0.02, enable_progress_bar = true)
results = read_results(sim)
gens=collect(get_components(Generator,sys));
gen_names = get_name.(gens);
# [2,11,23,27,28]
gen_labels= ["DE0 2 lignite", "DE0 4 biomass",   "DE0 2 biomass", "DE0 3 nuclear", "DE0 0 biomass", "DE0 1 nuclear", "DE0 1 lignite"];#["DE0 0 PHS",  "DE0 3 PHS", "DE0 4 PHS"];
# gen_labels= gen_names[gen_labels_i];
# "DE0 4 PHS", "DE0 0 ror", "DE0 3 PHS", "DE0 4 ror", "DE0 1 ror", "DE0 2 PHS", "DE0 3 ror", "DE0 2 onwind", "DE0 2 offwind-dc", "DE0 4 offwind-dc", "DE0 0 solar", "DE0 3 solar", "DE0 1 solar", "DE0 2 offwind-ac", "DE0 3 onwind", "DE0 2 solar", "DE0 4 onwind", "DE0 4 offwind-ac", "DE0 1 onwind", "DE0 0 onwind", "DE0 4 solar", "DE0 0 hydro", "DE0 3 hydro", "DE0 0 nuclear", "DE0 2 lignite", "DE0 4 biomass", "DE0 3 biomass", "DE0 1 coal", "DE0 1 biomass", "DE0 0 CCGT", "DE0 1 CCGT", "DE0 0 coal", "DE0 2 coal", "DE0 2 biomass", "DE0 3 nuclear", "DE0 0 biomass", "DE0 1 nuclear", "DE0 1 lignite", "DE0 4 coal", "DE0 3 coal", "DE0 4 lignite"]
# ["DE0 0 PHS",  "DE0 2 offwind-dc","DE0 3 hydro", "DE0 0 nuclear", "DE0 2 lignite", "DE0 4 lignite"];
angle = [ get_state_series(results, (gen_labels[i], :δ)) for i=1:7];

freq = [ get_state_series(results, (gen_labels[i], :ω)) for i=1:7];
# freq = get_state_series(results, (gen_names[6], :ω));
p_ser = [get_activepower_series(results, gen_labels[i]) for i=1:7];

v_mag = [get_voltage_magnitude_series(results,i) for i=1:5];
# plot(angle, xlabel = "time", ylabel = "rotor angle [rad]");
# plot(freq, xlabel = "time", ylabel = "frequency [p.u.]", lw=2, label= permutedims(gen_labels))

plot(p_ser[5],lw=2,xlabel = "time in s", ylabel =  "active power in p.u.",fontsize=24, label=(gen_labels[5]))

plot(angle,lw=2,xlabel = "time in s", ylabel =  "angle in rad",fontsize=24, label=permutedims(gen_labels))

plot(freq,lw=2,xlabel = "time in s", ylabel =  "frequency in p.u.",fontsize=24, label=permutedims(gen_labels))
plot(v_mag,lw=2,xlabel = "time in s", ylabel =  "voltage magnitude in p.u.",fontsize=24, label="bus ".*["1" "2" "3" "4" "5"])
