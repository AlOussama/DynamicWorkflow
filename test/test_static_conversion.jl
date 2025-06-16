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

T_vec = 15 # 533#67#35#15 #1154# 3450 # 3457# 1# 2581 # 3457 #success#:4:1
N=length(T_vec)
conv = ones(N,1);
x_connection =2e-3;
# configure_logging(; console_level = Logging.Error)

sys = DynamicWorkflow.convert_system(src_net,T_vec,x_connection) # convert the pypsa network to a PowerSystems network

gens = collect(get_components(ThermalStandard,sys)); #get list of convential generators
gen_nom_bus = zeros(length(sys.bus_numbers)); # initialize series for the generation capacity of the buses 
gen_act_bus = zeros(length(sys.bus_numbers)); # initialize series for the actual generation at time step t of the buses
gen_max_bus = zeros(length(sys.bus_numbers)); # initialize series for the maximum generation at time step t of the buses

l = length(sys.bus_numbers);
[gen_max_bus[get_number(get_bus(x))]+=get_max_active_power(x) for x in gens]; # sum the maximum generation at time step t per bus
[gen_act_bus[get_number(get_bus(x))]+=get_active_power(x) for x in gens]; # sum the actual generation at time step t per bus 
[gen_nom_bus[get_number(get_bus(x))]+=get_rating(x) for x in gens];# sum the the rated generation capacity per bus
[set_bustype!(x,"PQ") for x in collect(get_components(ACBus,sys))]; # set all buses to PQ buses

ind_s = sortperm(gen_act_bus,rev=true); # sort the actual generation vectors in descending order
pv_buses = collect(get_buses(sys,Set(ind_s[1:20]))); # alternatively 1:l÷3*2
[set_bustype!(x,"PV") for x in pv_buses]; # set the 22- buses with the highest actural generation as pv buses


ind_s = sortperm(gen_max_bus,rev=true); # sort the maximum generation at time step t in descending order
ref_buses = collect(get_buses(sys,Set(ind_s[1]))); # set the bus with the highest possible generation to a slack bus
[set_bustype!(x,"REF") for x in ref_buses];


loads = collect(get_components(PowerLoad,sys));
loads_p = [get_active_power(l) for l in loads]; 
loads_i_sort = sortperm(loads_p,rev=false);

sum([get_active_power(g) for g in loads])*100
gens = collect(get_components(Generator,sys));
sum([get_active_power(g) for g in gens])*100


ac_sol = PowerFlows.solve_powerflow(ACPowerFlow(check_reactive_power_limits=true),sys)
histogram(ac_sol["bus_results"].Vm,nbins=30)