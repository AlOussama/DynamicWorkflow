using Revise
using DataFrames
using CSV
using .Pypsa2PS
using PowerSystems
using PowerFlows
using Test
using Plots
const PS=PowerSystems
const PF=PowerFlows
# using Pypsa2PS

# Path to your pypsa NetCDF file and config file
# nc_file = "data\\de_n32_solved.nc"              # Adjust path as needed
nc_file = "data\\81BusGrid.nc"              # Adjust path as needed
config_file = "mapping_config.csv"                # Adjust path as needed

# Read configuration DataFrame
config_df = CSV.read(config_file, DataFrame, missingstring="NULL")

# Build the PowerSystems System from the NetCDF file
sys = Pypsa2PS.convert_system(
    nc_file;
    snapshot=5,
    x_connection=0.2,
    power_factor=0.99,
    transformer_base_MVA=2000,
    config=config_df
)

println("System successfully built from $nc_file")


gens = collect(PS.get_components(Generator,sys));

# [set_available!(x,false) for x in gens if get_active_power(x)<1e-2];
# [set_magnitude!(get_bus(x),0) for x in gens if get_active_power(x)<1e-2];
# [remove_component!(ACBus,sys, get_name(get_bus(x))) for x in gens if get_active_power(x)<1e-2];


# [remove_component!(sys, x) for x in gens if get_active_power(x)<1e-2];


loads = collect(PS.get_components(PowerLoad,sys));
buses = collect(PS.get_components(ACBus,sys));
lines = collect(PS.get_components(Line,sys));
pg = [get_active_power(x) for x in gens];
pl = [get_active_power(x) for x in loads];
# @test isapprox(sum(pg),sum(pl))

pg_nom = [get_rating(x)*get_available(x) for x in gens];
ig_max = argmax(pg_nom);
g_slack = gens[ig_max];
b_slack = get_bus(g_slack);
get_bus_numbers
set_bustype!(b_slack,"REF")

pf_res = PF.solve_powerflow(ACPowerFlow(check_reactive_power_limits=true),sys)
fl_n = [get_active_power_flow(x)/get_rating(x) for x in lines];
pf_res_bus = pf_res["bus_results"]
histogram(pf_res_bus.Vm,nbins=50)