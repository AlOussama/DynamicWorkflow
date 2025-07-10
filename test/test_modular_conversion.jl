using Revise
using DataFrames
using CSV
using .Pypsa2PS
using PowerSystems
using PowerFlows
using Plots
const PS=PowerSystems
const PF=PowerFlows
# using Pypsa2PS

# Path to your pypsa NetCDF file and config file
nc_file = "data\\de_n32_solved.nc"              # Adjust path as needed
config_file = "mapping_config.csv"                # Adjust path as needed

# Read configuration DataFrame
config_df = CSV.read(config_file, DataFrame, missingstring="NULL")

# Build the PowerSystems System from the NetCDF file
sys = Pypsa2PS.convert_system(
    nc_file;
    snapshot=1,
    x_connection=0.2,
    power_factor=0.9,
    transformer_base_MVA=2000,
    config=config_df
)

println("System successfully built from $nc_file")


gens = collect(PS.get_components(Generator,sys));
loads = collect(PS.get_components(PowerLoad,sys));
buses = collect(PS.get_components(ACBus,sys));
pg = [get_active_power(x) for x in gens];
pl = [get_active_power(x) for x in loads];
@test isapprox(sum(pg),sum(pl))

pg_nom = [get_rating(x) for x in gens];
ig_max = argmax(pg_nom);
g_slack = gens[ig_max];
b_slack = get_bus(g_slack);
get_bus_numbers
set_bustype!(b_slack,"REF")

pf_res = PF.solve_powerflow(ACPowerFlow(check_reactive_power_limits=true),sys)
pf_res_bus = pf_res["bus_results"]
histogram(pf_res_bus.Vm,nbins=50)