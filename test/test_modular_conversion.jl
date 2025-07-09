using Revise
using DataFrames
using CSV
using .Pypsa2PS
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
    x_connection=0.1,
    power_factor=0.999,
    transformer_base_MVA=2000,
    config=config_df
)

println("System successfully built from $nc_file")