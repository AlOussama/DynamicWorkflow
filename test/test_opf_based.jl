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
using PowerSimulations
using Ipopt
# src_net= "data\\de_n5_solved2.nc"
# TODO : save the results to the system:  set_system!(results::OptimizationProblemResults, system) 
# https://github.com/NREL-Sienna/PowerSimulations.jl/blob/3cde51f76135c993ecbf3711099083e9f2394be8/src/simulation/simulation_problem_results.jl#L197C1-L197C11
src_net= "data\\de_n32_solved.nc"

# src_net= "data\\de120.nc"
# src_net= "data\\de40.nc"

# src_net= "data\\81BusGrid.nc"
# intialization success: 
# T =1, avr=proportional_avr ==> 1 
# T =1, avr=fixed_avr ==> 1 


T_vec = 1;# 15;# 15 #1154# 3450 # 3457# 1# 2581 # 3457 #success#:4:1
N=length(T_vec)
conv = ones(N,1);
x_transformer= 0.163; # from powerfactory template "255 MVA 380/19 kV": 0.064 for 100MVA base base and 0.163 for 250MVA base   0.1; # based on the transformer s_nom # based on pypsa assumption
s_transformer= 250; #  transformer base power 2000MVA based on pypsa assumption
# sys_base = 100;

# x_connection =x_transformer*sys_base/s_transformer;#0.5e-3;#0.5e-3;#0.5e-3;#1e-3;#1e-4;#1e-3;
# # configure_logging(; console_level = Logging.Error)
power_factor= 0.99; #0.95;


sys = DynamicWorkflow.convert_system(src_net,T_vec,x_transformer,power_factor,s_transformer) # convert the pypsa network to a PowerSystems network


transform_single_time_series!(sys, Dates.Hour(8760), # horizon
                                        Dates.Hour(8760), # interval
                                               );


# problem = DecisionModel(template, sys; optimizer = solver, horizon = Hour(1))

solver = optimizer_with_attributes(Ipopt.Optimizer)
template= template_economic_dispatch(network=NetworkModel(ACPPowerModel;#LPACCPowerModel;#
use_slacks=true,
))
# edit component models
set_device_model!(template, Line, StaticBranchBounds)
# set_device_model!(template, Transformer2W, StaticBranch) # no bounds needed for equivalent transformers
# set_device_model!(template, TapTransformer, StaticBranch)
set_device_model!(template,RenewableDispatch,RenewableConstantPowerFactor) #TODO model changed to constant power factor
problem = DecisionModel(template, sys; optimizer = solver, horizon = Hour(154))
build!(problem;output_dir = mktempdir())
PowerSimulations.solve!(problem)
res = OptimizationProblemResults(problem);
res_var = read_variables(res);
res_expr = read_expressions(res);
res_param=read_parameters(res)
v_m = res_var["VoltageMagnitude__ACBus"]


# Assume v_m is a DataFrame where the first column contains time stamps
time_v = v_m[!, 1]

# Plot each other column over the time column
plt = plot(title="Voltage Magnitude over Time", xlabel="Snapshot")
for col in names(v_m)[2:end]
    plot!(plt, time_v, v_m[!, col], label=col)
end
display(plt)




res_ren_pq_factor = res_var["ReactivePowerVariable__RenewableDispatch"][!,2:end]./res_var["ActivePowerVariable__RenewableDispatch"][!,2:end]


res_conv_pq_factor = res_var["ReactivePowerVariable__ThermalStandard"][!,2:end]./res_var["ActivePowerVariable__ThermalStandard"][!,2:end]


res_ren_active_pu = res_var["ActivePowerVariable__RenewableDispatch"][!,2:end]./res_param["ActivePowerTimeSeriesParameter__RenewableDispatch"][!,2:end]
plotFlag = true;
time_v = res_var["ReactivePowerVariable__RenewableDispatch"][!,1]
if plotFlag

    plt1 = plot(title="Renewable Q/P factor over Time", xlabel="Snapshot")
    for col in names(res_ren_pq_factor)[2:end]
        plot!(plt1, time_v, res_ren_pq_factor[!, col], label=col)
    end
    display(plt1)





    plt1 = plot(title="Conv Q/P factor over Time", xlabel="Snapshot")
    for col in names(res_conv_pq_factor)[2:end]
        plot!(plt1, time_v, res_conv_pq_factor[!, col], label=col)
    end
    display(plt1)




    plt1 = plot(title="Renewable active power in pu from the maximum available", xlabel="Snapshot")
    for col in names(res_ren_active_pu)#[41:50]
        plot!(plt1, time_v, res_ren_active_pu[!, col], label=col)
    end
    display(plt1)
end