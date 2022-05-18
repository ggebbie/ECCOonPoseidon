include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2, 
    PrettyTables
import CairoMakie as Mkie
import GeoMakie
using .OHC_helper
include(srcdir("config_exp.jl"))

runpath,diagpath = listexperiments(exprootdir())
# abbreviations for each experiment for labels, etc.
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments


OHC_datadir = joinpath(datadir(), "OHC_data") 

# basin_name="Pacific"
basin_name="Atlantic"

exp_type = "Total"
exp1 = "variable_density"

exp1_dir = joinpath(OHC_datadir, basin_name, exp_type, exp1)

levels_list = ["0-700", "700-2000", "2000-5500"]

level_idx = 0 
header = [key for key in values(shortnames)]

mean_data_array = zeros(3, length(header))
std_data_array = zeros(3, length(header))


for levels in levels_list
    
    level_idx += 1 
    exp1_filename = searchdir(exp1_dir, levels)[1]
    exp1_path = joinpath(exp1_dir, exp1_filename)

    println("Plotting the datasets")
    println(exp1_path[length(OHC_datadir)+1:end])

    @load exp1_path OHC_native
    exp1_OHC_native = OHC_native
    num_exps = length(keys(exp1_OHC_native))

    tecco = 1992+1/24:1/12:2018
    ylbl  = "OHC "*L"[ZJ]"


    exp_idx = 0 
    for (keys,values) in shortnames
        exp_idx += 1
        mean_data_array[level_idx, exp_idx] = mean(exp1_OHC_native[keys])
        std_data_array[level_idx, exp_idx] = std(exp1_OHC_native[keys])

    end

    
end

#can save to latex!
# println(pretty_table(data_array; header = header, row_names = levels_list))

isdir(plotsdir("OHC_stats")) ? nothing : mkdir(plotsdir("OHC_stats"))
outputfile = plotsdir("OHC_stats/"*"stats"*basin_name*"_"*exp_type*"_"*exp1*"_"*".txt")
open(outputfile, "w") do f
    write(f, "Means \n")
    pretty_table(f, mean_data_array; header = header, row_names = levels_list)
    write(f, "Stds \n")
    pretty_table(f, std_data_array; header = header, row_names = levels_list)

end

