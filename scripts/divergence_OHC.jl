include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2
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

basin_name="Atlantic"
# basin_name="Atlantic"

exp_type = "Total"
exp1 = "const_density"
exp2 = "variable_density"

exp1_dir = joinpath(OHC_datadir, basin_name, exp_type, exp1)
exp2_dir = joinpath(OHC_datadir, basin_name, exp_type, exp2)

levels_list = ["0-700", "700-2000", "2000-5500"]

for levels in levels_list
    exp1_filename = searchdir(exp1_dir, levels)[1]
    exp2_filename = searchdir(exp2_dir, levels)[1]

    exp1_path = joinpath(exp1_dir, exp1_filename)
    exp2_path = joinpath(exp2_dir, exp2_filename)

    println("Comparing the two datasets")
    println(exp1_path[length(OHC_datadir)+1:end])
    println(exp2_path[length(OHC_datadir)+1:end])

    @load exp1_path OHC_native
    exp1_OHC_native = OHC_native
    @load exp2_path OHC_native
    exp2_OHC_native = OHC_native

    tecco = 1992+1/24:1/12:2018
    ylbl  = "OHC "*L"[ZJ]"
    fig, ax = subplots(figsize=(10,5))
    for (keys,values) in shortnames
        diff = exp1_OHC_native[keys] - exp2_OHC_native[keys] 
        ax.plot(tecco,diff,"-"*marks[keys],label = shortnames[keys])
    end
    ax.grid(true)
    ax.set_title(basin_name * " OHC "*levels* "m" * 
                "\n ρ - ρ(z) ")
    ax.legend()
    ax.set_ylabel(ylbl)
    ax.set_xlabel("calendar years")

    isdir(plotsdir("OHC_diff")) ? nothing : mkdir(plotsdir("OHC_diff"))
    outputfile = plotsdir("OHC_diff/"*"diff_"*basin_name*"_"*exp_type*"_"*exp1*"_vs_"*exp2*"_"*levels*".pdf")
    fig.savefig(outputfile)
end