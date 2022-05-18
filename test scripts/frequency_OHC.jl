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

# basin_name="Pacific"
basin_name="Atlantic"

exp_type = "Total"
exp1 = "variable_density"

exp1_dir = joinpath(OHC_datadir, basin_name, exp_type, exp1)

levels_list = ["0-700", "700-2000", "2000-5500"]


for levels in levels_list
    exp1_filename = searchdir(exp1_dir, levels)[1]
    exp1_path = joinpath(exp1_dir, exp1_filename)

    println("Plotting the datasets")
    println(exp1_path[length(OHC_datadir)+1:end])

    @load exp1_path OHC_native
    exp1_OHC_native = OHC_native
    num_exps = length(keys(exp1_OHC_native))

    tecco = 1992+1/24:1/12:2018
    ylbl  = "OHC "*L"[ZJ]"

    fig, axs = subplots(num_exps, figsize=(8,10))
    ax_idx = 0
    for (keys,values) in shortnames
        ax_idx += 1
        Fs =  1.0 #1 month = 1/12 yrs
        x = 1.0 .* exp1_OHC_native[keys]
        (fft_freq, fft_abs) = do_FFT(x, Fs)

        ax = axs[ax_idx]
        ax.plot(fft_freq,fft_abs,"-"*marks[keys],label = shortnames[keys])
        ax.grid(true)
        ax.legend()
        ax.set_ylabel(ylbl)
    end
    fig.suptitle(basin_name * " OHC Frequency "*levels* "m")
    axs[end].set_xlabel("frequency (1/months)")

    isdir(plotsdir("OHC_freq")) ? nothing : mkdir(plotsdir("OHC_freq"))
    outputfile = plotsdir("OHC_freq/"*"freq_"*basin_name*"_"*exp_type*"_"*exp1*"_"*levels*".pdf")
    fig.savefig(outputfile)
end