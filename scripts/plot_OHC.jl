#plot ocean heat content trends for a specific depth-level

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2
import CairoMakie as Mkie
import GeoMakie
using .OHC_helper
using GoogleDrive,NCDatasets, NetCDF

include(srcdir("config_exp.jl"))
pygui(false)
runpath,diagpath = listexperiments(exprootdir())
# abbreviations for each experiment for labels, etc.
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments

OHC_datadir = joinpath(datadir(), "OHC_data") 

basin_name="Pacific"
# basin_name="Atlantic"
exp2 = "const_density"

exp1_dir = joinpath(OHC_datadir, "OHC", exp1)

levels_list = ["0-700", "700-2000", "2000-5500"]
levels_vals = [(0, 700), (700,2000), (2000, 5500)]

OHC_G19, GH19_time = get_GH19()

ECCObaseidx = findall(tecco .== 1995.0416666666667)
G19baseidx =  findall(GH19_time .== 1995)

fig, ax = plt.subplots(2, 3, figsize=(20,9))
ylbl  = "OHC"*L"[ZJ]"

tecco = collect(1992+1/24:1/12:2018)
fcycle = 1 # units: yr^{-1}
overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)

for j in 1:length(levels_list)
    # returns "full" second
    #full means that seasonal pattern has not been removed 
    levels = levels_list[j]
    lev_pair = levels_vals[j]
    lev_idx = findall( -lev_pair[2] .< z[:].< -lev_pair[1])
    exp1_filename = searchdir(exp1_dir, basin_name)[1]
    exp1_path = joinpath(exp1_dir, exp1_filename)
    
    println("Comparing the two datasets")
    println(exp1_path[length(OHC_datadir)+1:end])

    @load exp1_path OHC_native
    exp1_OHC_native = OHC_native

    exps = [exp1_OHC_native, exp1_OHC_native]

    for i in 1:length(exps) 

        ax[i, j].plot(GH19_time[end-5:end], OHC_G19[levels][end-5:end] .- baseline,
        label = "GH19 Global", color = "black", alpha = 0.4, linestyle="dashed", 
        linewidth=2.5)
    end 
    for (keys,values) in shortnames
        for i in 1:length(exps) 
            OHC_vals = sum(exps[i][keys][lev_idx, :], dims = 1)[:]
            if i == 1
                ax[i, j].plot(tecco,values .- values[ECCObaseidx],"-"*marks[keys],label = shortnames[keys], 
                linewidth=1, markersize = 2)        
            else
                OHC_mean = mean(OHC_vals)
                #add mean back in after remove seasonal
                OHC_vals =  remove_seasonal(OHC_vals,Ecycle,Fcycle) .+ val_mean
                ax[i, j].plot(tecco,values .- values[ECCObaseidx],"-"*marks[keys],label = shortnames[keys], 
                linewidth=1, markersize = 2)  
            end
            ax[i,j].set_xlim(tecco[1], tecco[end])
            ax[i,j].set_ylabel(ylbl)
            ax[i,j].grid(true)
            i == 1 ? suffix = "" : suffix = "interannual"
            ax[i, j].set_title(levels* "m " * suffix)
        end
    end
    ax[2,j].set_xlabel("calendar years")
end
fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.12)
ax[4].legend(loc="lower left",bbox_to_anchor=(0.1, -0.3), ncol=4)
fig.suptitle(basin_name * " Ocean Heat Content ", size = 20)
fig.tight_layout()
isdir(plotsdir("OHC")) ? nothing : mkdir(plotsdir("OHC"))
outputfile = plotsdir("OHC/"*basin_name*"_original_and_noseasonal_all" * ".pdf")
fig.savefig(outputfile)
close("all")