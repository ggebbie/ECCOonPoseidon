#plot ocean heat content trends for a specific depth-level

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2
# import CairoMakie as Mkie
# import GeoMakie
using .OHC_helper
using GoogleDrive,NCDatasets, NetCDF

include(srcdir("config_exp.jl"))
pygui(false)
runpath,diagpath = listexperiments(exprootdir())
# abbreviations for each experiment for labels, etc.
shortnames = expnames()
shortnames["iter129_bulkformula"] = "iter129"
shortnames["iter0_bulkformula"] = "iter0"
skip_exp = ["noIA", "129ff"] #IA is not interesting, 129ff follows 129bf fairly well 

marks = expsymbols()
nexp = length(shortnames) # number of experiments

OHC_datadir = joinpath(datadir(), "OHC_data") 

# basin_name="Pacific"
basin_name="Atlantic"
exp1 = "const_density"

exp1_dir = joinpath(OHC_datadir)

levels_list = ["0-700", "700-2000", "2000-5500"]
levels_vals = [(0, 700), (700,2000), (2000, 5500)]

OHC_G19, GH19_time = get_GH19()

tecco = collect(1992+1/24:1/12:2018)
ECCObaseidx = findall(tecco .== 1995.0416666666667)[1]
G19baseidx =  findall(GH19_time .== 1995)[1]
# ECCObaseidx = findall(tecco .== 2005.0416666666667)[1]
# G19baseidx =  findall(GH19_time .== 2005)[1]

fig, ax = plt.subplots(2, 3, figsize=(20,9))
ylbl  = "OHC"*L"[ZJ]"

fcycle = 1 # units: yr^{-1}
overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)
#colorblind colors 
CB_color_cycle = ["#377eb8", "#ff7f00", "#4daf4a",
"#f781bf", "#a65628", "#984ea3",
"#999999", "#e41a1c", "#dede00"]
for j in 1:length(levels_list)
    # returns "full" second
    #full means that seasonal pattern has not been removed 
    levels = levels_list[j]
    lev_pair = levels_vals[j]
    lev_idx = findall( -lev_pair[2] .< z[:].< -lev_pair[1])
    exp1_filename = searchdir(exp1_dir, basin_name)[2]
    exp1_path = joinpath(exp1_dir, exp1_filename)
    
    println("Comparing the two datasets")
    println(exp1_path[length(OHC_datadir)+1:end])

    @load exp1_path OHC_native
    exp1_OHC_native = OHC_native

    exps = [exp1_OHC_native, exp1_OHC_native]
    for i in 1:length(exps) 

        ax[i, j].plot(GH19_time[G19baseidx:end], OHC_G19[levels][G19baseidx:end] .- OHC_G19[levels][G19baseidx],
        label = "GH19 Global", color = "black", alpha = 0.75, linestyle="dashed", 
        linewidth=2.5)
    end 
    k = 0
    for (keys,values) in shortnames
        k+=1
        if values ∉ skip_exp
            for i in 1:length(exps) 
                OHC_vals = sum(exps[i][keys][lev_idx, :], dims = 1)[:] .* 1e-21
                if i == 1
                    println(values)
                    ax[i, j].plot(tecco,OHC_vals .- OHC_vals[ECCObaseidx],"-"*marks[keys],label = shortnames[keys], 
                    linewidth=2, markersize = 2, color = CB_color_cycle[k])        
                else
                    OHC_mean = mean(OHC_vals) 
                    #add mean back in after remove seasonal
                    OHC_vals =  remove_seasonal(OHC_vals,Ecycle,Fcycle) .+ OHC_mean
                    ax[i, j].plot(tecco,OHC_vals .- OHC_vals[ECCObaseidx],"-"*marks[keys],label = shortnames[keys], 
                    linewidth=2, markersize = 2, color = CB_color_cycle[k])        
                end
                ax[i,j].set_xlim(tecco[1], tecco[end])
                ax[i,j].set_ylabel(ylbl)
                ax[i,j].grid(true)
                i == 1 ? suffix = "" : suffix = "interannual"
                ax[i, j].set_title(levels* "m " * suffix)
            end
        end
    end
    ax[2,j].set_xlabel("calendar years")
end
fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.12)
ax[4].legend(loc="lower left",bbox_to_anchor=(0.1, -0.3), ncol=4)
fig.suptitle(basin_name * " Ocean Heat Content ", size = 20)
fig.tight_layout()
isdir(plotsdir("OHC")) ? nothing : mkdir(plotsdir("OHC"))
outputfile = plotsdir(basin_name*"_original_and_noseasonal_all" * ".pdf")
println("saving plots at.. " * outputfile)
fig.savefig(outputfile)
close("all")