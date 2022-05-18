#plots multiple levels of OHC 
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
exp1 = "const_density"

exp1_dir = joinpath(OHC_datadir, "OHC", exp1)
GH19_levels_list = ["0-700", "700-2000", "2000-6000"]
levels_list = ["0-700", "700-2000", "2000-6000"]
levels_vals = [(0, 700), (700,2000), (2000, 6000)]

# levels_list = ["2000-4000", "4500-bottom "]
# levels_vals = [(2000, 4000), (4500,6000)]

GH19_file = datadir("oceanheatcontent_GH19.nc")
GH19_time = ncread(GH19_file, "time")
OHC_GH19_700 = ncread(GH19_file, "H700")
OHC_GH19_mid = ncread(GH19_file, "Hmid")
OHC_GH19_deep = ncread(GH19_file, "Hdeep")
OHC_G19 = Dict("0-700" => OHC_GH19_700, "700-2000" =>OHC_GH19_mid, "2000-6000" => OHC_GH19_deep)

tecco = collect(1992+1/24:1/12:2018)
tJohnson = collect(2005+1/24:1/12:2015)

base_date = 2005 
ECCObaseidx = findall( base_date - 1/12 .< tecco .< base_date + 1/12)[1]
G19baseidx =  findall( base_date - 1/12 .< GH19_time .< base_date + 1/12)[1]
Johnsonbaseidx =  findall(tJohnson .== base_date)

fcycle = 1 # units: yr^{-1}
overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)

# my_ylims = Array{Any,1}()
my_ylims =  [(-70.37037236960289, 37.78475976166067), (-34.93573308110659, 18.957294703238247),(-9.903848550000001, 8.50626955)]
skip_exp = ["noIA", "129ff"]

shortnames["iter129_bulkformula"] = "iter129"
shortnames["iter0_bulkformula"] = "iter0"

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

    exps = [exp1_OHC_native]
    Aₑ = 5.1e14
    # month * ()
    fig, ax = plt.subplots(1, 1, figsize=(8,10))

    ylbl  = "OHC"*L"[ZJ]"
    CB_color_cycle = ["#377eb8", "#ff7f00", "#4daf4a",
    "#f781bf", "#a65628", "#984ea3",
    "#999999", "#e41a1c", "#dede00"]
    for i in 1:length(exps) 
        if levels in GH19_levels_list
            cropGH19t = GH19_time[end-5:end]
            cropGH19OHC = OHC_G19[levels][end-5:end]
            cropbase = OHC_G19[levels][G19baseidx]
            # cropbase = mean(cropGH19OHC)ß
            ax.plot(cropGH19t, cropGH19OHC.- cropbase,
            label = "GH19 Global", color = "black", alpha = 0.4, linestyle="dashed", 
            linewidth=2.5)
        end
    end 
    println(levels)
    k = 0
    for (keys,values) in shortnames
        k +=1
        if values ∉ skip_exp
            for i in 1:length(exps) 
                exp = exps[i]
                values = 1e-21 * sum(exp[keys][lev_idx, :], dims = 1)[:]

                prev_mean = mean(values)
                values =  remove_seasonal(values,Ecycle,Fcycle) .+ prev_mean
                ax.plot(tecco,values .- values[ECCObaseidx],"-"*marks[keys],label = shortnames[keys], 
                linewidth=2, markersize = 5, color = CB_color_cycle[k])
                println(keys)
                ax.set_xlim(tecco[1], tecco[end])
                ax.set_ylabel(ylbl)
                ax.grid(true)
                i == 1 ? suffix = "" : suffix = "(interannual)"
                ax.set_title(levels* "(m) " * suffix)
            end
        end
    end

    fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.12)
    legend_size = Dict("size" => 12)
    ax.legend(loc="center",bbox_to_anchor=(0.5, -0.1), ncol=5, prop=legend_size)
    ax.set_xlabel("calendar years")

    fig.suptitle(basin_name * " Ocean Heat Content ", size = 20)
    fig.tight_layout()
    
    isdir(plotsdir("OHC")) ? nothing : mkdir(plotsdir("OHC"))
    outputfile = plotsdir("OHC/"*basin_name*"_noseasonal_" *levels* ".pdf")
    fig.savefig(outputfile)
    println(outputfile)
    close("all")
end
