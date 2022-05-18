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
# basin_name="Pacific"
basin_name="Atlantic"
exp1 = "const_density"
exp1_dir = joinpath(OHC_datadir, "OHC", exp1)

GH19_levels_list = ["0-700", "700-2000", "2000-6000"]
levels_list = ["0-700", "700-2000", "2000-6000"]
levels_vals = [(0, 700), (700,2000), (2000, 6000)]

# levels_list = ["4000-5500"]
# levels_vals = [(4000,5500)]

GH19_file = datadir("oceanheatcontent_GH19.nc")
GH19_time = ncread(GH19_file, "time")
OHC_GH19_700 = ncread(GH19_file, "H700")
OHC_GH19_mid = ncread(GH19_file, "Hmid")
OHC_GH19_deep = ncread(GH19_file, "Hdeep")
OHC_G19 = Dict("0-700" => OHC_GH19_700, "700-2000" =>OHC_GH19_mid, "2000-6000" => OHC_GH19_deep)
tecco = collect(1992+1/24:1/12:2018)

ECCObaseidx = findall(tecco .== 1995.0416666666667)
G19baseidx =  findall(GH19_time .== 1995)

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
    Aₑ = 5.1e14
    fig, ax = plt.subplots(2, 1, figsize=(10,8))

    ylbl  = L"[W / m^2 ]"

    for i in 1:length(exps) 
        if levels in GH19_levels_list
            cropGH19t = GH19_time[end-5:end]
            cropGH19OHC = OHC_G19[levels][end-5:end]
            cropbase = OHC_G19[levels][G19baseidx]
            # cropbase = mean(cropGH19OHC)

            # ax[i].plot(cropGH19t, cropGH19OHC.- cropbase,
            # label = "GH19 Global", color = "black", alpha = 0.4, linestyle="dashed", 
            # linewidth=2.5)
        end
    end 
    for (keys,values) in shortnames
        for i in 1:length(exps) 
            exp = exps[i]
            values = sum(exp[keys][lev_idx, :], dims = 1)[:]
            if i == 1
                dts = (tecco[3:end] .- tecco[1:end-2]) .* 365 .* 24 .* 60 .* 60 
                dOHC = (1e21 / Aₑ) * (values[3:end] - values[1:end-2])
                OHC_dt =  dOHC ./ dts
                #scale by surface area of Earth
                if j == 1
                    anomalies = findall(OHC_dt .> 0.5)
                end
                ax[i].plot(tecco[2:end-1],OHC_dt,"-"*marks[keys],label = shortnames[keys], 
                linewidth=1, markersize = 2)       
            else
                values =  remove_seasonal(values,Ecycle,Fcycle) 
                dts = (tecco[13:end] .- tecco[1:end-12]) .* 365 .* 24 .* 60 .* 60 
                dOHC = (1e21 / Aₑ) .* (values[13:end] - values[1:end-12])
                OHC_dt =  dOHC ./ dts
                ax[i].plot(tecco[7:end-6],OHC_dt,"-"*marks[keys],label = shortnames[keys], 
                linewidth=1, markersize = 2)
            end
            # ax[i,j].legend()
            ax[i].set_xlim(tecco[1], tecco[end])
            ax[i].set_ylabel(ylbl)
            ax[i].grid(true)
            i == 1 ? suffix = "" : suffix = "interannual"
            ax[i].set_title(levels* "m " * suffix)
        end
    end
    fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.12)
    ax[end].legend(loc="lower left",bbox_to_anchor=(0, -0.3), ncol=6)
    ax[2].set_xlabel("calendar years")

    fig.suptitle(basin_name * " Ocean Heat Uptake ", size = 20)
    fig.tight_layout()
    isdir(plotsdir("OHC_dt")) ? nothing : mkdir(plotsdir("OHC_dt"))
    outputfile = plotsdir("OHC_dt/"*"dt_"*basin_name*"_original_and_noseasonal_" *levels* ".pdf")
    fig.savefig(outputfile)
    close("all")
end
