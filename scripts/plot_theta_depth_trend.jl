#plot temperature depth plots
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

theta_datadir = joinpath(datadir(), "θ_data") 
OHC_datadir = joinpath(datadir(), "OHC_data") 

basin_name="Pacific"
# basin_name="Atlantic"

exp1 = "const_density"
exp1_dir = joinpath(theta_datadir, "OHC", exp1)
basin_vol_pth = searchdir(OHC_datadir, basin_name*"Level")[1]
@load joinpath(OHC_datadir, basin_vol_pth) level_volumes
basin_volumes = level_volumes

tecco = collect(Float64, 1992+1/24:1/12:2018)

fig, ax = plt.subplots(1, 1, figsize=(9,5))

xlbl  = L"ΔT [°C]"
ax.set_xlabel(xlbl, fontsize = 13)
ax.set_ylabel("Depth [m]", fontsize = 13)
ax.grid(true)

fcycle = 1 # units: yr^{-1}
overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)

# returns "full" second
#full means that seasonal pattern has not been removed 
exp1_filename = searchdir(exp1_dir, basin_name)[2]
exp1_path = joinpath(exp1_dir, exp1_filename)

println("Comparing the two datasets")
println(exp1_path[length(theta_datadir)+1:end])

@load exp1_path θ_native
exp1_theta_native = θ_native

exps = exp1_theta_native

skip_exp = ["noIA", "129ff"]


shortnames["iter129_bulkformula"] = "iter129"
shortnames["iter0_bulkformula"] = "iter0"

k = 0
CB_color_cycle = ["#377eb8", "#ff7f00", "#4daf4a",
                  "#f781bf", "#a65628", "#984ea3",
                  "#999999", "#e41a1c", "#dede00"]
for (keys,values) in shortnames
    k +=1
    if values ∉ skip_exp
        slope = zeros(Float64, length(z))
        intercep = zeros(Float64, length(z))
        errs = zeros(Float64, 2, length(z),)
        offset = 1
        for level in 1:length(z)
            y_true = Real.(exps[keys][level, offset:end] )
            t = tecco[offset:end]
            X = ones(length(t), 2)
            X[:, 2] .= t
            # println(y_true)
            Beta = (X' * X ) \ (X' * y_true)
            intercep[level] = Beta[1]
            slope[level] = Beta[2] 
            y_pred = Beta[1] .+ (Beta[2] .* t)
            SE = standard_error(t, y_true, y_pred)
            (err1, err2) = conf_int(Beta[2], SE)
            errs[1, level] = err1 
            errs[2, level] = err2 
        end
    
    println(keys)
    println("errors for: ", z[3])
    println(errs[:, 3])
    println("Estimated Change")
    println(slope[3])

    ax.plot(slope .* 26, z, marker=string(marks[keys]) ,label = shortnames[keys],
                markersize = 5, linewidth = 2, color = CB_color_cycle[k])
    end
end
ax.set_xlim(-0.1, 0.25)
ax.set_ylim(-6000, 0)

fig.subplots_adjust(top=0.8, left=0.1, right=0.9, bottom=0.12)
legend_size = Dict("size" => 12)
ax.legend(loc="center",bbox_to_anchor=(0.5, -0.3), ncol=4, prop=legend_size)
fig.suptitle(basin_name * " Ocean Temp. Change [1992-2018]", size = 15)
fig.tight_layout()
outputfile = plotsdir("OHC/theta_depth_plot_"*basin_name*"_original_and_noseasonal" * ".pdf")
fig.savefig(outputfile)
close("all")