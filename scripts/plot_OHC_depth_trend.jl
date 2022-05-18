#plot OHC depth plots
include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF
using .OHC_helper
# import CairoMakie as Mkie
# import GeoMakie



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

# exp1 = "const_density"
exp1_dir = joinpath(OHC_datadir)
basin_vol_pth = searchdir(OHC_datadir, basin_name*"LevelVolumes")[1]
@load joinpath(OHC_datadir, basin_vol_pth) level_volumes
basin_volumes = level_volumes

tecco = collect(Float64, 1992+1/24:1/12:2018)

fig, ax = plt.subplots(1, 1, figsize=(9,5))

xlbl  = L"ΔH [ZJ]"
ax.set_xlabel(xlbl, fontsize = 13)
ax.set_ylabel("Depth [m]", fontsize = 13)
ax.grid(true)

fcycle = 1 # units: yr^{-1}
overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)

# returns "full" second
#full means that seasonal pattern has not been removed 
exp1_filename = searchdir(exp1_dir, "full_OHC_"*basin_name)[1]
exp1_path = joinpath(exp1_dir, exp1_filename)

println("Comparing the two datasets")
println(exp1_path[length(theta_datadir)+1:end])

@load exp1_path OHC_native
exp1_OHC_native = OHC_native

exps = exp1_OHC_native

skip_exp = ["noIA", "129ff"]
Aₑ = 5.1e14
ρ = 1029
cₚ = 3085

shortnames["iter129_bulkformula"] = "iter129"
shortnames["iter0_bulkformula"] = "iter0"
k = 0

#colorblind colors 
CB_color_cycle = ["#377eb8", "#ff7f00", "#4daf4a",
"#f781bf", "#a65628", "#984ea3",
"#999999", "#e41a1c", "#dede00"]

for (keys,values) in shortnames
    k +=1
    if values ∉ skip_exp
        slopes = zeros(Float64, length(z))
        intercepts = zeros(Float64, length(z))
        errs = zeros(Float64, 2, length(z),)
        offset = 1
        for level in 1:length(z)
            y_true = Real.(exps[keys][level, offset:end] ) .* 1e-21
            t = tecco[offset:end]
            X = ones(length(t), 2)
            X[:, 2] .= t
            β = (X' * X ) \ (X' * y_true)
            intercepts[level] = β[1]
            slopes[level] = β[2] 

            y_pred = intercepts[level] .+ (slopes[level] .* t)
            SE = standard_error(t, y_true, y_pred)
            (err1, err2) = conf_int(slopes[level], SE)
            errs[1, level] = err1 
            errs[2, level] = err2 
        end

        ax.plot(slopes .* 26, z, marker=string(marks[keys]) ,label = shortnames[keys],
                markersize = 5, linewidth = 2, color = CB_color_cycle[k])
    end
end
println(ax.get_xlim())
ax.set_ylim(-6000, 0)

fig.subplots_adjust(top=0.8, left=0.1, right=0.9, bottom=0.12)
legend_size = Dict("size" => 12)
ax.legend(loc="center",bbox_to_anchor=(0.5, -0.3), ncol=4, prop=legend_size)
fig.suptitle(basin_name * " Ocean Heat Content Change [1992-2018]", size = 15)
fig.tight_layout()
outputfile = plotsdir("OHC/OHC_depth_plot_"*basin_name*"_original_and_noseasonal" * ".pdf")
fig.savefig(outputfile)
close("all")