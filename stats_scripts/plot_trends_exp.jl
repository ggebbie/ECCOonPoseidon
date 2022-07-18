#plot OHC depth plots
include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF
using .OHC_helper

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

exp1_dir = joinpath(OHC_datadir)
basin_vol_pth = searchdir(OHC_datadir, basin_name*"LevelVolumes")[1]
@load joinpath(OHC_datadir, basin_vol_pth) level_volumes
basin_volumes = level_volumes

tecco = collect(Float64, 1992+1/24:1/12:2018)

fig, ax = plt.subplots(1, 1, figsize=(9,5))

xlbl  = L"Experiment"
ax.set_xlabel(xlbl, fontsize = 13)
ax.set_ylabel("ΔH [ZJ]", fontsize = 13)
ax.grid(true)

fcycle = 1 # units: yr^{-1}
overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)

# returns "full" second
#full means that seasonal pattern has not been removed 
exp1_filename = searchdir(exp1_dir, "full_OHC_"*basin_name)[1]
exp1_path = joinpath(exp1_dir, exp1_filename)

println("Comparing the two datasets")

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
E,F = trend_matrices(tecco)
perc_diff = Dict()
for (keys,values) in shortnames
    k +=1
    if values ∉ skip_exp
        slopes = zeros(Float64, length(z))
        intercepts = zeros(Float64, length(z))
        errs = zeros(Float64, 2, length(z),)
        uplvl = -2000; botlvl = -3000
        # lvl1 = -4000; lvl2 = -5000
        lvls = findall( botlvl .<= z[:].<= uplvl)
        y_true = sum(Real.(exps[keys][lvls, :] ), dims = 1) .* 1e-21
        β = F * y_true' 
        ts = tecco .- mean(tecco)
        y_pred =  β[1] .+ (β[2] .* ts)
        num = sum((y_true .- y_pred).^2)
        denom = sum(ts.^2)
        SE = sqrt(num * inv(length(ts) - 2) * inv(denom))
        t = β[2]/SE
        println(SE * 26)
        ax.errorbar(values, β[2]* 26, yerr = SE * 26 , fmt = "o")
        perc_diff[values] = β[2]* 26
        
    end
end
# fig.subplots_adjust(top=0.8, left=0.1, right=0.9, s=0.12)
# legend_size = Dict("size" => 12)
# ax.legend(loc="center",bbox_to_anchor=(0.5, -0.3), ncol=4, prop=legend_size)
# fig.suptitle(basin_name * " Ocean Heat Content Change [1992-2018]", size = 15)
# fig.tight_layout()
perc_change_dict = Dict()

for (keys,values) in shortnames
    if values ∉ skip_exp
        baseline = perc_diff["iter0"]
        perc_diff[values] = (perc_diff[values]) / baseline
    end
end

for (keys,values) in shortnames
    if values ∉ skip_exp
        baseline = perc_diff["iter0"]
        perc_change_dict[values] = rel_perc_diff(perc_diff[values], baseline)
    end
end

println("fraction change from iter0");println(perc_diff)
println("perc change from iter0");println(perc_change_dict)

ax.set_title("Estimated ΔH_{2-3km} over 26 years")
outputfile = plotsdir("OHCstats/ΔOHC_2to3km"*basin_name * ".pdf")
fig.savefig(outputfile)
# close("all")