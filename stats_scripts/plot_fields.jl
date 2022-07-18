
# depthlbl = string(abs(round(z[lvl_idx], digits = 2)))


# this script will compile and save a basinwide average of OHC 
# for a particular basin (basin_name) 

include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf
using .OHC_helper
using PyPlot   # important!
using PyCall
@pyimport seaborn as sns
sns.set()
pygui(false)
cm = pyimport("cmocean.cm")
colorway = cm.balance;

include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
tecco = collect(Float64, 1992+1/24:1/12:2018)

runpath,diagpath = listexperiments(exprootdir())
shortnames = expnames()
ignore_list = ["noIA", "129ff"]
shortnames = reduce_dict(shortnames, ignore_list)

marks = expsymbols()
nexp = length(shortnames) # number of experiments

ocean_mask = wet_pts(Γ)
msk = Float32.(ocean_mask)
cell_depths = get_cell_depths(msk, ΔzF, Γ.hFacC)

uplvl = -2000; botlvl = -3000
tt3km = findall( botlvl .<= z[:].<= uplvl)

filedir = "ECCO_vars/"
filename = "THETA_AFLX";
expname = "iter0_bulkformula"
@load datadir(filedir*filename*"_"*expname*".jld2") var_exp
iter0_init = copy(var_exp[1][:, 50]); var_exp = nothing 
# pre-allocate β, linear trends

proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
fig, ax = plt.subplots(1, 1, figsize=(6,5), 
                        subplot_kw=Dict("projection"=> proj0))
cf = plot_field!(iter0_init, λ, ϕ, ax, colorway)
ax.set_extent((110, -70, -55, 60))
# lvls_str = string(abs(uplvl)) * "to" * string(abs(botlvl))
ax.set_title("Level 50 Vertical Advective Temp. Flux 0bf")
# depthlbl = string(abs(round(z[lvl_idx], digits = 2)))
cbar = fig.colorbar(cf, label = L"^\circ"*"C")

#Saving the Figureå
folder = "/home/ameza/ECCOonPoseidon/plots/OHC Climatologies/" * expname * "/"
mkpath(folder[1:end-1])
# lvls_str = "lvl" * string(lvl_idx)
outputfile = plotsdir(folder * "vertAFlux" * expname * ".pdf")

savefig(outputfile)
println("saving plots at.. " * outputfile)

close("all")