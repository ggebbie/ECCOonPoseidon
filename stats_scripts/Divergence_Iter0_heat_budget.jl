
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
sns.set(); pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;


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

uplvl = -2000; botlvl = -3000
tt3km = findall( botlvl .<= z[:].<= uplvl)
region = "NPAC"; 
ocean_mask = wet_pts(Γ)
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ; region)/tmp/recove
msk = PAC_msk;

cell_depths = get_cell_depths(msk, ΔzF, Γ.hFacC)
cell_volumes = get_cell_volumes(area, cell_depths);
""" get the bottom heating """
lvls = tt3km
bottom_lev = cell_volumes[:, lvls[end] + 1]
has_bottom = similar(bottom_lev)
has_bottom[findall(bottom_lev .> 0) ] .= 0.0
has_bottom[findall(bottom_lev .== 0) ] .= inv(1029 * 3994)
has_bottom = has_bottom .* msk
#geothermal fluxes are not adjusted
GTF = has_bottom .* get_geothermalheating(γ)
"""Compare vertical convergence of pot. temperature"""
filedir = "ECCO_vars/"
filename1 = "DeepVertHeatConv";
filename2 = "DeepTHETAZ"
filename3 = "DeepLateralADV"
filename4 = "DeepLateralDF"

θC, θz = Dict(), Dict()
lvls = tt3km
for (keys,values) in shortnames
    θC[keys] = Float64[]
    θz[keys] = Float64[]
    expname = keys; println(keys)
    @load datadir(filedir*filename2*"_"*expname*".jld2") var_exp
    LHS = copy(var_exp); var_exp = nothing
    @load datadir(filedir*filename1*"_"*expname*".jld2") var_exp
    VertF = copy(var_exp); var_exp = nothing
    @load datadir(filedir*filename3*"_"*expname*".jld2") var_exp
    HorADV = copy(var_exp); var_exp = nothing
    @load datadir(filedir*filename4*"_"*expname*".jld2") var_exp
    HorDIF= copy(var_exp); var_exp = nothing

    for tt in 1:length(tecco)
        push!(θz[keys], sum(LHS[tt] .* msk)) #total tendency
        tmp = -(VertF[tt] ./ area)
        tmp = tmp .- GTF
        tmp = (tmp .+ HorADV[tt+312]) .+ HorDIF[tt+312]
        push!(θC[keys], sum(tmp) )
    end
end

dθz_LHS,dθz_RHS, resids, means = Dict(), Dict(), Dict(), Dict()
for (keys,values) in shortnames
    dθz_LHS[keys] =  (θz[keys][3:end] - θz[keys][1:end-2]) ./ (2 * 2.62e6)
    dθz_RHS[keys] =  θC[keys][2:end-1] 
    resids[keys] = dθz_LHS[keys] .- dθz_RHS[keys]
    means[keys * "_LHS"] = mean(dθz_LHS[keys]);
    means[keys * "_RHS"] = mean(dθz_RHS[keys])
end

# tecco_off = tecco[1:end-1] .+  (diff(tecco) ./2)
tecco_off = tecco[2:end-1]
fig, ax = subplots(2, 1, figsize = (12, 7))
ax[1].set_title("Divergence From 0bf (" * region * ", z=2-3km)")
ax[2].set_title(" Total " * region * " Vert. Heat Flux (LHS), z=2-3km")
plot_div_0bf!(dθz_LHS, tecco_off, shortnames, ignore_list, ax[1]; ylabel =  L"Δ ^\circ C m / s")
plot_ts!(dθz_LHS, tecco_off, shortnames, ignore_list, ax[2]; ylabel =  L" ^\circ C m / s", baseline = 0)
tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "TotθVFlxLHS_" * region * "_2km3km.png")

fig, ax = plt.subplots(2, 1, figsize = (12, 7))
ax[1].set_title("Divergence From 0bf (" * region * ", z=2-3km)")
ax[2].set_title(" Total " * region * " Vert. Heat Flux (RHS), z=2-3km")
plot_div_0bf!(dθz_RHS, tecco_off, shortnames, ignore_list, ax[1]; ylabel =  L"Δ ^\circ C m / s")
plot_ts!(dθz_RHS, tecco_off, shortnames, ignore_list, ax[2]; ylabel =  L" ^\circ C m / s", baseline = 0)
tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "TotθVFlxRHS_" * region * "_2km3km.png")

fig, ax = plt.subplots(1, figsize = (12, 7))
ax.set_title( region * " Vert. Heat Flux Residual (LHS - RHS), z=2-3km")
plot_ts!(resids, tecco_off, shortnames, ignore_list, ax; ylabel =  L" ^\circ C m / s", baseline =  0)
# tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "TotθVFlxResid_" * region * "_2km3km.png")

θbar = Dict()
θbarLHS, θbarRHS = Dict(),Dict()

for (keys,values) in shortnames
    TD = sum(cell_depths[:, lvls])
    θbar[keys] = θz[keys] ./ TD
    start = θz[keys][2] / TD
    θbarLHS[keys] =  cumsum(vcat(start, dθz_LHS[keys].* (2.62e6 / TD)))
    θbarRHS[keys] =  cumsum(vcat(start, dθz_RHS[keys].* (2.62e6 / TD)))
end

fig, ax = plt.subplots(1, 2, figsize = (12, 5))
ax[1].set_title( region *  " θ̄ Reconstruction (LHS),  z=2-3km")
ax[2].set_title( region * " θ̄ Reconstruction (RHS),  z=2-3km")
# plot_ts!(θbar, tecco, shortnames, ignore_list, ax[1]; ylabel =  L" ^\circ C")
plot_ts!(θbarLHS, tecco[2:end], shortnames, ignore_list, ax[1]; ylabel =  L" ^\circ C")
plot_ts!(θbarRHS, tecco[2:end], shortnames, ignore_list, ax[2]; ylabel =  L" ^\circ C")
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θReconstruct_" * region * "_2km3km.png")

close("all")

