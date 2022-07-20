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
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ; region)
msk = PAC_msk;

cell_depths = get_cell_depths(msk, ΔzF, Γ.hFacC)
cell_volumes = get_cell_volumes(area, cell_depths);
""" get the bottom heating """
lvls = tt3km
bottom_lev = cell_volumes[:, lvls[end] + 1]
has_bottom = similar(bottom_lev)
has_bottom[findall(bottom_lev .> 0) ] .= 0.0
has_bottom[findall(bottom_lev .== 0) ] .= 1.0
has_bottom = has_bottom .* msk
#geothermal fluxes are not adjusted
GTF = has_bottom .* get_geothermalheating(γ)
"""Compare vertical convergence of pot. temperature"""
filedir = "ECCO_vars/"
# fname1 = "Deep_THETA"
fname1 = "Deep_STHETA"
fname2 = "THETA_BUDG";

θC, θCAR, θCAH, θCDH, θCDR, θz = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
budg_names = ["AdvR", "AdvH", "DiffH", "DiffR"]
lvls = tt3km
tot_vol = Float32(sum(cell_volumes[:, lvls]))
vol_weight(x) = Float32(sum(x .* msk) / tot_vol)
for (key,values) in shortnames
    expname = key; println(key)
    @time @load datadir(filedir*fname1*"_"*expname*".jld2") var_exp
    sθ = var_exp
    @time @load datadir(filedir*fname2*"_"*expname*".jld2") var_exp
    HBUDG = var_exp; 
    for var in [θz, θC, θCAH, θCAR, θCDH, θCDR]
        var[expname] = Float64[]
    end
    for tt in 1:length(tecco)
        tot = @views HBUDG["AdvH"][tt] .+ HBUDG["DiffH"][tt] .- 
                    (HBUDG["AdvR"][tt] .+ HBUDG["DiffZ"][tt] .- GTF )

        AdvH =  vol_weight(HBUDG["AdvH"][tt]);push!(θCAH[expname],AdvH)
        AdvR =  vol_weight(-HBUDG["AdvR"][tt]);push!(θCAR[expname],AdvR)
        DiffH = vol_weight(HBUDG["DiffH"][tt]);push!(θCDH[expname],DiffH)
        DiffR = vol_weight(-HBUDG["DiffZ"][tt]);push!(θCDR[expname],DiffR)
        
        avg1 = vol_weight(tot)
        avg2 = mean(sθ[tt]; weights = cell_volumes[:, lvls])

        push!(θC[expname], Float32.(avg1))
        push!(θz[expname], Float32.(avg2))
    end
    @time GC.gc(true) #garbage collecting 
end

dθz_LHS,dθz_RHS, resids, means = Dict(), Dict(), Dict(), Dict()
Δt₁ = Float32.((tecco[3:end] .- tecco[1:end-2]) .* (12 * 30.5 * 86400))
central_diff(x, dt) = (x[3:end] .- x[1:end-2]) ./ (dt)
Δt₂ = mean(Float32.((tecco[2:end] .- tecco[1:end-1]) .* (12 * 30.5 * 86400)))
fwd_diff(x, dt) = (x[2:end] .- x[1:end-1]) ./ (dt)
fwd_mean(x)=(x[2:end].+x[1:end-1])./2
perc_diff(x, y) = 100 * (x - y)/y
resid(x, y) = x - y
fcycle = 1 # units: yr^{-1}
# for removing seasonal cycle from monthly averaged timeseries
overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco[2:end],overtones)

metric_name = "PercDiffPruned"
# metric_units = L"Δ ^\circ C / s"
metric_units = "% Difference"
remove_seasonal_with_mean(x, Ecycle, Fcycle) = remove_seasonal(x,Ecycle,Fcycle) .+ mean(x)
for (keys,values) in shortnames
    dθz_LHS[keys] =  Float32.(fwd_diff(θz[keys], Δt₂))
    dθz_RHS[keys] =  Float32.(fwd_mean(θC[keys]))
    dθz_LHS[keys] = remove_seasonal_with_mean(dθz_LHS[keys],Ecycle,Fcycle)
    dθz_RHS[keys] = remove_seasonal_with_mean(dθz_RHS[keys],Ecycle,Fcycle)
    # resids[keys] = resid.(dθz_LHS[keys],dθz_RHS[keys])
    resids[keys] = prune(perc_diff.(dθz_RHS[keys],dθz_LHS[keys]))
    σ=std(prune(resid.(dθz_LHS[keys],dθz_RHS[keys])))/std(dθz_LHS[keys])
    println(values," ratio of STD ",string(σ))
end

tecco_off = tecco[1:end-1] .+  (diff(tecco) ./2)
# tecco_off = tecco[2:end-1]
fig, ax = subplots(2, 2, figsize = (24, 8))
ax[2].set_title("Divergence From 0bf (" * region * ", z=2-3km)")
ax[1].set_title(" Total " * region * " Heat Flux (LHS), z=2-3km")
ax[4].set_title("Divergence From 0bf (" * region * ", z=2-3km)")
ax[3].set_title(" Total " * region * " Heat Flux (RHS), z=2-3km")
plot_div_0bf!(dθz_LHS, tecco_off, shortnames, ignore_list, ax[2]; ylabel =  L"Δ ^\circ C / s")
plot_ts!(dθz_LHS, tecco_off, shortnames, ignore_list, ax[1]; ylabel =  L" ^\circ C / s", baseline = 0)
plot_div_0bf!(dθz_RHS, tecco_off, shortnames, ignore_list, ax[4]; ylabel =  L"Δ ^\circ C / s")
plot_ts!(dθz_RHS, tecco_off, shortnames, ignore_list, ax[3]; ylabel =  L" ^\circ C / s", baseline = 0)
tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "TotθFlx_" * region * "_2km3km.png")

fig, ax = plt.subplots(1, figsize = (12, 7))
ax.set_title( region * " Heat Flux " * metric_name *", z=2-3km")
plot_ts!(resids, tecco_off, shortnames, ignore_list, ax; ylabel =  metric_units, baseline =  0)
tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "TotθFlx" * metric_name * "_" * region * "_2km3km.png")

slice(i) =pycall(pybuiltin("slice"), PyObject, i)
slice(i, j) =pycall(pybuiltin("slice"), PyObject, i, j)
fig = plt.figure( figsize = (12, 7))
gs = fig.add_gridspec(4,5)
element(i,j) = get(gs, (i,j))
a = fig.add_subplot(element(slice(0,4), slice(2,5)))
a.boxplot(collect(values(resids)), showmeans=true)
a.set_xticklabels(keys(resids))
a.set_title(region * " Heat Flux Budget " * metric_name * ", z=2-3km")
a.set_ylabel(metric_units)

axs = [fig.add_subplot(element(0,slice(0,2))), fig.add_subplot(element(1,slice(0,2))), 
fig.add_subplot(element(2,slice(0,2))), fig.add_subplot(element(3,slice(0,2)))]
i = 0 
for (keys,values) in shortnames
    i+=1
    axs[i].hist(resids[keys], bins = 25, alpha = 0.3)
    axs[i].set_ylabel(L"Frequency")
    axs[i].set_xlabel(metric_units)
    axs[i].set_title(keys)
end
tight_layout(); fig.savefig(plotsdir() * "/OHC_Divergence/" * "TotθFlx" * metric_name * "Stats_" * region * "_2km3km.png")

fig, ax = subplots(4, 1, figsize = (10, 12))
i = 0 
for (keys,values) in shortnames
    i+=1
    j = 0
    for var in [θCAR, θCAH, θCDH, θCDR]
        j+= 1
        ax[i].hist(var[keys], label=budg_names[j], alpha = 0.3)
    end
    ax[i].set_title(values)
    ax[i].set_xlabel(L"\circ C / s")
    ax[i].legend()
end
fig.suptitle("Distribution of Monthly Heat Budget Terms in " * region * ", z = 2-3km")
tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "HeatBudgDist" * region * "_2km3km.png")

θbar = Dict()
θbarLHS, θbarRHS = Dict(),Dict()
Δt = mean(Float32.((tecco[2:end] .- tecco[1:end-1]) .* (12 * 30 * 86400)))
reconstruct(x0, dx) = cumsum(vcat(x0, dx .* Δt))

for (keys,values) in shortnames
    start = θz[keys][2]
    θbarLHS[keys] =  reconstruct(start, dθz_LHS[keys])
    θbarRHS[keys] =  reconstruct(start, dθz_RHS[keys])
end
fig, ax = subplots(2, 3, figsize = (12, 9))
ax[1].set_title(region *  " θ̄,  z=2-3km")
plot_ts!(θbarLHS, tecco, shortnames, ignore_list, ax[1]; ylabel =  L" ^\circ C")
plot_ts!(θbarRHS, tecco, shortnames, ignore_list, ax[end]; ylabel =  L" ^\circ C")
ax[end].set_title("Sum of all terms")

i = 0
z = Dict()
for var in [θCAR, θCAH, θCDH, θCDR]
    i+=1
    temp_dict = Dict()
    for (keys,values) in shortnames 
        start = θz[keys][2]
        temp_dict[keys] =  reconstruct(start/4, var[keys][2:end])
    end
    z[budg_names[i]] = temp_dict["iter129_bulkformula"]
    plot_ts!(temp_dict, tecco, shortnames, ignore_list, ax[i+1]; ylabel =  L" ^\circ C")
    ax[i+1].set_title(budg_names[i])
end

fig.suptitle(region * " θ̄ Reconstructed from \n heat budget terms")
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θReconstructAll_" * region * "_2km3km.png")

fig, ax = subplots(1, 1, figsize = (5, 5))
start = 0.0 
GTF_tot = sum(GTF) ./ tot_vol
GTF_temp = reconstruct(start, fill(GTF_tot, length(tecco[2:end])))
ax.plot(tecco, GTF_temp, color="black")
ax.set_ylabel(L" ^\circ C")
ax.set_title("Geothermal Heat Flux")
fig.savefig(plotsdir() * "/OHC_Divergence/" * "GeothermalFlx" * region * "_2km3km.png")
close("all")
