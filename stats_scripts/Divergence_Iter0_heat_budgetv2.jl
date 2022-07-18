
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
fname1 = "Deep_THETA"
fname2 = "THETA_BUDG";

θC, θCAR, θCAH, θCDH, θCDR, θz = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
budg_names = ["AdvR", "AdvH", "DiffH", "DiffR"]
lvls = tt3km
tot_vol = sum(cell_volumes[:, lvls])
GTF_tot = sum(GTF) ./ tot_vol
vol_weight(x) = sum(x .* msk) / tot_vol
for (key,values) in shortnames
    expname = key; println(key)
    @time @load datadir(filedir*fname1*"_"*expname*".jld2") var_exp
    θ = var_exp
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
        
        avg1 = sum(tot.* msk ) / tot_vol
        avg2 = mean(θ[tt], cell_volumes[:, lvls])

        push!(θC[expname], avg1)
        push!(θz[expname], avg2  )
    end
    @time GC.gc(true) #garbage collecting 
end
dθz_LHS,dθz_RHS, resids, means = Dict(), Dict(), Dict(), Dict()
central_diff(x, Δt) = (x[3:end] .- x[1:end-2]) ./ (2*Δt)
perc_diff(x, y) = (x - y)/y
resid(x, y) = x - y
metric_name = "Residual"
metric_units = L"Δ ^\circ C / s"
for (keys,values) in shortnames
    dθz_LHS[keys] =  central_diff(θz[keys], 2.62e6)
    dθz_RHS[keys] =  θC[keys][2:end-1]
    resids[keys] = resid.(dθz_LHS[keys],dθz_RHS[keys])
end

# tecco_off = tecco[1:end-1] .+  (diff(tecco) ./2)
tecco_off = tecco[2:end-1]
fig, ax = subplots(2, 1, figsize = (12, 7))
ax[1].set_title("Divergence From 0bf (" * region * ", z=2-3km)")
ax[2].set_title(" Total " * region * " Heat Flux (LHS), z=2-3km")
plot_div_0bf!(dθz_LHS, tecco_off, shortnames, ignore_list, ax[1]; ylabel =  L"Δ ^\circ C / s")
plot_ts!(dθz_LHS, tecco_off, shortnames, ignore_list, ax[2]; ylabel =  L" ^\circ C / s", baseline = 0)
tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "TotθFlxLHS_" * region * "_2km3km.png")

fig, ax = subplots(2, 1, figsize = (12, 7))
ax[1].set_title("Divergence From 0bf (" * region * ", z=2-3km)")
ax[2].set_title(" Total " * region * " Heat Flux (RHS), z=2-3km")
plot_div_0bf!(dθz_RHS, tecco_off, shortnames, ignore_list, ax[1]; ylabel =  L"Δ ^\circ C / s")
plot_ts!(dθz_RHS, tecco_off, shortnames, ignore_list, ax[2]; ylabel =  L" ^\circ C / s", baseline = 0)
tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "TotθFlxRHS_" * region * "_2km3km.png")

fig, ax = plt.subplots(1, figsize = (12, 7))
ax.set_title( region * " Heat Flux " * metric_name *", z=2-3km")
plot_ts!(resids, tecco_off, shortnames, ignore_list, ax; ylabel =  metric_units, baseline =  0)
tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "TotθFlxResid_" * region * "_2km3km.png")

fig1, axs = subplots(4, 1, figsize = (7, 15))
fig1.suptitle(region * " Heat Flux " * metric_name *" , z=2-3km")

for (keys,values) in shortnames, ax in axs
    ax.hist(resids[keys], bins = 10)
    ax.set_ylabel(L"Frequency")
    ax.set_xlabel(L"Frequency")
    ax.set_title(keys)
end
tight_layout(); fig1.savefig(plotsdir() * "/OHC_Divergence/" * "TotθFlxResidHis_" * region * "_2km3km.png")

fig2, ax2 = subplots(1, 1, figsize = (10, 5))
ax2.boxplot(collect(values(resids)), showmeans=true)
ax2.set_xticklabels(keys(resids))
ax2.set_title(region * " Heat Flux Budget Residual (%), z=2-3km")
ax2.set_ylabel(L"\%")
tight_layout(); fig2.savefig(plotsdir() * "/OHC_Divergence/" * "TotθFlxResidBox_" * region * "_2km3km.png")

fig, ax = subplots(4, 1, figsize = (7, 15))
for (keys,values) in shortnames
    j = 0
    for var in [θCAR, θCAH, θCDH, θCDR]
        j+= 1
        ax[i].hist(var[keys], label=budg_names[j])
    end
    ax[i].set_title("Distribution of Monthly Heat Budget Terms for " 
    * values * "in " * region * ", z = 2-3km")
    ax[i].set_xlabel(L"\circ C / s")
    ax[i].legend()
end
tight_layout(); savefig(plotsdir() * "/OHC_Divergence/" * "HeatBudgDist" * region * "_2km3km.png")

θbar = Dict()
θbarLHS, θbarRHS = Dict(),Dict()
reconstruct(x0, dx) = cumsum(vcat(x0, dx .* (2.62e6 )))

for (keys,values) in shortnames
    start = θz[keys][2]
    θbarLHS[keys] =  reconstruct(start, dθz_LHS[keys]) 
    θbarRHS[keys] =  reconstruct(start, dθz_RHS[keys])
end
fig, ax = subplots(2, 3, figsize = (12, 7))

ax[1].set_title(region *  " θ̄,  z=2-3km")
# ax[2].set_title( region * " θ̄ Reconstruction,  z=2-3km")
# plot_ts!(θbar, tecco, shortnames, ignore_list, ax[1]; ylabel =  L" ^\circ C")
plot_ts!(θbarLHS, tecco[2:end], shortnames, ignore_list, ax[1]; ylabel =  L" ^\circ C")
plot_ts!(θbarRHS, tecco[2:end], shortnames, ignore_list, ax[end]; ylabel =  L" ^\circ C")
ax[end].set_title("Sum of all terms")

i = 0
z = Dict()
for var in [θCAR, θCAH, θCDH, θCDR]
    i+=1
    temp_dict = Dict()
    for (keys,values) in shortnames 
        start = θz[keys][2]
        temp_dict[keys] =  reconstruct(start, var[keys][2:end-1])
    end
    z[budg_names[i]] = temp_dict["iter129_bulkformula"]
    plot_ts!(temp_dict, tecco[2:end], shortnames, ignore_list, ax[i+1]; ylabel =  L" ^\circ C")
    ax[i+1].set_title(budg_names[i])
end

fig.suptitle(region * " θ̄ Reconstructed from \n heat budget terms")
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θReconstructAll_" * region * "_2km3km.png")

fig, ax = subplots(1, 1, figsize = (5, 5))
start = 0.0 
GTF_temp = reconstruct(start, fill(GTF_tot, length(tecco[2:end-1])))
ax.plot(tecco[2:end], GTF_temp, color="black")
ax.set_ylabel(L" ^\circ C")
ax.set_title("Geothermal Heat Flux")
fig.savefig(plotsdir() * "/OHC_Divergence/" * "GeothermalFlx" * region * "_2km3km.png")
close("all")


