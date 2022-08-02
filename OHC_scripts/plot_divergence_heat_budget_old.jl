using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf
using PyPlot   # important!
using PyCall
@pyimport seaborn as sns
sns.set(); pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;

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
