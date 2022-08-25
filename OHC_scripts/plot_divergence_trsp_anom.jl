using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyCall
using DataFrames, GLM

ocean_mask = wet_pts(Γ)

filedir = "ECCO_vars/"
trspfname = "NTRSP"
nt = length(tecco)

lat_mask = MeshArray(γ, Int64)
ϕround = round.(ϕ, digits = 2)
lat_mask[findall(ϕround .== 30.11) ] .= 1;
lat_mask[findall(ϕround .!= 30.11) ] .= 0;
lat_mask[findall(ocean_mask .== 0)] .= 0;

meridion_trsp = Dict()
for expname in keys(shortnames)
    meridion_trsp[expname] = Float64[]
    @time for tt in 1:nt
        fdir = filedir * expname * "/" * trspfname * "_full_"*string(tt) *".jld2"
        NTRSP = load_object(datadir(fdir)) .* lat_mask
        push!(meridion_trsp[expname], sum(NTRSP) * 1e-6)
    end
    @time GC.gc() #garbage collecting 
end


labels_L = [L"\theta^{129}", L"\theta^{\Delta F}",L"\theta^{\Delta T}", L"\theta^{0}"]

fig, ax =  plt.subplots(1, 1, figsize =(15, 15))
for (i, key) in enumerate(keys(shortnames))
    expname = key
    baseline = meridion_trsp["iter0_bulkformula"]
    ax.plot(tecco, meridion_trsp[expname] .- baseline, color = colors[i], label = labels_L[i])
    ax.set_xlabel("time (years)")
    ax.set_ylabel("Sv")
end
ax.set_title("Volume Flux Anomaly through 30" * L"^\circ" * "N \n" *suffix)
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=4)
fig.tight_layout();
fig.savefig(plotsdir() * "/OHC_Divergence/MassTransport/" * "VolFlxAnom" * "_" * suffix * ".png")

fig, ax =  plt.subplots(1, 1, figsize =(15, 15))
for (i, key) in enumerate(keys(shortnames))
    expname = key
    baseline = meridion_trsp["iter0_bulkformula"]
    ax.plot(tecco, meridion_trsp[expname], color = colors[i], label = labels_L[i])
    ax.set_xlabel("time (years)")
    ax.set_ylabel("Sv")
end
ax.set_title("Volume Flux through 30" * L"^\circ" * "N \n" *suffix)
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=4)
fig.tight_layout();

# fig.tight_layout();
fig.savefig(plotsdir() * "/OHC_Divergence/MassTransport/" * "VolFlx" * "_" * suffix * ".png")

