using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyCall
using DataFrames, GLM

filedir = "ECCO_vars/"
trspfname = "NTRSP"
nt = length(tecco)

lat_mask = MeshArray(γ, Int64)
ϕround = round.(ϕ, digits = 2)
lat_mask[findall(ϕround .== 30.11) ] .= 1;
lat_mask[findall(ϕround .!= 30.11) ] .= 0;
lat_mask[findall(msk .== 0)] .= 0;

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

fig, ax =  plt.subplots(1, 1, figsize =(15, 7.5))
for (i, key) in enumerate(keys(shortnames))
    expname = key
    ax.plot(tecco, meridion_trsp[expname], color = colors[i], label = key)
    ax.set_xlabel("time (years)")
    ax.set_ylabel("Sv")
end
ax.set_title("Volume Flux Anomaly \n through 30" * L"^\circ" * "N \n" *suffix)
ax.legend()
# fig.tight_layout();
fig.savefig(plotsdir() * "/OHC_Divergence/MassTransport/" * "VolFlx" * "_" * region *  "_" * suffix * ".png")
