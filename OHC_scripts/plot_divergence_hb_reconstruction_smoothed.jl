using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyCall

cm = pyimport("cmocean.cm");colorway = cm.balance;
mpatches = pyimport("matplotlib.patches")
lines = pyimport("matplotlib.lines")
Line2D = lines.Line2D

θbar = Dict()
θbarLHS, θbarRHS = Dict(),Dict()
Δt = Δt₁

for key in keys(shortnames)
    start = θz[key][7]
    θbarLHS[key] =  reconstruct(start, dθz_LHS_Annual[key])
    θbarRHS[key] =  reconstruct(start, dθz_RHS_Annual[key])
end

fig, ax = subplots(1, 1, figsize = (12, 9))
ax.set_title(region *  " θ̄,  z=2-3km")
idx=round(Int,window_size/2)
plot_ts!(θbarLHS, tecco[idx+1:end-idx], shortnames, ignore_list, ax; ylabel =  L" ^\circ C", linestyle = "--")
# ax[end].set_title("Sum of all terms")
plot_ts!(θbarRHS, tecco[idx+1:end-idx], shortnames, ignore_list, ax; ylabel =  L" ^\circ C")
ax.get_legend().remove()
patches = [Line2D([0], [0], color=colors[i], label=collect(keys(shortnames))[i]) for i in 1:4]
push!(patches, Line2D([0], [0],color="k", label="Reconstruction"))
push!(patches, Line2D([0], [0],color="k", label="Approx θ", linestyle = "--"))
ax.legend(handles = patches)

fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θReconstructionApproxAnnual_" * region * suffix *  "_2km3km.png")
close("all")