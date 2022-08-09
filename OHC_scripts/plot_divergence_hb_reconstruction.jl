using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
cm = pyimport("cmocean.cm");colorway = cm.balance;
mpatches = pyimport("matplotlib.patches")
lines = pyimport("matplotlib.lines")
Line2D = lines.Line2D

θbar = Dict()
θbarLHS, θbarRHS = Dict(),Dict()
reconstruct(x0, dx) = cumsum(vcat(Float32(x0), dx .* Float32(Δt₂)))

for key in keys(shortnames)
    start = Float32(θz[key][1])
    θbarLHS[key] =  Float32.(reconstruct(start, dθz_LHS[key]))
    θbarRHS[key] =  Float32.(reconstruct(start, dθz_RHS[key]))
end

fig, ax = subplots(2, 3, figsize = (16, 14))
fig.suptitle(region * " θ̄ Reconstructed from \n heat budget terms")
start = 0.0 
GTF_tot = GTF
GTF_temp = reconstruct(start, fill(GTF_tot, length(tecco[2:end])))
ax[1].plot(tecco, GTF_temp, color="black")
ax[1].set_ylabel(L" ^\circ C")
ax[1].set_title("Geothermal Heat Flux")

for (i, var) in enumerate([θCAR, θCAH, θCDH, θCDR])
    temp_dict = Dict()
    for (keys,values) in shortnames 
        start = θz[keys][1]
        temp_dict[keys] =  reconstruct(start/4, var[keys][2:end])
    end
    ax[i+1].set_title(budg_names[i])
    plot_ts!(temp_dict, tecco, shortnames, ignore_list, ax[i+1]; ylabel =  L" ^\circ C")
end

ax[end].set_title(region *  " θ̄,  z=2-3km")
plot_ts!(θbarLHS, tecco, shortnames, ignore_list, ax[end]; ylabel =  L" ^\circ C", linestyle = "--")
# ax[end].set_title("Sum of all terms")
plot_ts!(θbarRHS, tecco, shortnames, ignore_list, ax[end]; ylabel =  L" ^\circ C")
map(x -> x.get_legend().remove(), ax[2:6])
# ax[4].get_legend().remove()
patches = [Line2D([0], [0], color=colors[i], label=collect(values(shortnames))[i]) for i in 1:4]
push!(patches, Line2D([0], [0],color="k", label="Reconstruction"))
push!(patches, Line2D([0], [0],color="k", label="Approx θ", linestyle = "--"))
ax[4].legend(handles = patches)
sns.move_legend(ax[4], "lower center", bbox_to_anchor=(.5, -.5), ncol=2)
fig.tight_layout()

fig.savefig(plotsdir() * "/OHC_Divergence/" * "θReconstructAllTerms_" * region * suffix * "_2km3km.png")
close("all")

fig, ax = subplots(1, 1, figsize = (12, 9))
ax.set_title(region *  " θ̄,  z=2-3km")
plot_ts!(θz, tecco, shortnames, ignore_list, ax; ylabel =  L" ^\circ C", linestyle = "--")
plot_ts!(θbarRHS, tecco, shortnames, ignore_list, ax; ylabel =  L" ^\circ C")
ax.get_legend().remove()
patches = [Line2D([0], [0], color=colors[i], label=collect(keys(shortnames))[i]) for i in 1:4]
push!(patches, Line2D([0], [0],color="k", label="Reconstruction"))
push!(patches, Line2D([0], [0],color="k", label="Actual", linestyle = "--"))
ax.legend(handles = patches)

fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θReconstruction_" * region * suffix * "_2km3km.png")
close("all")


fig, ax = subplots(1, 1, figsize = (12, 9))
ax.set_title(region *  " θ̄,  z=2-3km")
plot_ts!(θbarLHS, tecco, shortnames, ignore_list, ax; ylabel =  L" ^\circ C", linestyle = "--")
# ax[end].set_title("Sum of all terms")
plot_ts!(θbarRHS, tecco, shortnames, ignore_list, ax; ylabel =  L" ^\circ C")
ax.get_legend().remove()
patches = [Line2D([0], [0], color=colors[i], label=collect(keys(shortnames))[i]) for i in 1:4]
push!(patches, Line2D([0], [0],color="k", label="Reconstruction"))
push!(patches, Line2D([0], [0],color="k", label="Approx θ", linestyle = "--"))
ax.legend(handles = patches)

fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θReconstructionApprox_" * region * suffix * "_2km3km.png")
close("all")
####EXTRA STUFF 
t=Float32.(θz["iter129_bulkformula"][2:end].-θz["iter129_bulkformula"][1:end-1])/(Δt₂)
test=cumsum(vcat(θz["iter129_bulkformula"][1], t .* Δt₂))

Float32(θz["iter129_bulkformula"][224])+(t[224]*Float32(Δt₂))
θz["iter129_bulkformula"][225]

Gtf = -(θbarRHS["iter129_bulkformula"][end] .- θbarLHS["iter129_bulkformula"][end])
stoyr = 3.154e+7 / 1
len_in_sec = (tecco[end] - tecco[1]) * stoyr
x = Gtf/len_in_sec

a = θbarRHS["iter129_bulkformula"] .- θbarLHS["iter129_bulkformula"]