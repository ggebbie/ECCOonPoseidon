using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyCall


# include(srcdir("config_exp.jl"))

fig, axs = plt.subplots(1, 4,  sharex="row", figsize =(27, 14))
# fig.suptitle("Total Δθ̄ by each budget term ")
budget_colors = sns.color_palette("Paired")[7:11]
budget_colors[end] = (0, 0, 0)
extr = [0, 0]
labels_L = [L"\theta^{129}", L"\theta^{\Delta F}",  L"\theta^{\Delta T}", L"\theta^{0}"]
budg_terms = ["AdvH", "AdvR", "DiffH", "DiffR", "total"]
baseline = θ_budget_deltas["iter129_bulkformula"]


Dict(key1 => Dict(key2 => d[key1][key2] .- d["2"][key2] for 
key2 in keys(d[key1])) for key1 in keys(d))

for (i, key) in enumerate(keys(shortnames))
    expname = key
    axs[i].set_title(labels_L[i])
    j = 0 
    for contr in budg_terms
        j+=1
        axs[i].plot(θ_budget_deltas[key][contr], z[lvls], color = budget_colors[j], 
        label = contr)
    end
    (i < 2) && (    axs[i].set_ylabel(L"depth[m]")    )
    (i > 1) && (axs[i].tick_params(labelleft=false))
    axs[i].set_xlabel("Net " * L"^\circ" * "C")

    # axs[i].set_xlim([-0.8, 0.8])
end
# axs[4].legend(bbox_to_anchor=(-1.15, -0.15), loc="upper center", ncol = 5,borderaxespad= 0.2)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Depth_Plots/" * "DepthθBudgChange_" * region * suffix * ".png")
close("all")

#Now see if there is any correlation 
#(plot surface warming as a trends as a function of deep cooling)