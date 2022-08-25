
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions, NonNegLeastSquares
using PyCall
using DataFrames, GLM
#df = DataFrame(θz); df = df .- df["iter0_bulkformula] #initiate dataframe & do mean about iter0
#
# θsurf_noszn = ref_month_1( θsurf_noszn,"iter0_bulkformula")
# θdeep_noszn = ref_month_1( θz,"iter0_bulkformula")
reference_text = ""
# df = DataFrame(θdeep_noszn)
df = DataFrame(θz);
fig, ax = plt.subplots(1, figsize = (15, 15))
fig.suptitle("Average N. Pac. Potential Temperature \n (2000 - 3000 m)")
i = 0
ax.plot(tecco, df.iter129_bulkformula, c = colors[1], label = L"\theta^{129}")
ax.plot(tecco, df.noinitadjust,  c = colors[2], label = L"\theta^{\Delta F}")
ax.plot(tecco, df.nosfcadjust,  c = colors[3], label = L"\theta^{\Delta T}")
ax.plot(tecco, df.iter0_bulkformula,  c = colors[4], label = L"\theta^{0}")

c1, c2= round.(coefs, digits = 2)
ax.set_xlabel("time")
ax.set_ylabel(L"^\circ"*"C")
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "θDeep_" * region * "_" * suffix * ".png",
dpi = 1000)
close("all")

labels_L = [L"\theta^{129}", L"\theta^{\Delta F}",  L"\theta^{\Delta T}", L"\theta^{0}"]
diff_df = df[[312], :]
diff_df .= Matrix(df[[312], :]) - Matrix(df[[1], :])
fig, ax = plt.subplots(1, figsize = (15, 15))
for (i, name) in enumerate(keys(shortnames))
    ax.barh(labels_L[i], diff_df[!, name], color= colors[i], label = labels_L[i])
end
# ax.set_yticks(1:4, labels=[L"\theta^{129}"])

ax.set_xlabel("Δθ (" *L"^\circ"*"C)")
ax.set_title("Potential Temperature Change by experiment")
ax.legend()
# sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "θChange_" * region * "_" * suffix * ".png", dpi = 1000)
