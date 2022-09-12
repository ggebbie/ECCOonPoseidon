using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
cm = pyimport("cmocean.cm");colorway = cm.balance;
mpatches = pyimport("matplotlib.patches")
lines = pyimport("matplotlib.lines")
Line2D = lines.Line2D


fig, ax = plt.subplots(1, 1, figsize = (12, 9))
ax.set_title(region *  " θ̄,  z=2-3km")
plot_ts!(θz, tecco, shortnames, ignore_list, ax; ylabel =  L" ^\circ C", linestyle = "-")
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "θAverage_" * region * suffix * ".png")
close("all")

