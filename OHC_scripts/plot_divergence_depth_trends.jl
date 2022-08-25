using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, JLD2,
NCDatasets, NetCDF, Printf
using .OHC_helper
using PyCall

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", 
            font_scale=1.50,
              palette = sns.color_palette("deep"));#sns.set_context("talk")
outdir = "ECCO_vars/"

filedir = "ECCO_vars/"
filename = "THETAs"
filename2 = "ETAN"
θDiff = Dict()
smush_depths =  Float32.(OHC_helper.smush(cell_depths))
smush_depths[findall(smush_depths .== 0)] = Inf32
inv_depths = 1 ./ smush_depths
for (key,values) in shortnames
    expname = key; println(key)
    @time θ = load_object(datadir(filedir*filename*"_"*expname*".jld2")); 
    @time ETAN = load_object(datadir(filedir*filename2*"_"*expname*".jld2"));  
    # var_exp  = Vector{MeshArrays.gcmarray{Float32, 1, Matrix{Float32}}}(undef, 0)
    sθz = zeros(50, 312)
    for tt in [1, 312]
        s1 = @views ((ETAN[tt] .* inv_depths) )
        s1 .+= 1
        s1 .*= msk 
        sθ = @views θ[tt] .* s1
        for k in 1:size(sθz, 1)
            sθz[k, tt] = volume_mean(sθ[:, k]; weights = cell_volumes[:, k])
        end
    end
    θDiff[expname] = sθz[:, 312] .- sθz[:, 1]
    @time GC.gc(true)
end


df = Float32.(DataFrame(θDiff))
fig, ax = plt.subplots(1, figsize = (15, 10))
fig.suptitle("Temperature Change in \n " * suffix * " " * region * reference_text)
i = 0
ax.plot(df.iter129_bulkformula, z, c = colors[1], label = L"\theta^{129}", marker = "o")
ax.plot(df.noinitadjust, z,  c = colors[2], label = L"\theta^{\Delta F}", marker = "o")
ax.plot(df.nosfcadjust, z,  c = colors[3], label = L"\theta^{\Delta T}", marker = "o")
ax.plot(df.iter0_bulkformula, z,  c = colors[4], label = L"\theta^{0}", marker = "o")

ax.set_xlabel("Δθ [" * L"^\circ" * "C]")
ax.set_ylabel("Depth [m]")
ax.set_xlim((-0.075, 0.075))
ax.set_ylim((-6000, -500))
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Depth_Plots/" * "Δθ" * region *".png",
dpi = 1000)

df = Float32.(DataFrame(θDiff))
fig, ax = plt.subplots(1, figsize = (15, 15))
fig.suptitle("Avg. North Pacific Potential Temperature Change")
i = 0
ax.plot(df.iter129_bulkformula, z, c = colors[1], label = L"ECCO", marker = "o")
ax.set_xlabel("Δθ [" * L"^\circ" * "C]")
ax.set_ylabel("Depth [m]")
ax.set_xlim((-0.075, 0.075))
ax.set_ylim((-6000, -500))
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Depth_Plots/" * "Δθ_ECCO" * region *".png",
dpi = 1000)


