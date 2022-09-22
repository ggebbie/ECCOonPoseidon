
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions, NonNegLeastSquares
using PyCall
using DataFrames, GLM

df = DataFrame(θz);
df_anom = round.(Float64.(df .- df.iter129_bulkformula), digits = 8)
ols = lm(@formula(iter0_bulkformula ~ 0 + noinitadjust + nosfcadjust), df_anom)
ΔF, ΔT = round.(coef(ols), digits = 3)
iter0_predict = (ΔF .* df_anom.noinitadjust) .+ (ΔT .* df_anom.nosfcadjust) .+ 
                df.iter129_bulkformula 
iter129_predict = df.iter0_bulkformula .- (ΔF .* df_anom.noinitadjust) .- (ΔT .* df_anom.nosfcadjust) 
noinit_predict = (df_anom.noinitadjust) .+ df.iter129_bulkformula 
nosfc_predict = (df_anom.nosfcadjust) .+ df.iter129_bulkformula 

fig, ax = plt.subplots(1, figsize = (15, 10))
ax.set_title(L"θ^{0}" * " ≈  $ΔF ΔF + $ΔT ΔT + " *L"θ^{129}" * "\n")

# fig.suptitle("Applying NNLinear Regression using sensitivity Experiments \n " * suffix * " " * region * reference_text)
ax.plot(tecco, df.iter129_bulkformula, c = colors[1], label = L"\theta^{129}")
ax.plot(tecco, df.iter0_bulkformula,  c = colors[4], label = L"\theta^{0}")

ax.plot(tecco, iter129_predict, c = colors[1], linestyle = "--")
ax.plot(tecco, iter0_predict, c = colors[4], linestyle = "--")

# ax.fill_between(tecco, iter129_predict, df.iter129_bulkformula, color = "grey", 
# alpha = 0.2)
ax.plot(tecco, NaN.*tecco, c = "k", linestyle = "--", label = "prediction")

ax.set_xlabel("time")
ax.set_ylabel(L"^\circ"*"C")
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "θPredictionsAbout129_" * region * "_" * suffix * ".png",
dpi = 1000)

df_anom = round.(Float64.(df .- df.iter0_bulkformula), digits = 8)
ols = lm(@formula(iter129_bulkformula ~ 0 + noinitadjust + nosfcadjust), df_anom)
ΔF, ΔT = round.(coef(ols), digits = 3)
iter129_predict = (ΔF .* df_anom.noinitadjust) .+ (ΔT .* df_anom.nosfcadjust) .+ 
                df.iter0_bulkformula 
iter0_predict = df.iter129_bulkformula .- (ΔF .* df_anom.noinitadjust) .- (ΔT .* df_anom.nosfcadjust) 
noinit_predict = (df_anom.noinitadjust) .+ df.iter0_bulkformula 
nosfc_predict = (df_anom.nosfcadjust) .+ df.iter0_bulkformula 

# iter0_predict = Float32.(iter0_predict)

fig, ax = plt.subplots(1, figsize = (15, 10))
ax.set_title(L"θ^{129}" * " ≈  $ΔF ΔF + $ΔT ΔT + " *L"θ^{0}" * "\n")

# fig.suptitle("Applying NNLinear Regression using sensitivity Experiments \n " * suffix * " " * region * reference_text)
ax.plot(tecco, df.iter129_bulkformula, c = colors[1], label = L"\theta^{129}")
ax.plot(tecco, df.iter0_bulkformula,  c = colors[4], label = L"\theta^{0}")

ax.plot(tecco, iter129_predict, c = colors[1], linestyle = "--")
ax.plot(tecco, iter0_predict, c = colors[4], linestyle = "--")

# ax.fill_between(tecco, iter129_predict, df.iter129_bulkformula, color = "grey", 
# alpha = 0.2)
ax.plot(tecco, NaN.*tecco, c = "k", linestyle = "--", label = "prediction")

ax.set_xlabel("time")
ax.set_ylabel(L"^\circ"*"C")
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "θPredictionsAbout0_" * region * "_" * suffix * ".png",
dpi = 1000)




# fig, ax = plt.subplots(1, figsize = (17, 10))
# ax.set_title("residual of least squares estimator")
# iter129_resid = df.iter129_bulkformula .- iter129_predict
# noinit_resid = df.noinitadjust .- noinit_predict
# nosfc_resid = df.nosfcadjust .- nosfc_predict
# iter0_resid = df.iter0_bulkformula .- iter0_predict

# ax.scatter(tecco, iter129_resid, color = colors[1], label = L"\theta^{129}")
# ax.scatter(tecco, noinit_resid,color = colors[2], label = L"\theta^{\Delta F}")
# ax.scatter(tecco, nosfc_resid, color = colors[3], label = L"\theta^{\Delta T}")
# ax.scatter(tecco, iter0_resid, color = colors[4], label = L"\theta^{0}")

# ax.set_xlabel("time")
# ax.set_ylabel("" * L"^\circ"*"C")
# ax.legend()
# sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
# fig.tight_layout()
# fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "θPredictionsMLRResid_" * region * "_" * suffix * ".png",
# dpi = 1000)
