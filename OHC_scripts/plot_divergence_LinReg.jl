
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
# df = df .- df.iter0_bulkformula[1];
A = zeros(312, 2)
A .= Matrix{Float32}(df[!, ["noinitadjust",  "nosfcadjust"]]) .- df.iter129_bulkformula

b =  df.iter0_bulkformula .- df.iter129_bulkformula
b2 = df.nosfcadjust .- df.iter129_bulkformula
b1 = df.noinitadjust .- df.iter129_bulkformula

using Convex, SCS
params = Variable(2)
# params.value = [1.2; -0.2]
add_constraint!(params, params ≤ 1)
problem = minimize(sumsquares(A*params -b))
# sumsquares(A[:, 1]*params[1] -b1) + 
# sumsquares(A[:, 2]*params[2] -b2))
const MOI = Convex.MOI
solve!(problem, 
MOI.OptimizerWithAttributes(SCS.Optimizer,  
"eps_abs" => 1e-15, "eps_rel" => 1e-10),
warmstart = false)
coefs = Float32.(params.value)
iter129_predict = df.iter0_bulkformula .-  A * coefs

iter0_predict = A * coefs .+ df.iter129_bulkformula
noinit_predict = A[:, 1] * coefs[1] .+ df.iter129_bulkformula
nosfc_predict = A[:, 2] * coefs[2] .+ df.iter129_bulkformula

fig, ax = plt.subplots(1, figsize = (15, 10))
# fig.suptitle("Applying NNLinear Regression using sensitivity Experiments \n " * suffix * " " * region * reference_text)
ax.plot(tecco, df.iter129_bulkformula, c = colors[1], label = L"\theta^{129}")
ax.plot(tecco, df.noinitadjust,  c = colors[2], label = L"\theta^{\Delta F}")
ax.plot(tecco, df.nosfcadjust,  c = colors[3], label = L"\theta^{\Delta T}")
ax.plot(tecco, df.iter0_bulkformula,  c = colors[4], label = L"\theta^{0}")

ax.plot(tecco, iter129_predict, c = colors[1], linestyle = "--")
ax.plot(tecco, noinit_predict,c = colors[2], linestyle = "--")
ax.plot(tecco, nosfc_predict, c = colors[3], linestyle = "--")
ax.plot(tecco, iter0_predict, c = colors[4], linestyle = "--")

# ax.plot(tecco, df.noinitadjust .+ 0.004587216254992363, c = colors[1], linestyle = "-.")
# ax.plot(tecco, df.iter129_bulkformula .- 0.004587216254992363, c = colors[4], linestyle = "-.")

# ax.fill_between(tecco, iter129_predict, df.iter129_bulkformula, color = "grey", 
# alpha = 0.2)
ax.plot(tecco, NaN.*tecco, c = "k", linestyle = "--", label = "prediction")

c1, c2= round.(coefs, digits = 2)
ax.set_title("iter0 ≈  $c1 ΔF + $c2 Δθ₀ + iter129")
ax.set_xlabel("time")
ax.set_ylabel(L"^\circ"*"C")
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "θPredictionsLIN_" * region * "_" * suffix * ".png",
dpi = 1000)


using GLM
df = DataFrame(θz);
df_anom = round.(Float64.(df .- df.iter129_bulkformula), digits = 8)
ols = lm(@formula(iter0_bulkformula ~ 0 + noinitadjust + nosfcadjust), df_anom)
ΔF, ΔT = round.(coef(ols), digits = 3)
iter0_predict = (ΔF .* df_anom.noinitadjust) .+ (ΔT .* df_anom.nosfcadjust) .+ 
                df.iter129_bulkformula 
iter129_predict = df.iter0_bulkformula .- (ΔF .* df_anom.noinitadjust) .- (ΔT .* df_anom.nosfcadjust) 
noinit_predict = (ΔF .* df_anom.noinitadjust) .+ df.iter129_bulkformula 
nosfc_predict = (ΔF .* df_anom.nosfcadjust) .+ df.iter129_bulkformula 

# iter0_predict = Float32.(iter0_predict)

fig, ax = plt.subplots(1, figsize = (15, 10))
ax.set_title("iter0 ≈  $ΔF ΔF + $ΔT Δθ₀ + iter129")

# fig.suptitle("Applying NNLinear Regression using sensitivity Experiments \n " * suffix * " " * region * reference_text)
ax.plot(tecco, df.iter129_bulkformula, c = colors[1], label = L"\theta^{129}")
ax.plot(tecco, df.noinitadjust,  c = colors[2], label = L"\theta^{\Delta F}")
ax.plot(tecco, df.nosfcadjust,  c = colors[3], label = L"\theta^{\Delta T}")
ax.plot(tecco, df.iter0_bulkformula,  c = colors[4], label = L"\theta^{0}")

ax.plot(tecco, iter129_predict, c = colors[1], linestyle = "--")
ax.plot(tecco, noinit_predict,c = colors[2], linestyle = "--")
ax.plot(tecco, nosfc_predict, c = colors[3], linestyle = "--")
ax.plot(tecco, iter0_predict, c = colors[4], linestyle = "--")

# ax.fill_between(tecco, iter129_predict, df.iter129_bulkformula, color = "grey", 
# alpha = 0.2)
ax.plot(tecco, NaN.*tecco, c = "k", linestyle = "--", label = "prediction")

ax.set_xlabel("time")
ax.set_ylabel(L"^\circ"*"C")
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "θPredictionsMLR_" * region * "_" * suffix * ".png",
dpi = 1000)

fig, ax = plt.subplots(1, figsize = (17, 10))
ax.set_title("residual of least squares estimator")
iter129_resid = df.iter129_bulkformula .- iter129_predict
noinit_resid = df.noinitadjust .- noinit_predict
nosfc_resid = df.nosfcadjust .- nosfc_predict
iter0_resid = df.iter0_bulkformula .- iter0_predict

ax.scatter(tecco, iter129_resid, color = colors[1], label = L"\theta^{129}")
ax.scatter(tecco, noinit_resid,color = colors[2], label = L"\theta^{\Delta F}")
ax.scatter(tecco, nosfc_resid, color = colors[3], label = L"\theta^{\Delta T}")
ax.scatter(tecco, iter0_resid, color = colors[4], label = L"\theta^{0}")

ax.set_xlabel("time")
ax.set_ylabel("" * L"^\circ"*"C")
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "θPredictionsMLRResid_" * region * "_" * suffix * ".png",
dpi = 1000)


using GLM
df = DataFrame(θz);
df_anom = round.(Float64.(df .- df.iter0_bulkformula), digits = 8)
ols = lm(@formula(iter129_bulkformula ~ 0 + noinitadjust + nosfcadjust), df_anom)
ΔF, ΔT = round.(coef(ols), digits = 3)
iter129_predict = (ΔF .* df_anom.noinitadjust) .+ (ΔT .* df_anom.nosfcadjust) .+ 
                df.iter0_bulkformula 
iter0_predict = df.iter129_bulkformula .- (ΔF .* df_anom.noinitadjust) .- (ΔT .* df_anom.nosfcadjust) 
noinit_predict = (ΔF .* df_anom.noinitadjust) .+ df.iter0_bulkformula 
nosfc_predict = (ΔF .* df_anom.nosfcadjust) .+ df.iter0_bulkformula 

# iter0_predict = Float32.(iter0_predict)

fig, ax = plt.subplots(1, figsize = (15, 10))
ax.set_title("iter129 ≈  $ΔF ΔF + $ΔT Δθ₀ + iter0")

# fig.suptitle("Applying NNLinear Regression using sensitivity Experiments \n " * suffix * " " * region * reference_text)
ax.plot(tecco, df.iter129_bulkformula, c = colors[1], label = L"\theta^{129}")
ax.plot(tecco, df.noinitadjust,  c = colors[2], label = L"\theta^{\Delta F}")
ax.plot(tecco, df.nosfcadjust,  c = colors[3], label = L"\theta^{\Delta T}")
ax.plot(tecco, df.iter0_bulkformula,  c = colors[4], label = L"\theta^{0}")

ax.plot(tecco, iter129_predict, c = colors[1], linestyle = "--")
ax.plot(tecco, noinit_predict,c = colors[2], linestyle = "--")
ax.plot(tecco, nosfc_predict, c = colors[3], linestyle = "--")
ax.plot(tecco, iter0_predict, c = colors[4], linestyle = "--")

# ax.fill_between(tecco, iter129_predict, df.iter129_bulkformula, color = "grey", 
# alpha = 0.2)
ax.plot(tecco, NaN.*tecco, c = "k", linestyle = "--", label = "prediction")

ax.set_xlabel("time")
ax.set_ylabel(L"^\circ"*"C")
ax.legend()
sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "θPredictionsMLR2_" * region * "_" * suffix * ".png",
dpi = 1000)




