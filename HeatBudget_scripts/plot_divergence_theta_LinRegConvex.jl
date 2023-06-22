
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions, NonNegLeastSquares
using PyCall
using DataFrames, GLM
using Convex, SCS
const MOI = Convex.MOI

tecco = 1992+1/24:1/12:2018; tecco = collect(tecco)
tecco95 = 2005 .> tecco .> 1995
tecco = tecco[tecco95]

df = DataFrame(θz)[tecco95, :];
ΔF = df.iter129_bulkformula .- df.nosfcadjust 
ΔT = df.iter129_bulkformula .- df.noinitadjust 
A = [ΔF ΔT ΔF.^2 ΔT.^2]
A_iter129 = deepcopy(A) 
b =  df.iter129_bulkformula .- df.iter0_bulkformula

params = Variable(4,  Positive())
add_constraint!(params, params <= 1)
problem = minimize(sumsquares(A*params -b))
solve!(problem, MOI.OptimizerWithAttributes(SCS.Optimizer,  
"eps_abs" => 1e-13, "eps_rel" => 1e-13), silent_solver = true)
coefs = round.(params.value, digits = 3)

fig, ax = plt.subplots(1, figsize = (15, 10))
ax.plot(tecco, ΔF, c = colors[2], label = L"\Delta F = θ^{\Delta T, \Delta F} - θ^{\Delta T}")
ax.plot(tecco, ΔT, c = colors[3], label = L"\Delta T = θ^{\Delta T, \Delta F} - θ^{\Delta F}")
ax.legend();sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
ax.set_xlabel("time");ax.set_ylabel(L"^\circ"*"C")
path = "/OHC_Divergence/Reconstructions/" * "Δθ_about129_" * region * "_" * suffix * ".png"
fig.tight_layout(); fig.savefig(plotsdir() * path, dpi = 1000)

iter0_predict = df.iter129_bulkformula .- (A_iter129 * coefs) 
iter129_predict = df.iter0_bulkformula .+  (A_iter129 * coefs)

fig, ax = plt.subplots(1, figsize = (15, 10))
ax.plot(tecco, df.iter129_bulkformula, c = colors[1], label = L"\theta^{\Delta T, \Delta F}")
ax.plot(tecco, df.iter0_bulkformula,  c = colors[4], label = L"\theta^{0}")
ax.plot(tecco, iter129_predict, c = colors[1], linestyle = "--")
ax.plot(tecco, iter0_predict, c = colors[4], linestyle = "--")
ax.plot(tecco, NaN.*tecco, c = "k", linestyle = "--", label = "prediction")

c1, c2, c3, c4 = round.(coefs, digits = 2)
# c1 = 0.68; c2 = 0.85
ax.set_title(L"θ^{0}" * " ≈  $c1 ΔF + $c2 ΔT + $c3 ΔF^2 + $c4 ΔT^2+ " *L"θ^{\Delta T, \Delta F}" * "\n")

ax.set_xlabel("time");ax.set_ylabel(L"^\circ"*"C")
ax.legend();sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
path = "/OHC_Divergence/Reconstructions/" * "θPredCVX_about129_" * region * "_" * suffix * ".png"
fig.tight_layout(); fig.savefig(plotsdir() * path, dpi = 1000)

ΔF = df.noinitadjust .- df.iter0_bulkformula
ΔT = df.nosfcadjust .- df.iter0_bulkformula 
b =  df.iter129_bulkformula .- df.iter0_bulkformula

A = [ΔF ΔT ΔF.^2 ΔT.^2]

params = Variable(4,  Positive())
add_constraint!(params, params <= 1)
problem = minimize(sumsquares(A*params -b))
solve!(problem, MOI.OptimizerWithAttributes(SCS.Optimizer,  
"eps_abs" => 1e-13, "eps_rel" => 1e-13), silent_solver = true)
coefs = Float32.(params.value)
c1, c2, c3, c4 = round.(coefs, digits = 2)

iter0_predict = df.iter129_bulkformula .- (A * coefs) 
iter129_predict = df.iter0_bulkformula .+  (A * coefs)

fig, ax = plt.subplots(1, figsize = (15, 10))
ax.plot(tecco, ΔF, c = colors[2], label = L"\Delta F = θ^{\Delta F} - θ^{0}")
ax.plot(tecco, ΔT, c = colors[3], label = L"\Delta T = θ^{\Delta T} - θ^{0}")
ax.legend();sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
ax.set_xlabel("time");ax.set_ylabel(L"^\circ"*"C")
path = "/OHC_Divergence/Reconstructions/" * "Δθ_about0_" * region * "_" * suffix * ".png"
fig.tight_layout(); fig.savefig(plotsdir() * path, dpi = 1000)

fig, ax = plt.subplots(1, figsize = (15, 10))
ax.plot(tecco, df.iter129_bulkformula, c = colors[1], label = L"\theta^{\Delta T, \Delta F}")
ax.plot(tecco, df.iter0_bulkformula,  c = colors[4], label = L"\theta^{0}")
ax.plot(tecco, iter129_predict, c = colors[1], linestyle = "--")
ax.plot(tecco, iter0_predict, c = colors[4], linestyle = "--")
ax.plot(tecco, NaN.*tecco, c = "k", linestyle = "--", label = "prediction")

ax.set_title(L"θ^{\Delta T, \Delta F}" * " ≈  $c1 ΔF + $c2 ΔT + $c3 ΔF^2 + $c4 ΔT^2+ " *L"θ^{0}" * "\n")
ax.set_xlabel("time"); ax.set_ylabel(L"^\circ"*"C")
ax.legend(); sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
path = "/OHC_Divergence/Reconstructions/" * "θPredCVX_about0_" * region * "_" * suffix * ".png"
fig.tight_layout(); fig.savefig(plotsdir() * path, dpi = 1000)


fig, ax = plt.subplots(1, figsize = (15, 10))
ax.plot(tecco, A[:, 1], c = colors[2], linestyle = "--", label = L"θ^{\Delta F} - θ^{0}")
ax.plot(tecco, A[:, 2], c = colors[3], linestyle = "--", label = L"θ^{\Delta T} - θ^{0}")
ax.plot(tecco, A_iter129[:, 1], c = colors[2], label = L" θ^{\Delta T, \Delta F} - θ^{\Delta T}")
ax.plot(tecco, A_iter129[:, 2], c = colors[3], label = L" θ^{\Delta T, \Delta F} - θ^{\Delta F}")

ax.legend();sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
ax.set_xlabel("time");ax.set_ylabel(L"^\circ"*"C")
path = "/OHC_Divergence/Reconstructions/" * "Δθ_comp_" * region * "_" * suffix * ".png"
fig.tight_layout(); fig.savefig(plotsdir() * path, dpi = 1000)

A = [A_iter129; A]
b = [b; b]

params = Variable(4,  Positive())
add_constraint!(params, params <= 1)
problem = minimize(sumsquares(A*params -b))
solve!(problem, MOI.OptimizerWithAttributes(SCS.Optimizer,  
"eps_abs" => 1e-13, "eps_rel" => 1e-13), silent_solver = true)
coefs = round.(params.value, digits = 3)
c1, c2, c3, c4 = round.(coefs, digits = 2)

iter0_predict = [df.iter129_bulkformula;df.iter129_bulkformula]  .- (A * coefs) 
iter129_predict = [df.iter0_bulkformula; df.iter0_bulkformula] .+  (A * coefs)

nt = Int(length(iter129_predict)/2)

fig, ax = plt.subplots(1, figsize = (15, 10))
ax.plot(tecco, df.iter129_bulkformula, c = colors[1], label = L"\theta^{\Delta T, \Delta F}")
ax.plot(tecco, df.iter0_bulkformula,  c = colors[4], label = L"\theta^{0}")
ax.plot(tecco, iter129_predict[1:nt], c = colors[1], linestyle = "--")
ax.plot(tecco, iter0_predict[1:nt], c = colors[4], linestyle = "--")
ax.plot(tecco, iter129_predict[nt+1:2*nt], c = colors[1], linestyle = "-.")
ax.plot(tecco, iter0_predict[nt+1:2*nt], c = colors[4], linestyle = "-.")
ax.plot(tecco, NaN.*tecco, c = "k", linestyle = "--", label = "prediction")
ax.plot(tecco, NaN.*tecco, c = "k", linestyle = "-.", label = "prediction")

ax.set_title(L"θ^{\Delta T, \Delta F}" * " ≈  $c1 ΔF + $c2 ΔT + $c3 ΔF^2 + $c4 ΔT^2+ " *L"θ^{0}" * "\n")
ax.legend();sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, -.15), ncol=5)
ax.set_xlabel("time");ax.set_ylabel(L"^\circ"*"C")
path = "/OHC_Divergence/Reconstructions/" * "θPredCVX_aboutboth_" * region * "_" * suffix * ".png"
fig.tight_layout(); fig.savefig(plotsdir() * path, dpi = 1000)