using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyCall
using DataFrames, GLM
#df = DataFrame(θz); df = df .- df["iter0_bulkformula] #initiate dataframe & do mean about iter0
#
cm = pyimport("cmocean.cm");colorway = cm.balance;
θsurf_noszn = ref_month_1( θsurf_noszn,"iter0_bulkformula")
θdeep_noszn = ref_month_1( θz,"iter0_bulkformula")
reference_text = "\n Baseline: iter0_bf"
# df = DataFrame(θdeep_noszn)
df = DataFrame(θz); df = df .- df.iter0_bulkformula;
#initiate dataframe & do mean about iter0

fm1 = @formula(iter129_bulkformula ~ iter0_bulkformula + noinitadjust + nosfcadjust)
linearRegressor1 = lm(fm1, Float64.(df))
fm2 = @formula(iter0_bulkformula ~ iter129_bulkformula + noinitadjust + nosfcadjust)
linearRegressor2 = lm(fm2, Float64.(df))
iter129_predict = predict(linearRegressor1)
strings = string.(round.(coef(linearRegressor1), digits = 3))
form1 = "iter129_ref = " * strings[2] * "iter0_ref +" * strings[3] * "sfc_adjust_ref + " * strings[4] * "initadjust_ref +" *
string(strings[1])
iter0_predict = predict(linearRegressor2)
strings = string.(round.(coef(linearRegressor2), digits = 3))
form2 = "iter0_ref = " * strings[2] * "iter129_ref +" * strings[3] * "sfc_adjust_ref + " * strings[4] * "initadjust_ref +" *
string(strings[1])

fig, axs = plt.subplots(2, 1, figsize = (12, 8))
fig.suptitle("Applying Linear Regression using sensitivity Experiments \n Deep " * region * reference_text)
axs[1].plot(tecco, df[!, "iter129_bulkformula"], c = "#1f77b4", label = "data")
axs[1].plot(tecco, iter129_predict, c = "#1f77b4", linestyle = "--", label = "prediction")
axs[1].set_title(form1)
axs[2].plot(tecco, df[!, "iter0_bulkformula"], c = "#d62728", label = "data")
axs[2].plot(tecco, iter0_predict, c = "#d62728", linestyle = "--", label = "prediction")
axs[2].set_title(form2)
axs[1].set_xlabel("time")
axs[1].set_ylabel("ºC")
axs[2].set_xlabel("time")
axs[2].set_ylabel("ºC")
axs[1].legend()
axs[2].legend()

fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θDeepPredictionsiter0_" * region * suffix * ".png")

df = DataFrame(θsurf_noszn)
fm1 = @formula(iter129_bulkformula ~ iter0_bulkformula + noinitadjust + nosfcadjust)
linearRegressor1 = lm(fm1, df)
fm2 = @formula(iter0_bulkformula ~ iter129_bulkformula + noinitadjust + nosfcadjust)
linearRegressor2 = lm(fm2, df)
iter129_predict = predict(linearRegressor1)
strings = string.(round.(coef(linearRegressor1), digits = 3))
form1 = "iter129 = " * strings[2] * "iter0 +" * strings[3] * "sfc_adjust + " * strings[4] * "initadjust"
iter0_predict = predict(linearRegressor2)
strings = string.(round.(coef(linearRegressor2), digits = 3))
form2 = "iter0 = " * strings[2] * "iter129 +" * strings[3] * "sfc_adjust + " * strings[4] * "initadjust"

fig, axs = plt.subplots(2, 1, figsize = (12, 7))
fig.suptitle("Applying Linear Regression using sensitivity Experiments \n Surf " * region * reference_text)
axs[1].plot(tecco, θsurf_noszn["iter129_bulkformula"], c = "#1f77b4", label = "data")
axs[1].plot(tecco, iter129_predict, c = "#1f77b4", linestyle = "--", label = "prediction")
axs[1].set_title(form1)
axs[2].plot(tecco, θsurf_noszn["iter0_bulkformula"], c = "#d62728", label = "data")
axs[2].plot(tecco, iter0_predict, c = "#d62728", linestyle = "--", label = "prediction")
axs[2].set_title(form2)
axs[1].set_xlabel("time")
axs[1].set_ylabel("ºC")
axs[2].set_xlabel("time")
axs[2].set_ylabel("ºC")
axs[1].legend()
axs[2].legend()

fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θSurfPredictions_" * region * suffix * ".png")

#Now see if there is any correlation 
#(plot surface warming as a trends as a function of deep cooling)