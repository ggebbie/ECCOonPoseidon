include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf
using .OHC_helper
import PyPlot: plt  # important!
using PyCall
# using Pandas
@pyimport pandas as pd
@pyimport seaborn as sns
#sns = Seaborn 
sns.set_theme(); pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;

sns.set_theme("poster")

# Load an example dataset
tips = pd.DataFrame(Dict("total_bill" => 1:100, "tip" => 1:100, "smoker" => iseven.(collect(1:100))))
tips = tips.rename_axis("time")

# Create a visualization
p = sns.relplot(
    data=tips,
    x="total_bill", y="tip",
    hue="smoker", style="smoker",
)
savefig(plotsdir("seaborn_test.png"))