#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings,Plots
import PyPlot as plt
using PyCall
import NaNMath as nm
cm = pyimport("cmocean.cm");colorway = cm.balance;

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = OHC_helper.wet_pts(Γ)

regular_mask = LLCcropC(ocean_mask,γ)


regular_mask[regular_mask .== 0] .= NaN

#do not need to add anything to the lons, this function takes care of it 
regular_λ = LLCcropC(λ,γ)
regular_ϕ = LLCcropC(ϕ,γ)

fig, ax = plt.subplots( figsize=(15,10))
ax.scatter(regular_λ[:], regular_ϕ[:], regular_mask[:])
fig
