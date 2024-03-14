include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson,
    BenchmarkTools, LaTeXStrings, PyCall, DataFrames
import NaNMath as nm
import PyPlot as plt
import NumericalIntegration

cumul_integrate = NumericalIntegration.cumul_integrate

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
gs = pyimport("matplotlib.gridspec");
gridspec = gs.GridSpec
#blue, red, #green, orangeDataFrames

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 

uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)
tecco = 1992+1/24:1/12:2018

include(srcdir("plot_and_dir_config.jl"))
sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);

exps =  ["iter129_bulkformula", "iter0_bulkformula", "only_init", "only_wind"]
Hadvection = Dict(); diffusion = Dict(); GTF = Dict()
Vadvection = Dict()
ext_advection = Dict()
int_advection = Dict()
true_advection = Dict()
total = Dict()

temps = Dict()
for (i, expname) in enumerate(exps)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus_wextra" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    temps[expname] = vars["θ"]
    uθ = vars["VθSouth"]
    # uθ = vars["uvθ"]

    ∇wθ = vars["wθBot"] .- vars["wθTop"]
    Hadvection[expname] = (uθ)
    Vadvection[expname] = (∇wθ)
    ext_advection[expname] = uθ .+ ∇wθ
    int_advection[expname] = vars["wθ"] .+ vars["uvθ"] .- ext_advection[expname]
    true_advection[expname] = vars["wθ"] .+ vars["uvθ"]
    diffusion[expname] = (vars["κxyθ"] .+ vars["κzθ"])
    GTF[expname] = (vars["GTH"])
    total[expname] = GTF[expname] .+ diffusion[expname] .+ true_advection[expname]
end

for (i, expname) in enumerate(exps)
    error = abs.(int_advection[expname])
    truth = abs.(true_advection[expname])
    println(mean(error./ truth))
end
