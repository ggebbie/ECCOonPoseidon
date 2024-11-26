include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour, PyCall, 
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings, DataFrames
import NaNMath as nm
# import PyPlot as plt

include(srcdir("config_exp.jl"))
@pyimport seaborn as sns;

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))

@pyimport cmocean.cm as cmos


tecco = Float32.(collect(1992+1/24:1/12:2018))
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)

region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ctrl_vol = sum(cell_volumes[:, lvls])

vars =  ["iter0_bulkformula",
        "only_init", "only_kappa", "only_wind", "only_buoyancy"]
WEul = Dict(); WBol = Dict()
ΨEul = Dict();
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()
z_ref = findall( -2800 .<= -z[:].<= -2300)[1]

for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2"))
    Ψ_EulBol = read_file["Ψ_exp_timeseries"]
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_Eul = read_file["Ψ_exp_timeseries"]
    ΨEul[expname] = 1 .* Ψ_Eul
    ϕ_avg = read_file["ϕ_avg"]
    ϕ_ref = findall( 23 .<= ϕ_avg .<= 40)

    Ψ_Bol = Ψ_EulBol .- Ψ_Eul

    WEul[expname] = -Float32.(Ψ_Eul[z_ref, ϕ_ref[1], :][:])
    WBol[expname] = -Float32.(Ψ_Bol[z_ref, ϕ_ref[1], :][:])

    fname = datadir("native/" * expname * "NPAC_THETA_budget_Wθ_eul_bol_" * suffix * ".jld2")
    vars = jldopen(fname)["dθ"]

    wθ_resid[expname] = Float32.(vars["wθBot"] .- vars["wθTop"])
    wθ_bol[expname] = Float32.(vars["wθBot_bol"] .- vars["wθTop_bol"])
    wθ_eul[expname] = Float32.(wθ_resid[expname] .- wθ_bol[expname])

end

effect_exp =  ["only_wind", "only_kappa", "only_init"]

[WEul[expt] .-= WEul["iter0_bulkformula"] for expt in effect_exp]
[WBol[expt] .-= WBol["iter0_bulkformula"] for expt in effect_exp]

[wθ_resid[expt] .-= wθ_resid["iter0_bulkformula"] for expt in effect_exp]
# [wθ_eul[expt] .-= wθ_eul["iter0_bulkformula"] for expt in effect_exp]
# [wθ_bol[expt] .-= wθ_bol["iter0_bulkformula"] for expt in effect_exp]

vars =  ["only_wind",  "only_init", "only_kappa"]
var(x) = std(x; corrected = false)^2
r2(y, ỹ) = 1 - (mean( (y .- ỹ ).^ 2) / var(y))

lin_regress(x, y) = trend_matrices(x; remove_mean = false)[2] * y[:]

R² = DataFrame(Dict("Eulerian" => zeros(3), "Bolus" => zeros(3), "expname" => Array{String}(undef, 3)))
EulRegCoef = DataFrame(Dict("a" => zeros(3), "b" => zeros(3), "expname" => Array{String}(undef, 3)))
BolRegCoef = DataFrame(Dict("a" => zeros(3), "b" => zeros(3), "expname" => Array{String}(undef, 3)))

EulReg = DataFrame(Dict("only_wind" => zeros(312), "only_init" => zeros(312), "only_kappa" => zeros(312)))
BolReg = 1 .* EulReg

fig, ax = plt.subplots(1, 3, figsize = (15, 5), sharex = true, sharey = true)

for (i, expname) in enumerate(vars)
    println(expname)

    EulRegCoef.expname[i] = expname
    BolRegCoef.expname[i] = expname

    x1 = (100 * 86400 * 365 * 10) * WEul[expname] ./ ctrl_vol
    x2 = (100 * 86400 * 365 * 10) * WBol[expname] ./ ctrl_vol

    y1 = (100 * 86400 * 365 * 10) * wθ_resid[expname][:]
    # ax[1].plot(tecco, cumsum(vcat(0, y1 .* (100 * 2.628e+6)))[1:end-1], c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    ax[1].plot(tecco, y1, c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)

    a, b = lin_regress(x1, y1)
    EulRegCoef.a[i] = a
    EulRegCoef.b[i] = b
    ỹ1 = a .+ (b .* x1); EulReg[!, expname] .= ỹ1

    println("Eulerian")
    println("slope = ", b, " ")
    R².Eulerian[i] = 100 * r2(y1, ỹ1)
    println("R² = ", R².Eulerian[i])
    # ax[2].plot(tecco, cumsum(vcat(0, ỹ1 .* (100 * 2.628e+6)))[1:end-1], c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    ax[2].plot(tecco, ỹ1, c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)

    a, b = lin_regress(x2, y1)
    BolRegCoef.a[i] = a
    BolRegCoef.b[i] = b
    ỹ2 = a .+ (b .* x2); BolReg[!, expname] .= ỹ2
    println("Bolus")
    println("slope = ", b, " ")
    R².Bolus[i] = 100 * r2(y1, ỹ2)
    println("R² = ", R².Bolus[i])
    # ax[3].plot(tecco, cumsum(vcat(0, ỹ2 .* (100 * 2.628e+6)))[1:end-1], c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    ax[3].plot(tecco, ỹ2, c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
end
fig

EulRegCoef

BolRegCoef
R²
# EulReg_cum = 1 .* EulReg
# [EulReg_cum[:, i] .= cumsum(vcat(0, EulReg[:, i].* (2.628e+6)))[1:end-1]  for i in 1:3 ]

# BolReg_cum = 1 .* BolReg
# [BolReg_cum[:, i] .= cumsum(vcat(0, BolReg[:, i].* (100 * 2.628e+6)))[1:end-1]  for i in 1:3 ]

# wθ_df_cum = DataFrame(wθ_resid)
# [wθ_df_cum[:, i] .= cumsum(vcat(0, wθ_df_cum[:, i].* (100 * 2.628e+6)))[1:end-1]  for i in 1:6 ]

# wθ_df_cum, EulReg_cum, BolReg_cum
# wθ_df_cum
# ax[1].legend(frameon = false)