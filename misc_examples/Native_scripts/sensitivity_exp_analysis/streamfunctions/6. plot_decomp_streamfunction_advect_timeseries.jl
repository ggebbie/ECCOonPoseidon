include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour, PyCall, 
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings, DataFrames, Statistics
import NaNMath as nm
import PyPlot as plt

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
area_sum = sum(area .* PAC_msk)

vars =  ["iter0_bulkformula",
        "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]
WEul = Dict(); WBol = Dict()
ΨEul = Dict(); dθdz = Dict(); 
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()
z_ref = findall( -3400 .<= -z[:].<= -2000)

for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2"))
    Ψ_EulBol = read_file["Ψ_exp_timeseries"]
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_Eul = read_file["Ψ_exp_timeseries"]
    ΨEul[expname] = 1 .* Ψ_Eul
    ϕ_avg = read_file["ϕ_avg"]
    ϕ_ref = findall( 23 .<= ϕ_avg .<= 40)

    Ψ_Bol = Ψ_EulBol .- Ψ_Eul

    WEul[expname] = -Float32.(Ψ_Eul[z_ref, ϕ_ref[1], :])
    WBol[expname] = -Float32.(Ψ_Bol[z_ref, ϕ_ref[1], :])

    WEul[expname] = (WEul[expname][1:end-1, :]  .+ WEul[expname][2:end, :]) ./ 2
    WEul[expname] = sum(WEul[expname] .* ΔzF[z_ref[1:end-1]], dims = 1)[:] ./ ctrl_vol

    WBol[expname] = (WBol[expname][1:end-1, :]  .+ WBol[expname][2:end, :]) ./ 2
    WBol[expname] = sum(WBol[expname] .* ΔzF[z_ref[1:end-1]], dims = 1)[:] ./ ctrl_vol

    fname = datadir("native/" * expname * "NPAC_THETA_budget_Wθ_eul_bol_" * suffix * ".jld2")
    vars = jldopen(fname)["dθ"]

    wθ_resid[expname] = Float32.(vars["wθBot"] .- vars["wθTop"])
    wθ_bol[expname] = Float32.(vars["wθBot_bol"] .- vars["wθTop_bol"])
    wθ_eul[expname] = Float32.(wθ_resid[expname] .- wθ_bol[expname])

    fname = datadir(region * "_temperature_gradient_sens_exps.jld2")
    dθdz[expname] = Float32.(jldopen(fname)["adjust_exps"][expname])
end

effect_exp =  ["only_wind", "only_kappa", "only_init"]

[WEul[expt] .-= WEul["iter0_bulkformula"] for expt in effect_exp]
[WBol[expt] .-= WBol["iter0_bulkformula"] for expt in effect_exp]
[wθ_resid[expt] .-= wθ_resid["iter0_bulkformula"] for expt in effect_exp]

[dθdz[expt] .= (dθdz[expt] + dθdz["iter0_bulkformula"]) ./2 for expt in effect_exp]

vars =  ["only_wind",  "only_init", "only_kappa"]
var(x) = std(x; corrected = false)^2
r2(y, ỹ) = 1 - (mean( (y .- ỹ ).^ 2) / var(y))
cov_frac(y, ỹ) = cov(ỹ,y ) / (std(y) * std(y))

R² = DataFrame(Dict("Eulerian" => zeros(3), "Bolus" => zeros(3), 
"Residual" => zeros(3), "ϵ" => zeros(3),
"expname" => Array{String}(undef, 3)))
EulRegCoef = DataFrame(Dict("a" => zeros(3), "b" => zeros(3), "expname" => Array{String}(undef, 3)))
BolRegCoef = DataFrame(Dict("a" => zeros(3), "b" => zeros(3), "expname" => Array{String}(undef, 3)))

EulReg = DataFrame(Dict("only_wind" => zeros(312), "only_init" => zeros(312), "only_kappa" => zeros(312)))
BolReg = 1 .* EulReg

fig, ax = plt.subplots(1, 3, figsize = (15, 5), sharex = true, sharey = true)
for (i, expname) in enumerate(vars)

    println(expname)
    println(" ")

    BolRegCoef.expname[i] = expname
    R².expname[i] = expname

    b = mean(dθdz["iter0_bulkformula"])

    Δt = (tecco[2] - tecco[1] )/ 10
    x1 = -(86400 * 365 * 10) * WEul[expname] #
    x2 = -(86400 * 365 * 10) * WBol[expname]

    y1 = (86400 * 365 * 10) * wθ_resid[expname][:]

    println("Δθ = ", sum(Δt .* y1))

    # ax[1].plot(tecco, cumsum(vcat(0, y1 .* (100 * 2.628e+6)))[1:end-1], c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    ax[1].plot(tecco, y1, c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)

    BolRegCoef.b[i] = b
    EulRegCoef.b[i] = b

    ỹ = (b .* x1) .+ (b .* x2); 
    resid = y1 .- ỹ 
    R².ϵ[i] = cov_frac(y1, resid)
    println("Full Rsq")
    println("R² = ", cov(y1,ỹ ) / (std(y1) * std(y1)))
    R².Residual[i] = cov_frac(y1, ỹ)

        
    ỹ1 = (b .* x1); EulReg[!, expname] .= ỹ1
    println("Eulerian")
    println("slope = ", b, " ")
    R².Eulerian[i] = cov_frac(y1, ỹ1)
    println("R² = ", R².Eulerian[i])
    println("Δθ = ", sum(Δt .* ỹ1))

    # ax[2].plot(tecco, cumsum(vcat(0, ỹ1 .* (100 * 2.628e+6)))[1:end-1], c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    ax[2].plot(tecco, ỹ1, c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)

    # a, b = lin_regress(x2, y1)
    # BolRegCoef.a[i] = a
    # BolRegCoef.b[i] = b
    ỹ2 = (b .* x2); BolReg[!, expname] .= ỹ2
    println("Bolus")
    println("slope = ", b, " ")
    R².Bolus[i] = cov_frac(y1, ỹ2)
    println("R² = ", R².Bolus[i])
    println("Δθ = ", sum(Δt .* ỹ2))

    # ax[3].plot(tecco, cumsum(vcat(0, ỹ2 .* (100 * 2.628e+6)))[1:end-1], c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    ax[3].plot(tecco, ỹ2, c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)

    println(" ")
end
fig
R²
EulRegCoef
EulBolRegCoef
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