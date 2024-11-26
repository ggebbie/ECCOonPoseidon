include("../../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour, PyCall, 
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings, DataFrames, Statistics
import NaNMath as nm
import PyPlot as plt
using DSP
include(srcdir("config_exp.jl"))
@pyimport seaborn as sns;

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))

@pyimport cmocean.cm as cmos

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/(7*12), fs = 1),Butterworth(4))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end


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
z_ref = [z_ref[1], z_ref[end] + 1]
for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2"))
    Ψ_EulBol = read_file["Ψ_exp_timeseries"]
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_Eul = read_file["Ψ_exp_timeseries"]
    ΨEul[expname] = 1 .* Ψ_Eul
    ϕ_avg = read_file["ϕ_avg"]
    ϕ_ref = findall( 23 .<= ϕ_avg .<= 75)

    Ψ_Bol = Ψ_EulBol .- Ψ_Eul

    WEul[expname] = -Float32.(Ψ_Eul[z_ref, ϕ_ref[1], :])
    WBol[expname] = -Float32.(Ψ_Bol[z_ref, ϕ_ref[1], :])
end


Hadvection = Dict(); diffusion = Dict(); GTF = Dict()
Vadvection = Dict()
Topadvection = Dict()
Botadvection = Dict()
exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
              "only_init", "only_kappa", "only_wind"]
temps = Dict()
for (i, expname) in enumerate(exps)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    temps[expname] = vars["θ"]
    uθ = vars["VθSouth"]
    ∇wθ = vars["wθBot"] .- vars["wθTop"]
    Hadvection[expname] = (uθ)
    Vadvection[expname] = (∇wθ)
    Botadvection[expname] = vars["wθBot"]
    Topadvection[expname] = -vars["wθTop"]

    diffusion[expname] = (vars["κxyθ"] .+ vars["κzθ"])
    GTF[expname] = (vars["GTH"])
end


fig, ax = plt.subplots(1, 2)
expname = "only_kappa"
ax[1].plot(-WBol[expname][1, :] .- WEul[expname][1, :])
ax[2].plot(Topadvection[expname])
fig

fig, ax = plt.subplots(1, 2)
ax[1].scatter(-WBol[expname][1, :] .- WEul[expname][1, :], Topadvection[expname])
fig

X = (WBol[expname][2, :] .+ WEul[expname][2, :]) ./ (sum(area .* PAC_msk))
b = Botadvection[expname]
a = pinv(X) * b


fig, ax = plt.subplots()
ax.plot(Topadvection[expname])
ax.plot(a .* X)
corr(Topadvection[expname], X)
fig