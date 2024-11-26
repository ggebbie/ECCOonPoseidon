include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall
import PyPlot as plt 
include(srcdir("config_exp.jl"))

include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

region = "NPAC"
NPAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

include(srcdir("plot_and_dir_config.jl"))

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); wherenan = (!isnan).(ϕ_avg .* 1)

vars =  ["only_kappa",  "iter0_bulkformula"]
ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_PAC_only_init" *".jld2"))["ϕ_avg"]
ϕ_ref = findall( 23 .<= ϕ_avg .<= 40)
Ψs = Dict()
z_ref = findall( -2700 .<= -z[:].<= -2300)[1]
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()
PMOC_max_lat = findall( -10 .<= ϕ_avg .<= -5)[3]
PMOC_max_d = findall( -4200  .<= -z[:].<= -3900)[end]
PMOC_max = Dict()
for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    # Ψs[expname] = Ψ_exp[z_ref, ϕ_ref, :][:]
    Ψs[expname] = -(Ψ_exp[[38, 40],  ϕ_ref[1], :])
    PMOC_max[expname] = -(Ψ_exp[PMOC_max_d,  PMOC_max_lat, :])
    fname = datadir("native/" * expname * "NPAC_THETA_budget_Wθ_eul_bol_" * suffix * ".jld2")
    vars = jldopen(fname)["dθ"]
    wθ_resid[expname] =  vars["wθBot"] .- vars["wθTop"]
    wθ_bol[expname] = vars["wθBot_bol"] .- vars["wθTop_bol"]
    wθ_eul[expname] = wθ_resid[expname] .- wθ_bol[expname]

end

PMOC_max["only_kappa"] .= PMOC_max["only_kappa"] .- PMOC_max["iter0_bulkformula"]
Ψs["only_kappa"] .= Ψs["only_kappa"] .- Ψs["iter0_bulkformula"]
wθ_resid["only_kappa"] .= wθ_resid["only_kappa"] .- wθ_resid["iter0_bulkformula"]

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/30, fs = 1),Butterworth(2))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end

function mult_regress(x, t)
    nt = length(t)
    # 2 = y-intercept and 2 trends
    nβ = 3
    # would be neat to return whatever type goes in
    E = zeros(Float32,nt,nβ)
    E[:,1] = ones(Float32,nt,1)
    E[:,2] .= t
    E[:,3] .= x
    F = (E'*E)\E' # least squares estimator
    return E,F
end


ϕ_avg = zonal_average(ϕ, PAC_msk .* area);
ϕ_ref = findall( 23 .<= ϕ_avg .<= 40)

fig,axs=plt.subplots(4, 1, figsize = (10, 15), sharex = true)
wind_color = exp_colors["only_wind"]
axs[1].plot(tecco, 1e-6 .* PMOC_max["only_kappa"], c= wind_color, linewidth = 3);
axs[1].set_ylabel("∇ × " * L"\mathbf{\tau_{wind}}'" * " \n [N per m³]", fontweight = "bold")
axs[2].plot(tecco,-wθ_resid["only_kappa"] , c= wind_color, linewidth = 3);
axs[2].set_ylabel("∇ × " * L"\mathbf{\tau_{wind}}'" * " \n [N per m³]", fontweight = "bold")
axs[3].plot(tecco, 1e-6 .* Ψs["only_kappa"][1, :],label = "z = 150 m", c= wind_color, linewidth = 3);
axs[4].plot(tecco, 1e-6 .* Ψs["only_kappa"][2, :],label = "z = 2500 m", c= wind_color, linewidth = 3);
# axs[5].plot(tecco, 1e-6 .* Ψs["only_wind"][3, :],label = "z = 2500 m", c= wind_color, linewidth = 3);
# axs[6].plot(tecco, 1e-6 .* Ψs["only_wind"][4, :],label = "z = 3900 m", c= wind_color, linewidth = 3);
[a.set_ylim(-0.7, 0.7) for a in axs[3:end]]
[a.set_ylabel(L"\mathbf{\psi_{wind}'}" * "\n [Sv]", fontweight = "bold") for a in axs[3:end]]
[a.legend(frameon = false) for a in axs[2:end]]

fig