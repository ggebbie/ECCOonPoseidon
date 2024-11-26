include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings
import PyPlot as plt 
import NaNMath as nm
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

normalize(x, dims) = (x .- mean(x, dims = dims)) ./ std(x, dims = dims)
normalize(x) = (x[:] .- mean(x[:])) ./ std(x[:])

corr(x, y) = cov(normalize(x), normalize(y)) / (var(normalize(x)) * var(normalize(y)))
corr(x, y, dims) = cov(normalize(x, dims), normalize(y, dims); dims = dims) / (var(normalize(x); dims = dims) * var(normalize(y); dims = dims))

function regress_matrix(x)

    nt = length(x)

    # would be neat to return whatever type goes in
    E = zeros(Float32,nt,1)
    E[:, 1] .= x
    
    F = (E'*E)\E' # least squares estimator

    return F
end

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/3, fs = 1),Butterworth(2))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist);
    nϕ = length(ϕ_avg)
    curlτ = zeros(Float32, nϕ, nt);

    @time for tt = 1:nt
        println(tt)

        Tname = datafilelist_τ[tt]
        τx, τy = extract_ocnTAU(diagpath, expname , Tname, γ)
        τcurl = MeshArrays.curl(τx, τy, Γ)
        curlτ[:, tt] .= zonal_average(τcurl, PAC_msk .* area)[:]

    end

    return curlτ
end

τ_mean = Dict()
τ_mean["only_wind"]   = get_transports(diagpath, "only_wind", γ)
τ_mean["iter0_bulkformula"]   = get_transports(diagpath, "iter0_bulkformula", γ)

τ_mean["only_wind"] .-= τ_mean["iter0_bulkformula"]

vars =  ["only_wind",  "iter0_bulkformula"]
ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_PAC_only_init" *".jld2"))["ϕ_avg"]
ϕ_ref = findall( 23 .<= ϕ_avg .<= 40)
Ws = Dict(); ΨEul = Dict()
z_ref = findall( -3300 .<= -z[:].<= -2000)
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    ΨEul[expname] = 1 .* Ψ_exp
    
    tmp = -(Ψ_exp[1:end-1,  :, :]  .+ Ψ_exp[2:end,  :, :]) ./ 2
    tmp = tmp[:, ϕ_ref[1], :]
    tmp1 = tmp[15, :]
    tmp2 = sum(tmp[z_ref[1:end-1], :] .* ΔzF[z_ref[1:end-1]], dims = 1)[:] ./ sum(ΔzF[z_ref[1:end-1]])
    Ws[expname] = vcat(tmp1', tmp2')
    # Ws[expname] = -(Ψ_exp[[15, 40],  ϕ_ref[1], :] .- Ψ_exp[[15, 40],  ϕ_ref[end], :])

    fname = datadir("native/" * expname * "NPAC_THETA_budget_Wθ_eul_bol_2to3.jld2")
    vars = jldopen(fname)["dθ"]
    wθ_resid[expname] =  vars["wθBot"] .- vars["wθTop"]
    wθ_bol[expname] = vars["wθBot_bol"] .- vars["wθTop_bol"]
    wθ_eul[expname] = wθ_resid[expname] .- wθ_bol[expname]

end

Ws["only_wind"] .-= Ws["iter0_bulkformula"]
wθ_resid["only_wind"] .-= wθ_resid["iter0_bulkformula"]
ΨEul["only_wind"] .-= ΨEul["iter0_bulkformula"]

VAR = -1e-6 .* ΨEul["only_wind"][z_ref, ϕ_ref, :]
Y = reshape(1 .* VAR, :, 312)
Y = Y .- mean(Y, dims = 2)
Y[isnan.(Y)] .= 0.0

s = svd(Y; full = true)
PC1 = vec(s.U[:, 1]' * Y)

# function EOF1()

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); ϕ_avg_tau = 1 .* ϕ_avg
ϕ_ref = findall( 23 .<= ϕ_avg .<= 65)

fig,axs=plt.subplots(5, 1, figsize = (10, 15), sharex = true)

wind_color = exp_colors["only_wind"]
curl = mean(τ_mean["only_wind"][ϕ_ref, :], dims = 1)[:]
axs[1].plot(tecco,low_pass(curl), c= wind_color, linewidth = 3,label = "WIND Effect");
mult = 100 * 86400 * 365 * 10 
axs[1].set_ylabel("∇ × " * L"\mathbf{\tau_{wind}}" * " \n [N per m³]", fontweight = "bold")
axs[2].plot(tecco,mult * wθ_resid["only_wind"] , c= wind_color, linewidth = 3,label = "WIND Effect");
axs[2].set_ylabel( L"-\mathbf{∇(wθ)}" * " \n [cK per decade]", fontweight = "bold")
axs[3].plot(tecco, 1e-6 .* Ws["only_wind"][1, :],label = "WIND Effect", c= wind_color, linewidth = 3);
axs[4].plot(tecco, 1e-6 .* Ws["only_wind"][2, :],label =  "WIND Effect", c= wind_color, linewidth = 3);
axs[1].set_ylim(-2e-13, 2e-13)
axs[2].set_ylim(-2, 2); axs[2].invert_yaxis()
axs[3].set_ylabel(L"\mathbf{W |_{z = 150}}" * "\n [Sv]", fontweight = "bold")
axs[4].set_ylabel(L"\mathbf{W |_{z = 2500}}" * "\n [Sv]", fontweight = "bold")
[a.set_ylim(-4.5, 4.5) for a in axs[3:end-1]]
# [a.legend(frameon = false) for a in axs[2:end]]
axs[5].plot(tecco, -PC1[:], c= "k", linewidth = 3);
axs[5].set_ylabel("PC 1", fontweight = "bold")
# axs[5].set_ylim(-6, 6)
fig

axs[5].set_xlabel("time", fontweight = "bold")
fig.savefig(plotsdir("native/sensitivity_exps/WindEffectTau.png"), bbox_inches = "tight")

corr(normalize(PC1), normalize(curl))
corr(normalize(PC1), normalize(low_pass(curl[:])))
corr(normalize(PC1), normalize(PC1))
corr(normalize(PC1), normalize(Ws["only_wind"][1, :]))
corr(normalize(PC1), normalize(Ws["only_wind"][2, :]))
corr(normalize(PC1), -normalize(wθ_resid["only_wind"]))

nt  = length(curl)
freqs = 1 ./ reverse(FFTW.rfftfreq(nt)[2:end]) 
freqs ./= 12
Δt = 1; N = nt; T = N*Δt

spectral_density(x, N, T) = reverse((2*T * inv(N^2)) .* abs2.(rfft(x .- mean(x))[2:end] ))
rolling_spectral_density(x, N, T) = rollmean(spectral_density(x, N, T), 3)
perc_var_spectral_dens(x, N, T) = rolling_spectral_density(x, N, T) ./ sum(rolling_spectral_density(x, N, T))
fig,axs=plt.subplots(5, 1, figsize = (10, 15), sharex = true)

roll_freqs = rollmean(freqs, 3)
wind_color = exp_colors["only_wind"]
axs[1].plot(roll_freqs,perc_var_spectral_dens(curl[:], N, T), c= wind_color, linewidth = 3);
axs[1].set_ylabel("∇ × " * L"\mathbf{\tau_{wind}}'" * " \n [N per m³]", fontweight = "bold")
axs[2].plot(roll_freqs,perc_var_spectral_dens(-wθ_resid["only_wind"], N, T)  , c= wind_color, linewidth = 3);
axs[2].set_ylabel( L"-\mathbf{∇(wθ)}'" * " \n [N per m³]", fontweight = "bold")
axs[3].plot(roll_freqs, perc_var_spectral_dens(Ws["only_wind"][1, :], N, T) ,label = "z = 150 m", c= wind_color, linewidth = 3);
axs[4].plot(roll_freqs, perc_var_spectral_dens(Ws["only_wind"][2, :], N, T),label = "z = 2500 m", c= wind_color, linewidth = 3);
axs[5].plot(roll_freqs, perc_var_spectral_dens(PC1[:], N, T), c= "k", linewidth = 3);
axs[3].set_ylabel(L"\mathbf{W |_{z = 150}}" * "\n [Sv]", fontweight = "bold")
axs[4].set_ylabel(L"\mathbf{W |_{z = 2500}}" * "\n [Sv]", fontweight = "bold")
[ax.set_ylim(-0.01, 0.15) for ax in axs]
[ax.set_xscale("log") for ax in axs]
axs[5].set_ylabel("PC 1", fontweight = "bold")
axs[5].set_xlabel("Period [years]", fontweight = "bold")
fig.savefig(plotsdir("native/sensitivity_exps/WindEffectTauSpectra.png"), bbox_inches = "tight")
fig


objective = 1e-6 .*ΨEul["only_wind"]
objective = objective .- mean(objective, dims = 3)
F = (PC1'*PC1)\PC1' # least squares estimator

regress = reshape(objective , :, 312)* F' 
regress = reshape(regress, size(objective)[1:2]...)

reconstruction = zeros(size(objective))
for i in 1:size(reconstruction, 1), j in 1:size(reconstruction, 2)
    reconstruction[i, j, :] .= vec(regress[i, j] .* PC1)
end

denom = sum((objective .- mean(objective, dims = 3)).^2, dims = 3)[:, :, 1]
num = sum((reconstruction .- objective).^2, dims = 3)[:, :, 1]
Rsq = 100 .* (1 .- (num ./ denom))


objective = 1 .* τ_mean["only_wind"]
for i in 1:size(objective, 1)
    objective[i, :] .= low_pass(objective[i, :] )
end
objective = objective .- mean(objective, dims = 2)
regress = reshape(objective , :, 312)* F' 
regress = reshape(regress, size(objective)[1]...)
reconstruction = zeros(size(objective))
for i in 1:size(reconstruction, 1)
    reconstruction[i, :] .= vec(regress[i] .* PC1)
end

denom = sum((objective .- mean(objective, dims = 2)).^2, dims = 2)[:, 1]
num = sum((reconstruction .- objective).^2, dims = 2)[:, 1]
Rsq_tau = 100 .* (1 .- (num ./ denom))


read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_only_init.jld2"))
ϕ_avg = read_file["ϕ_avg"]

fig = plt.figure()

# Add the subplots
# The arguments for add_subplot are (number of rows, number of columns, index)
ax1 = fig.add_subplot(5, 1, 1)
ax2 = fig.add_subplot(5, 1, (2, 5), sharex=ax1)

vmax = round(nm.maximum(abs.(PC1_regress))) / 2
levels = 0:10:100
CM = ax2.contourf(ϕ_avg, z, Rsq, vmin = 0, vmax = 100, cmap = "Spectral_r", levels = levels)
ax2.invert_yaxis()
fig

ax1.plot(ϕ_avg_tau, Rsq_tau); ax1.set_ylim(0, 50)
fig.colorbar(CM, orientation = "horizontal")
fig