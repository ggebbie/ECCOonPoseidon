#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    LaTeXStrings, PyCall, RollingFunctions
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));
sns.set_style("white")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

tecco = 1992+1/24:1/12:2018; tecco = Float32.(tecco)

region = "NPAC"; 
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"


nt = 312
labels = ["Initial 0", "Initial 129"]


fcycle = 1 # units: yr^{-1}
overtones= 10; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)
E,F = trend_matrices(tecco)


window_size = 12
window_pad = Int(window_size/2)
extract_interannual(x) =  vcat(fill(NaN32, window_pad), 
                            rollmean(x[:], window_size), fill(NaN32, window_pad-1))
mean_interannual(T, x) = mean(T, extract_interannual(x))

line_fit(x) = F * x[:]
compute_trend(x) = line_fit(x)[2]
round_trend(x) = string(round(27 * compute_trend(x), digits = 3))
line_fit(x, t) = line_fit(x)[1] .+ line_fit(x)[2] .*(t .- mean(t))
detrend(x, t) = x .- line_fit(x, t) .+ mean(x)

extract_interannual(x, t) = remove_seasonal(detrend(x, t), E, F) .+ mean(x)


fig, ax = plt.subplots(2, 3, figsize = (15, 8), sharey = true)
colors = ["red", "blue"]


for (i, var) in enumerate(["iter0_bulkformula", "iter129_bulkformula"])
    println("-----")
    println(var)
    fname = datadir("native/" * region * "_" * var * "_TRSP" * ".jld2")
    transports = load(fname)["transports"]
    Vin = 1e-6 .* transports.Vin
    Wtop = -1e-6 .* transports.Wtop
    Wbot = 1e-6 .* transports.Wbot
    axes = ax[i, :]

    # Vout = withoutBering .* transports_dict[var].Vout

    #plot original
    axes[1].plot(tecco, Vin, c = colors[i], alpha = 0.2); 
    axes[2].plot(tecco, Wtop, c = colors[i], alpha = 0.2); 
    axes[3].plot(tecco, Wbot, c = colors[i], alpha = 0.2); 
    println("Vin Final: ", Vin[end])

    axes[1].plot(tecco, extract_interannual(Vin), label = var * "\n trend: "  * round_trend(Vin), c = colors[i], alpha = 1, linewidth = 2)
    axes[2].plot(tecco, extract_interannual(Wtop), label = var * "\n trend: "  * round_trend(Wtop), c = colors[i], alpha = 1, linewidth = 2)
    axes[3].plot(tecco, extract_interannual(Wbot), label = var * "\n trend: "  * round_trend(Wbot), c = colors[i], alpha = 1, linewidth = 2)
    println()
    # axes[4].plot(tecco, Vin .+ Wtop .+ Wbot, label = var, c = colors[i], alpha = 1); 

    axes[1].set_title( L"V_{in}"); axes[2].set_title( L"W_{out}"); ; axes[3].set_title( L"W_{in}"); 
    # axes[4].set_title( "Convergence"); 
    fig.suptitle(region * " Volume Fluxes [z = 2 - 3 km]")
end

[a.set_ylim(-7, 7) for a in ax]
[a.axhline(0, linestyle = "--", alpha = 0.5, color = "k", zorder = 0) for a in ax]

[a.legend(frameon = false) for a in ax]
fig