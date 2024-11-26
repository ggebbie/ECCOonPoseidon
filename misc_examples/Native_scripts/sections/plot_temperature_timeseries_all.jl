tecco = 1992+1/24:1/12:2018; tecco = collect(tecco)
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)
nz = 50; nt = length(tecco)

using PyCall
@pyimport seaborn as sns;

sns.set_theme(context = "notebook", style = "ticks",
              palette = sns.color_palette("colorblind"));

fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5))
volume_weight_column(x, ΔV) =  sum(x .* ΔV, dims =1) ./ sum(ΔV)


region = "NPAC"
fnames(section, expname) = datadir(region * "_" * expname * "_" * section * "_THETA_levels" * ".jld2")

expname = "iter129_bulkformula"
section = "P01"
include("plot_temperature_timeseries_and_trend.jl")
include("plot_temperature_spectra.jl")
section = "P02"
include("plot_temperature_timeseries_and_trend.jl")
include("plot_temperature_spectra.jl")
section = "P16"
include("plot_temperature_timeseries_and_trend.jl")
include("plot_temperature_spectra.jl")


