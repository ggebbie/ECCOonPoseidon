include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall
import PyPlot as plt 

include(srcdir("config_exp.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
 
ocean_mask = wet_pts(Γ)

region = "NPAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region, extent = "full")

ϕ_avg = zonal_average(ϕ, area .* PAC_msk)
ϕ_avg = ϕ_avg[isfinite.(ϕ_avg)]

#use same method to contruct the ϕ arr (just in case)
LC=LatitudeCircles(ϕ_avg,Γ); nl = length(LC)
lons = zeros(Float32, 360, nl); 
for ll in 1:nl
    IntegralPath = LC[ll]
    tabC=IntegralPath.tabC
    nC = size(tabC,1)
    for k=1:nC
        (a,i1,i2,w)=tabC[k,:]
        lons[k, ll] = round(Float32(λ.f[a][i1,i2] .* PAC_msk.f[a][i1,i2]), digits = 1)
    end
end
lons[lons .< 0.0] .+=360; lons[lons .== 0.0] .= NaN
lons_avg = zeros(360)
for ll = 1:360
    lons_avg[ll] = sum( ) / sum()
end
unique_lons = sort(unique(lons))[3:end]
lons_list = zeros(length(unique_lons)); 

include(srcdir("plot_and_dir_config.jl"))

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, mask)
    area = readarea(γ)
    area_mask = area .*mask
    LC=LatitudeCircles(ϕ_avg,Γ); nl = length(LC)

    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist)
    x = zeros(360, nl, nt)

    @time for tt = 1:nt
        println(tt)
        fname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        u, v, w = extract_eulerian_velocities(diagpath, expname, fname, γ)
        w_trsp = w[:, 43] .* area_mask;
        for ll in 1:nl
            IntegralPath = LC[ll]
            tabC=IntegralPath.tabC
            nC = size(tabC,1)
            for k=1:nC
                (a,i1,i2,w)=tabC[k,:]
                x[k, ll, tt] = w_trsp.f[a][i1,i2]
            end
        end

    end

    return x
end


expname = "only_sfc"
Whov_sfc = get_transports(diagpath, expname, γ, PAC_msk); 

expname = "iter0_bulkformula"
Whov_i0 = get_transports(diagpath, expname, γ, PAC_msk); 

expname = "only_kappa"
Whov_mix = get_transports(diagpath, expname, γ, PAC_msk); 

expname = "only_init"
Whov_init = get_transports(diagpath, expname, γ, PAC_msk); 

effects = [Whov_init .- Whov_i0, Whov_mix .- Whov_i0, Whov_sfc .- Whov_i0]
function meridional_sum(x)
    Wsum = zeros(length(unique_lons), 312)
    for (i, lon) in enumerate(unique_lons)
        lons_locs = findall(lon .== lons)
        for tt = 1:312
            x_vec = x[:, :, tt]
            Wsum[i, tt] .+= sum(x_vec[lons_locs])
        end
    end
    return Wsum
end

Whovs = [1e-6 .* meridional_sum(W) for W in effects]
Whovs[3]
pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass_filter(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/36, fs = 1),Butterworth(2))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, 312) .+ signal_mean
end

fig, ax = plt.subplots(1, 3,sharex = "col", figsize = (15.5, 5))
titles = ["INIT Effect", "MIXING Effect", "FORCING Effect"]
ax[1].set_ylabel("Time", fontweight = "bold")

for i in 1:3
    Whov_filt =  1 .* Whovs[i]
    Whov_filt[Whov_filt .== 0.0] .= NaN
    nl = size(Whov_filt, 1)    
    nan_mask = isfinite.( 1 .* Whov_filt[:, 100])

    ax[i].pcolormesh(unique_lons[nan_mask], tecco, Whov_filt[nan_mask, :]', cmap = cmo.balance, vmin = -.1, vmax = .1)
    Whov_filt[Whov_filt .== 0.0] .= NaN
    
    ax[i].set_title(titles[i])
    ax[i].set_xlabel("Longitude", fontweight = "bold")

    # ax[i, 2].set_xlim(-4, 4)
end
fig.subplots_adjust(wspace=0.4)
fig.savefig(plotsdir("native/sensitivity_exps/MeanEulHovmoller_Longitude_unfilt.png"), bbox_inches = "tight")
fig

fig, ax = plt.subplots(1, 3,sharex = "col", figsize = (17, 5))
titles = ["INIT Effect", "MIXING Effect", "FORCING Effect"]
ax[1].set_ylabel("Time", fontweight = "bold")

for i in 1:3
    Whov_filt =  1 .* Whovs[i]
    nl = size(Whov_filt, 1)    
    for ll = 1:nl
        Whov_filt[ll, :] .= low_pass_filter(Whov_filt[ll, :])
    end
    Whov_filt[Whov_filt .== 0.0] .= NaN
    nan_mask = isfinite.(sum(Whov_filt, dims = 2)[:, 1])
    ax[i].set_title(titles[i])
    ax[i].pcolormesh(unique_lons[nan_mask], tecco, Whov_filt[nan_mask, :]', cmap = cmo.balance, vmin = -.1, vmax = .1)
    ax[i].set_xlabel("Longitude", fontweight = "bold")
end
fig.subplots_adjust(wspace=0.3)
fig.savefig(plotsdir("native/sensitivity_exps/MeanEulHovmoller_Longitude_filt.png"), bbox_inches = "tight")

fig


[findall(wherelon), :]

# fig, ax = plt.subplots(3, 2, sharey = true,sharex = "col", figsize = (10, 20))
# titles = ["INIT Effect", "MIXING Effect", "FORCING Effect"]
# for i in 1:3
#     Whov_filt =  (Whovs[i][findall(wherelon), :])
#     nl = size(Whov_filt, 1)    
#     for ll = 1:nl
#         Whov_filt[ll, :] .= low_pass_filter(Whov_filt[ll, :])
#     end
#     ax[i, 1].set_title(titles[i])
#     ax[i, 1].pcolormesh(lons_avg, tecco, Whov_filt', cmap = cmo.balance, vmin = -.1, vmax = .1)

#     ax[i, 2].set_title(titles[i])
#     ax[i, 2].plot(sum(Whov_filt, dims = 1)[1, :], tecco)
#     # ax[i, 2].set_xlim(-4, 4)
# end
# fig.tight_layout()
# fig
