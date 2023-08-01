
module OHC_helper
export patchvolume, calc_OHC, get_cell_volumes, standardize, 
       OHC_outputpath, do_FFT, calc_θ_bar, standard_error, conf_int,
       lin_reg, get_cell_depths, get_GH19, plot_patch, extract_3d_var, 
       compute_β, get_sfcfrc_fnames, vec, get_diff_arrs, wet_pts,
       rm_vals, PAC_mask, plot_div_0bf!, plot_ts!, div_ma,
       smush, levels_2_string, level_mean, plot_field!,
       reduce_dict, vert_int_ma,  get_geothermalheating, calc_UV_conv,
       plot_resids!, prune, central_diff, fwd_diff, fwd_mean, 
       perc_diff, resid, slice, element_gs, volume_mean, slice_mavec, 
       get_trend, ma_horiz_avg, remove_anomaly, ma_zonal_avg, ma_zonal_sum,
       nanmaximum, nanminimum, load_object_compress, load_and_calc_UV_conv,
       diff_ma_vec, level_timeseries!, total_level_change, 
       extract_θRbudget, extract_θHbudget, calc_UV_conv3D!, exch_UV_cs3D,
       filter_heat_budget, extract_sθ, depth_average, ThroughFlowDim, extract_meridionalΨ̄,
       getUVek!, extract_meridionalΨ2, extract_meridionalΨ̄_diff, extract_meridionalΨEK,
       UVtoUEVN3D, UVtoENTransport,extract_sX, sum_vertical, inv_ma,findlatlon, findmin, 
       tofaces, UVtoTrsp, velocity2center3D, extract_ocnTAU, wrap_λ, 
       extract_heatbudgetH, extract_heatbudgetR, LLCcropC360,
       calc_Wconv3D!, get_min_lat, velocity2centerfast, get_max_lat, 
       extract_velocities_and_bolus, within_lon, get_cs_and_sn, rotate_UV_native,
       get_msk, get_ϕ_max_min_mask, interpolate_to_vertical_faces!, region_mask, 
       sum_heat_flux_profile
       
export densityJMD95


export RMSE, RelDiff
using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, CodecZlib, DrWatson, FFTW, NetCDF,
    Printf, PyCall, RollingFunctions
import PyPlot as plt

import NaNMath as nm 
import Base: vec
import Statistics: mean, std
import Base: sum
@pyimport seaborn as sns
@pyimport pandas as pd

function mesharray_sum(X::MeshArrays.gcmarray{T, 2, Matrix{T}}, dims) where T<:Real
    if dims == 1
        nz = size(X, 2)
        tmp = zeros(T, nz)
        [tmp[k] = Float32(sum(cell_volumes[:, k])) for k=1:nz]
        return tmp

    elseif dims == 2
        tmp = T.(similar(X)); fill!(tmp, 0.0f0)
        for ijk in eachindex(X)
            tmp.f[ijk[1]] .+= X.f[ijk] #ijk[1] is the face
        end
        return tmp 
    end
end

mean(x::Vector, weights::Vector) = sum(x .* weights) / sum(weights)
std(x::Vector,  weights::Vector) = sqrt(inv(sum(weights)) * sum((x .- mean(x, weights)).^2))
RMSE(Δ) = sqrt(sum(Δ.^2))
RMSE(x, y) = sqrt(sum((x.-y).^2))
RelDiff(x, y) = (x-y)/ (abs(x) + abs(y))
notfinite(x) = isnan(x) || isinf(x)
central_diff(x, dt) = (x[3:end] .- x[1:end-2]) ./ (dt)
fwd_diff(x, dt) = (x[2:end] .- x[1:end-1]) ./ (dt)
fwd_mean(x)=(x[2:end].+x[1:end-1])./2
perc_diff(x, y) = 100 * (x - y)/y
resid(x, y) = x - y
slice(i) =pycall(pybuiltin("slice"), PyObject, i)
slice(i, j) =pycall(pybuiltin("slice"), PyObject, i, j)
element_gs(i,j, gs) = get(gs, (i,j))
remove_anomaly(data) = Dict(key => data[key] .- mean(data[key]) for key in keys(data))
nanmaximum(a) = maximum(filter(!isnan,a)) 
nanminimum(a) = minimum(filter(!isnan,a)) 
import Base: extrema
import Base: findmin

"""
inv_ma(ma::MeshArray)

Inverts any mesharray grid type object. Zero values will be filled with fillval. 
```
"""
function inv_ma(ma::MeshArray; fillval = Inf)
    inv_ma = deepcopy(ma)
    for a in eachindex(ma)
        ma.f[a][iszero.(ma.f[a])] .= fillval
        inv_ma.f[a] .= inv.(ma.f[a])
    end
    
    return inv_ma
end

function extrema(d::Dict)
    extremas = map(x-> nm.extrema(x), values(d))
    clims = (minimum(minimum.(extremas)), maximum(maximum.(extremas)))
    return clims
end


function div_ma(ma1::MeshArray, ma2::MeshArray; fillval = nothing)
    temp = similar(ma1)
    for ff in eachindex(ma1)
        temp[ff] .= ma1[ff] ./ ma2[ff]
        if !isnothing(fillval)
            temp[ff][notfinite.(temp[ff])] .= fillval
        end
    end
    return temp 
end

function volume_mean(x::MeshArray; weights::MeshArray)
    numerator = [0.0]
    denom =[0.0]
    for idx in eachindex(x)
        numerator .+= sum(x.f[idx] .* weights.f[idx])
        denom .+= sum(weights.f[idx])
    end
    return numerator[1]/denom[1]

end


"""
    function do_FFT(cell_area, cell_depths)
    get FFT of some signal 
# Arguments
- `y`: signal 
- `Fs`: sampling freqeuncy 
# Output
- `fft_freq`: Positive frequenceies which occur in y 
- 'X_mag': Amplitude of only Positive frequecies in y 
"""
function do_FFT(y, Fs)
    N = length(y)
    fft_freq = fftfreq(N, Fs)
    X_mag = abs.(fft(y)) ./ N
    #fft_freq is symmetric, return antisymmetric
    return (fft_freq[1:Int(N/2)], 2*X_mag[1:Int(N/2)])
end

"""
    function lin_reg(x, y_true)
    get slope for a set of data
# Arguments
- `x`: input data 
- `y_true`: output data

# Output
- `β`: estimated slope and intercept 
"""
function lin_reg(x, y_true)
    X = ones(length(x), 2)
    X[:, 2] .= x
    β = (X' * X ) \ X' * y_true
    return β
end


"""
    function standard_error(x, y_true, y_pred)
    get standard error for a set of data and predictions
# Arguments
- `x`: input data 
- `y_true`: output data
- `y_pred`: estimated output 

# Output
- `SE`: standard error of prediction 
"""

function standard_error(x, y_true, y_pred)
    denom = sum((x .- mean(x)).^2)
    num = sum((y_pred .- y_true).^2)
    n = length(x)
    df = n - 2
    SE = sqrt(num / (denom * df))
    return SE
end

"""
    function conf_int(Β, SE)
    get confidence interval at 95& confidence 
    for a slope parameter Β
# Arguments
- `Β`: slope parameter
- `SE`: standard error for slope parameter
# Output
- `H`: volume-integrated ocean heat content 
"""

function conf_int(Β, SE)
    #for 95% conf int
    t = 1.96
    my_conf_int = (-t * SE, t * SE)
    return my_conf_int

end

function get_GH19()
    GH19_file = datadir("oceanheatcontent_GH19.nc")
    GH19_time = ncread(GH19_file, "time")
    OHC_GH19_700 = ncread(GH19_file, "H700")
    OHC_GH19_mid = ncread(GH19_file, "Hmid")
    OHC_GH19_deep = ncread(GH19_file, "Hdeep")
    OHC_G19 = Dict("0-700" => OHC_GH19_700, "700-2000" =>OHC_GH19_mid, "2000-5500" => OHC_GH19_deep)

    return OHC_G19, GH19_time

end


function depth_weighted_mean(ma, γ::gcmgrid, 
    cell_depths::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    #This should be rewritten to match volume_weighted_avg
    num_dims = ndims(ma)
    if num_dims > 1
        var = MeshArray(γ,Float32)
        ztot = MeshArray(γ,Float32)
        fill!(var,0.0f0) 
        fill!(ztot,0.0f0) 
        max_lvl = size(ma, 2)
        for k in 1:max_lvl
            for ff in eachindex(var)
                var[ff] .+= ma[ff, k] .* cell_depths[ff, k]
                ztot[ff] .+= cell_depths[ff, k]
            end
        end
        weighted_mean = MeshArray(γ,Float32)
        fill!(weighted_mean,0.0f0) 
        for ff in eachindex(var)
            weighted_mean[ff] = var[ff] ./ ztot[ff]
        end
        
        return weighted_mean
    else 
        println("Only one level provided, no mean can be computed")
        return ma
    end
end

function get_sfcfrc_fnames(variation)
    eccodrive = "/batou/eccodrive/files/Version4/Release4/"
    if lowercase(variation) == "adjust"
        fileroot = "input_forcing/"
    elseif lowercase(variation) == "unadjust"
        fileroot = "other2/input_forcing_unadjusted/"
    else 
        return NaN
    end

    filelist_adj1 = searchdir(eccodrive*fileroot,"1992") 
    head = length("eccov4r4") + 2
    tail =  -length("1992") - 1
    keys = [key[head:end + tail] for key in filelist_adj1]
    file_list = Dict()
    for key in keys 
        file_list[key] = (eccodrive * fileroot) .* searchdir(eccodrive*fileroot,key) 
    end
    return file_list
end

function get_diff_arrs(γ)
    eccodrive = "/batou/eccodrive/files/Version4/Release4/"
    fileroot = "input_init/"
    filelist = searchdir(eccodrive*fileroot,"total_") 
    head = length("total_") + 1
    tail =  -length("_r009bit11.bin")
    keys = [key[head:end + tail] for key in filelist]
    nz = 50
    file_list = Dict()
    for key in keys 
        println(key)
        filelist = (eccodrive * fileroot) .* searchdir(eccodrive*fileroot,key) 
        adj_fname = filter(x -> occursin("data",x),filelist)[1]
        unadj_fname = filter(x -> occursin("bin",x),filelist)[1]
        xx = γ.read(adj_fname,MeshArray(γ,Float32,nz)) #first guess 
        total = γ.read(unadj_fname,MeshArray(γ,Float32,nz)) #adjustment
        file_list[key * "_unadj"] = total
        file_list[key * "_adj"] = (total .+ xx)
    end
    return file_list, keys
end

function rm_vals(dict, new_vals)
    new_dict = Dict()
    for (keys, values) in dict
        values ∈ new_vals && (new_dict[keys] = values)
    end
    return new_dict
end

function PAC_mask(Γ, basins, basin_list, ϕ, λ; region = "", 
    extent = "model-defined", include_bering = true)
    ocean_mask = wet_pts(Γ)
    basin_name="Pacific"
    # IndID = findall(basin_list.=="Indian")[1]
    basinID=findall(basin_list.==basin_name)[1]
    Okhotsk = findall(basin_list.=="Okhotsk Sea")[1]      
    basinID2 = findall(basin_list.== "Japan Sea")[1]   
    basinID3 = findall(basin_list.== "East China Sea")[1]   
    basinID4 = findall(basin_list.== "South China Sea")[1]   
    basinID5 =  findall(basin_list.== "Bering Sea")[1]   
    PAC_msk=similar(basins)
    bounds = [0.0, 0.0]
    if region == "NPAC"
        # bounds[1] = 51.0
        bounds[1] = 70.0
        bounds[2] = 23.0
    elseif region == "ESPAC"
        bounds[1] = 17.0
        bounds[2] = -35.0
    elseif region == "NPAC30"
        bounds[1] = 70.0
        bounds[2] = 30.0  
    elseif region == "SPAC"
        bounds[1] = 25.0
        bounds[2] = -40.0
    elseif region == "EPAC"
        bounds[1] = 25.0
        bounds[2] = -25.0
    elseif region == "PAC56"
        bounds[1] = 65.0
        bounds[2] = -56.0
    elseif region == "PAC"
        bounds[1] = 70.0
        bounds[2] = -40.0
    else 
        bounds[1] = 70.0
        bounds[2] = -40.0
    end
    full_extent = (extent == "full")

    for ff in 1:5
        above_SO = (bounds[2] .< ϕ[ff] .< bounds[1]) #make this go to bering strait
        marginal_seas = (basins[ff] .== Okhotsk) .+ (basins[ff] .== basinID2) .+ 
                        (basins[ff] .== basinID3) .+  (basins[ff] .== basinID4)
        Okhotsk_sea = (basins[ff] .== Okhotsk)
        SouthChinaSea = 0 
        BeringSea     = (basins[ff] .== basinID5) .* include_bering
        model_PAC     = (basins[ff] .== basinID) 
        temp = model_PAC .+ BeringSea .+ Okhotsk_sea .+ (full_extent .* marginal_seas)
        temp = ocean_mask[ff] .* temp .* above_SO 
        PAC_msk[ff] .= temp 
    end

    PAC_msk[findall(PAC_msk .== 0)] = 0f0
    PAC_msk[findall(PAC_msk .> 0)] = 1f0 

    return PAC_msk
end

function sum_vertical(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, γ) where T<:Real
    sum_temp = MeshArray(γ,Float32)
    fill!(sum_temp, 0.0f0)
    for ijk in eachindex(ma)
        sum_temp.f[ijk[1]] .+= ma.f[ijk] #ijk[1] is the face
    end
    return sum_temp 
end

function plot_field!(var, λ, ϕ, ax, color)
    projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
    ax.coastlines()
    ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    cf = Vector{Any}(undef ,1)
    # min_val = minimum(MeshArrays.mask(var ,Inf))
    max_val = maximum(MeshArrays.mask(var,-Inf))
    min_val = -max_val
    for ff in 1:length(var)
            cf[1] = ax.pcolormesh(λ[ff], ϕ[ff], var[ff],shading="nearest", 
            transform=projPC, rasterized = true, vmin = min_val, vmax = max_val,
            cmap = color)
    end
    return cf[1]
end

function vert_int_ma(ma, Δz, lvls, ts)
    intVar =Vector{MeshArrays.gcmarray{Float64, 1, Matrix{Float64}}}(undef, length(ts))
    for tt in 1:length(ts)
        intVar[tt] = 0.0e0 .* similar(ma[1][:, 1])
        intVar[tt] .= 0.0e0
        for lvl in lvls 
            tmp = ma[tt][:, lvl] .* Δz[:, lvl]
            intVar[tt] .+= tmp
        end
    end
    return intVar
end

function get_theta_init(γ)
    fname = "/batou/eccodrive/files/Version4/Release4/input_init/xx_theta.0000000129.data";
    θadjust = read_bin(fname,Float32,γ)
    
    #returns in terms of θ/s 
    return θadjust
end

function get_trend(var,tecco,F)
    β = [0.0]
    for tt in 1:length(tecco)
        β .+= F[2,tt] * var[tt]
    end
    return β[1]
end



function ma_zonal_sum(ma::MeshArrays.gcmarray{T, 1, Matrix{T}}) where T<:Real
    temp = zeros(270)
    for ff in eachindex(ma) 
	    if (ff==1) || (ff == 2)
            temp[1:270] .+= @views sum( ma[ff], dims = 1)[:] 
        elseif ff == 4 || ff == 5 
            temp[1:270] .+= @views sum( reverse(ma[ff], dims = 1), dims = 2)[:] 
        # elseif ff == 3
        #     temp[271:end] .+= @views sum( reverse(ma[ff], dims = 1), dims = 2)[:]
        end
    end
    return temp 
end

function ma_zonal_sum(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    nz = size(ma, 2)
    temp = zeros(nz, 270)
    for lvl in 1:nz
        temp[lvl, :] .= ma_zonal_sum(ma[:, lvl])
    end

    return temp 
end

function ma_zonal_avg(ma::MeshArrays.gcmarray{T, 2, Matrix{T}},
    weights ::MeshArrays.gcmarray{S, 2, Matrix{S}}) where {T<:Real, S<:Real}
    nlev = size(ma, 2)
    temp_num = zeros(nlev, 270)
    temp_denom = zeros(nlev, 270)

    for lvl in 1:nlev
        temp_num[lvl, :] .= ma_zonal_sum(ma[:, lvl] .* weights[:, lvl])
        temp_denom[lvl, :] .= ma_zonal_sum(weights[:, lvl])
    end
    return temp_num ./ temp_denom
end

function integratedθ(dθ::Vector{T}, θ₀) where T<:Real
    avgdθ = fwd_mean(dθ) #put dθ onto beginning of the month 
    dt = 2.628f6 # one month time step 
    return cumsum(vcat(θ₀, avgdθ .* dt))
end


function extract_sθ(expname::String,diagpath::Dict{String, String}, 
    γ::gcmgrid, fnameS::String, fnameθ::String, 
    inv_depths::MeshArrays.gcmarray{T, 1, Matrix{T}}) where T 
    θ = γ.read(diagpath[expname]*fnameθ,MeshArray(γ,Float32,50))
    ETAN = γ.read(diagpath[expname]*fnameS,MeshArray(γ,Float32,1))
    sθ = similar(θ)

    s1 = ETAN .* inv_depths
    s1 .+= 1
    sθ = θ .* s1
    return sθ
end

"""
    function extract_meridionalΨek
    extract the ekman stream function by reading multiple files
    approximation from Jayne 2001 
"""
function extract_meridionalΨEK(expname,diagpath, Γ, γ, Hinv, finv, mask)
    #exf_zflux_set1
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)

    ρ0 = MeshArray(γ,Float32); fill!(ρ0, 1027)
    hinv = MeshArray(γ,Float32); fill!(hinv, inv(10))
    surfthick = Γ.hFacC[:,1]; surfthick[findall(surfthick .== 0)] = Inf32
    surfcell_inv = 1 ./ surfthick ; surfcell_inv = Float32.(surfcell_inv)

    hinv = hinv .* surfcell_inv

    U=0.0*Γ.hFacW #varmeta not quite correct
    V=0.0*Γ.hFacS #varmeta not quite correct
    fill!(U, 0); fill!(V, 0)
    for tt=1:nt
        fname = datafilelist[tt];
        EXF = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,15))
        ETAN = EXF[:, 1]
        τx = EXF[:, 14]; τy = EXF[:, 15]
        s = ETAN .* Hinv; s .+= 1 # scaling factor 
        u = MeshArray(γ,Float32,50); v = MeshArray(γ,Float32,50);
        getUVek!(u, v, τx, τy, finv, Hinv, hinv);
        u = u .* s; v = v .* s; 
        for i in eachindex(U)
            U[i]=U[i]+(u[i]/nt)
            V[i]=V[i]+(v[i]/nt)
        end
    end
    U = U .* Γ.hFacW; V = V .* Γ.hFacS #also need to multiply by s
    (Utr,Vtr)=UVtoTransport(U,V,Γ)
    Utr = Utr .* mask; Vtr = Vtr .* mask;

    ov=Array{Float64,2}(undef,nl,nz)
    for z=1:nz
        UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
        [ov[l,z]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
    end
    ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
    ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
    #for Jayne's definition, we do not have to switch the value of ψ
    return ψ̄
end

"""
    function getUeVe!
    computes velocities resulting Ekman Transport and compensating return flow
    U, V : empty arrats to be filled
    τx, τy : (unrotated) windstress on sea surface 
    finv : 1/f
    Hinv : 1 / (total depth of column at (x, y))
    hinv : 1 / (depth of surface cell at (x, y))

"""
function getUVek!(U::MeshArrays.gcmarray{T, 2, Matrix{T}}, V::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    τx::MeshArrays.gcmarray{T, 1, Matrix{T}}, τy::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    finv::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    Hinv::MeshArrays.gcmarray{T, 1, Matrix{T}}, hinv::MeshArrays.gcmarray{T, 1, Matrix{T}}) where {T<:Real}
    for ff in 1:5
        τxff = τx.f[ff]; τyff = τy.f[ff]
        fρ0inv = finv[ff] .* inv(1030)
        Hinvff = Hinv.f[ff]; hinvff = hinv.f[ff]
        for k in 1:50
            V.f[ff, k] .= ((hinvff .* (k == 1)) .- Hinvff ) .* τxff
            V.f[ff, k] .=  -fρ0inv .* V.f[ff, k]
            U.f[ff, k] .= ((hinvff .* (k == 1)) .- Hinvff ).* τyff
            U.f[ff, k] .=  fρ0inv .* U.f[ff, k]
        end
    end
end

function UVtoUEVN3D(u::MeshArrays.gcmarray{T, 2, Matrix{T}},
    v::MeshArrays.gcmarray{T, 2, Matrix{T}},G::NamedTuple) where T<:Real
    nz = size(u, 2)

    cs=convert.(T,G.AngleCS)
    sn=convert.(T,G.AngleSN)

    (uex,vex)=exch_UV_cs3D(u,v);uC=similar(u); vC=similar(v)
    #interpolates to center of grid square
    @inbounds for iF=1:5, k = 1:nz
            tmp1=uex.f[iF, k][1:end-1,:]; tmp2=uex.f[iF, k][2:end,:] #averages between 
            uC.f[iF, k].=reshape(nanmean([tmp1[:] tmp2[:]],2),size(tmp1))
            tmp1=vex.f[iF, k][:,1:end-1]; tmp2=vex.f[iF, k][:,2:end]
            vC.f[iF, k].=reshape(nanmean([tmp1[:] tmp2[:]],2),size(tmp1))
    end

    return uC.*cs-vC.*sn, uC.*sn+vC.*cs
end

function UVtoENTransport(u::MeshArrays.gcmarray{T, 2, Matrix{T}},
    #this function does not work as expected. Does not return northward transport

    v::MeshArrays.gcmarray{T, 2, Matrix{T}},Γ::NamedTuple, γ::gcmgrid) where T<:Real
    Utr, Vtr = UVtoTransport(u, v, Γ); 
    Etr, Ntr = UVtoUEVN3D(Utr, Vtr, Γ)
    #remove mask
    for ijk in eachindex(v)
        Etr.f[ijk][(!isfinite).(Etr.f[ijk])] .= 0
        Ntr.f[ijk][(!isfinite).(Ntr.f[ijk])] .= 0
    end
    return Etr, Ntr
end

function depth_average(ds::MeshArray, Δh::MeshArray, H::MeshArray, γ)
    depth_avg = MeshArray(γ,Float32)
    fill!(depth_avg, 0.0)
    nz = size(ds, 2)
    for ff=1:5, k=1:nz
        depth_avg[ff] .+= (ds[ff, k] .* Δh[ff, k]) ./ H[ff]
    end
    return depth_avg
end


function velocity2center3D(u::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    v::MeshArrays.gcmarray{T, 2, Matrix{T}}, Γ) where T<:Real

        for iF in eachindex(u.f) 
            u.f[iF][Γ.hFacW.f[iF] .== 0] .= NaN
            v.f[iF][Γ.hFacS.f[iF] .== 0] .= NaN
        end
        (tmpU,tmpV)=exch_UV_cs3D(u,v); uC=T.(similar(u)); vC=T.(similar(v))
        #need to implement nanmean in this for correctness 
        for a in eachindex(u)
            (s1,s2)=size(u.f[a])
            tmpU1=view(tmpU.f[a],1:s1,1:s2)
            tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
            tmpV1=view(tmpV.f[a],1:s1,1:s2)
            tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
            uC.f[a] =(tmpU1+tmpU2) / 2 #breaks down at boundaries
            vC.f[a] =(tmpV1+tmpV2) / 2
        end
    return uC, vC
end

function wrap_λ(λ)
    λ_wrap = deepcopy(λ)
    [λ_wrap.f[ff][λ.f[ff] .< 0.0] .+= 360 for ff = 1:5]
    return λ_wrap
end

function pcolormesh_ma(ax, λ, ϕ, var, bounds, cmap, transformation, add_gridlines)
    CF = Any[]
    for ff = 1:5
        cf = ax.pcolormesh(λ[ff], ϕ[ff],  var[ff], 
                            vmin = -bounds, vmax = bounds,
                            transform=transformation, cmap = cmap, 
                            shading = "nearest", rasterized = false)
        push!(CF, cf)
    end
    if add_gridlines
        gl = ax.gridlines(crs=transformation, draw_labels=true, linewidth=2, 
        color="gray", alpha=0, linestyle="--")
        gl.top_labels = false; gl.right_labels = false
    end
    return CF[1]
end

function pcolormesh_ma_filt(ax, λ, ϕ, var, bounds, cmap, transformation, add_gridlines, filterfunc)
    CF = Any[]
    for ff = 1:5
        cf = ax.pcolormesh(λ[ff], ϕ[ff],  filterfunc(var[ff]), 
                            vmin = -bounds, vmax = bounds,
                            transform=transformation, cmap = cmap, 
                            shading = "nearest", rasterized = false)
        push!(CF, cf)
    end
    if add_gridlines
        gl = ax.gridlines(crs=transformation, draw_labels=true, linewidth=2, 
        color="gray", alpha=0, linestyle="--")
        gl.top_labels = false; gl.right_labels = false
    end
    return CF[1]
end

function get_min_lat(ϕ, mask)
    ϕ_mask_Inf = ϕ .* mask; 
    for ff = 1:5
        ϕ_mask_Inf.f[ff][ϕ_mask_Inf[ff] .== 0] .= Inf
    end
    return minimum(ϕ_mask_Inf)
end

function get_max_lat(ϕ, mask)
    ϕ_mask_Inf = ϕ .* mask; 
    [ϕ_mask_Inf.f[ff][ϕ_mask_Inf[ff] .== 0] .= -Inf for ff = 1:5]
    ϕ_mask_max = maximum(ϕ_mask_Inf)
    return ϕ_mask_max
end

function LLCcropC360(x, γ; modify_λ = false)
    xcrop = LLCcropC(x, γ)
    xwrap = vcat(xcrop[181:end, :], xcrop[1:180, :])
    if modify_λ
        xwrap[xwrap .< 0.0] .+= 360
    end
    return xwrap
end

function within_lon(λ, λ1, λ2)
    λ_mask = λ .> Inf
    for ff in eachindex(λ)
        λ_mask.f[ff] .= λ1 .< λ.f[ff] .< λ2
    end
    return λ_mask
end

function get_cs_and_sn(γ)
    cs = MeshArray(γ,Float32); sn = MeshArray(γ,Float32)
    fill!(cs, 0.0); fill!(sn, 0.0)
    cs.f[1] .= 1; cs.f[2] .= 1
    sn.f[4] .= -1; sn.f[5] .= -1
    return cs, sn
end

function rotate_UV_native(uvel::MeshArrays.gcmarray{T,N,Matrix{T}},
    vvel::MeshArrays.gcmarray{T,N,Matrix{T}},
    cs::MeshArrays.gcmarray{T,1,Matrix{T}}, sn::MeshArrays.gcmarray{T,1,Matrix{T}}) where T<:AbstractFloat where N
    zeros(T, )
    evel = T.(similar(uvel))
    nvel = T.(similar(vvel))
    for ff in eachindex(evel)
        evel.f[ff] .= uvel.f[ff] .* cs.f[ff[1]] .- vvel.f[ff] .*sn.f[ff[1]]
        nvel.f[ff] .= uvel.f[ff] .* sn.f[ff[1]] .+ vvel.f[ff] .*cs.f[ff[1]]
    end
    return evel,nvel
end



function get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)
    abs_dist(x, r) = abs(x) < r

    if region == "NPAC"
        PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
        region2 = "PAC56"
        PAC56_mask = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region2)
    
        ϕ_mask_min = Float32(get_min_lat(ϕ, PAC_msk)); ϕ_min_mask = ϕ .> -Inf
        ϕ_min_mask.f[3] .= 0.0
        [ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 23.8, 0.1) .* PAC56_mask.f[ff] for ff in 1:2]
        [ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 22.9, 0.1) .* PAC56_mask.f[ff] for ff in 4:5]

        ϕ_mask_max = Float32(get_max_lat(ϕ, PAC_msk)); ϕ_max_mask = ϕ .> -Inf; 
        # [ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 51.0, 0.1) .* PAC56_mask.f[ff] for ff in 1:2]
        # [ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 50.3, 0.1) .* PAC56_mask.f[ff] for ff in 4:5]
        return ϕ_min_mask, 0.0 .* ϕ_max_mask
    elseif region == "ESPAC"
        PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
        region2 = "PAC56"
        ϕ_min_mask.f[3] .= 0.0
        PAC56_mask = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region2)
        ϕ_mask_min = Float32(get_min_lat(ϕ, PAC_msk)); ϕ_min_mask = ϕ .> -Inf

        [ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- (-34.3), 0.1) .* PAC56_mask.f[ff] for ff in 1:2]
        [ϕ_min_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- (-35.1), 0.1) .* PAC56_mask.f[ff] for ff in 4:5]
        ϕ_mask_max = Float32(get_max_lat(ϕ, PAC_msk)); ϕ_max_mask = ϕ .> -Inf; 

        [ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 17.2, 0.1) .* PAC56_mask.f[ff] for ff in 1:2]
        [ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- 16.3, 0.1) .* PAC56_mask.f[ff] for ff in 4:5]
        return ϕ_min_mask, ϕ_max_mask
    end
    
end

function sum_heat_flux_profile(ds::MeshArray, ΔV)
    nz = size(ds, 2)
    vol_avg = zeros(Float32, nz)

    for ff=1:5, k=1:nz
        vol_avg[k] += Float32(sum(ds[ff, k])) ./ ΔV[k]
    end
    return vol_avg
end

end