
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
       get_msk, get_ϕ_max_min_mask, interpolate_to_vertical_faces!
       
export densityJMD95


export RMSE, RelDiff
using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, CodecZlib, DrWatson, FFTW, NetCDF,
    Printf, PyCall, RollingFunctions
import NaNMath as nm 
import CairoMakie as Mkie
import Base: vec
import Statistics: mean, std

@pyimport seaborn as sns
@pyimport pandas as pd
# colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "notebook", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));

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

"""
findlatlon(λ::MeshArray, ϕ::MeshArray, λ₀::T, ϕ₀T::)
# Arguments
- `λ`: latitude mesharray  
- `ϕ`: longitude mesharray
- `λ₀`: lat of interest
- `ϕ₀`: cell-volume (masked components should be zero) 
# Output
- ((face, level), (cartesian index of minimum within a face), (minimum value))


"""
function findlatlon(λ::MeshArray, ϕ::MeshArray, λ₀::T, ϕ₀::T) where T<:Real

    λ_dist = (λ .- λ₀); λ_dist = map.(x ->x^2, λ_dist);
    ϕ_dist = (ϕ .- ϕ₀); ϕ_dist = map.(x ->x^2, ϕ_dist);
    tot_dist = λ_dist .+ ϕ_dist
    return findmin(tot_dist)
end

"""
findmin(ma::MeshArray)

Finds the mininum of any MeshArray object. Iterates through each face and level, and finds their 
respective minimum index
Returns ((face, level), (cartesian index of minimum within a face), (minimum value))
```
"""

function findmin(ma::MeshArray)
    f_list = Any[Inf]
    i_list = Any[Inf]
    val_list = Any[Inf]
    for ijk in eachindex(ma)
        temp_val, temp_idx = findmin(ma[ijk])
        if val_list[1] > temp_val
            f_list[1] = ijk
            i_list[1] = temp_idx
            val_list[1] =  temp_val
        end
        
    end
    return (f_list[1], i_list[1], val_list[1])
end
"""
densityJMD95(Θ,Σ,Π,Π0) from MITgcmTools.jl/PhysicalOceanography.jl SeawaterDensity, From Gael Forget

Compute potential density (ρP), in situ density (ρI), and density
referenced to PREF (Π0 in decibars) from potential temperature (Θ in °C),
salinity (Σ in psu) and pressure (Π in decibars) according to the
UNESCO / Jackett & McDougall 1994 equation of state.

Credits: code based on a Matlab implementation by B. Ferron

Reference: https://www.jodc.go.jp/info/ioc_doc/UNESCO_tech/059832eb.pdf

Check value: ρI = `1041.83267kg/m^3` for Θ=`3°Celcius`, Σ=`35psu`, Π=`3000dbar`
```
(ρP,ρI,ρR) = SeaWaterDensity(3.,35.5,3000.)
isapprox(ρI,1041.83267, rtol=1e-6)
```
"""
function densityJMD95(Θ,Σ,Π,Π0=missing)

   #square root salinity
   sqrtΣ= sqrt.(Σ)
   #compute density pure water at atm pressure
   ZR1= ((((6.536332E-9*Θ .-1.120083E-6).*Θ .+1.001685E-4).*Θ
   .-9.095290E-3).*Θ .+6.793952E-2).*Θ .+999.842594
   #seawater density atm pressure
   ZR2= (((5.3875E-9*Θ .-8.2467E-7).*Θ .+7.6438E-5).*Θ
   .-4.0899E-3).*Θ .+0.824493
   ZR3= (-1.6546E-6*Θ .+1.0227E-4).*Θ .-5.72466E-3
   ZR4= 4.8314E-4

   #potential density (referenced to the surface)
   ρP= (ZR4*Σ + ZR3.*sqrtΣ + ZR2).*Σ + ZR1

   #add the compression terms
   ZE = (-3.508914E-8*Θ .-1.248266E-8).*Θ .-2.595994E-6
   ZBW= ( 1.296821E-6*Θ .-5.782165E-9).*Θ .+1.045941E-4
   ZB = ZBW + ZE .* Σ

   ZD = -2.042967E-2
   ZC = (-7.267926E-5*Θ .+2.598241E-3).*Θ .+0.1571896
   ZAW= ((5.939910E-6*Θ .+2.512549E-3).*Θ .-0.1028859).*Θ .-4.721788
   ZA = ( ZD*sqrtΣ + ZC).*Σ + ZAW

   ZB1= (-0.1909078*Θ .+7.390729).*Θ .-55.87545
   ZA1= ((2.326469E-3*Θ .+1.553190).*Θ .-65.00517).*Θ .+1044.077
   ZKW= (((-1.361629E-4*Θ .-1.852732E-2).*Θ .-30.41638).*Θ
   .+2098.925).*Θ .+190925.6
   ZK0= (ZB1.*sqrtΣ + ZA1).*Σ + ZKW

   #in situ density
   ρI = ρP ./ (1.0 .-Π./(ZK0-Π.*(ZA-Π.*ZB)))

   #density referenced to level Π0
   if !ismissing(Π0)
      ρR = ρP ./ (1.0 .-Π0./(ZK0-Π0.*(ZA-Π0.*ZB)))
   else
      ρR = ρP
   end

   return ρR
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
    function calc_OHC(diagpath, expname, datafilelist,γ,
        masked_volume, nz)  
    get weight for rectangle region.
# Arguments
- `diagpath`: variable of interest
- `expname`: experiment of interest 
- `datafilelist`: lat
- `γ`: ECCO Grid 
- `masked_volume`: cell-volume (masked components should be zero) 
- `nz` = number of depth levels in γ
# Output
- `H`: volume-integrated ocean heat content 
"""

function calc_OHC(diagpath, expname, datafilelist,γ,
    masked_volume, temp_grid, nz)  
    
    tt = 0
    nc = 2
    # c_p = 3991.86795711963
    cₚ = 3850.0
    ρ = 1029
    H = Array{Float64,2}(undef, nz, 0)

    for fname in datafilelist
        tt += 1
        x = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,nc*nz))
        θ = 273.15 .+ x[:, 1:nz] # load potential temperature
        # S = x[:, (nz+1):end] #load salinity 
        H_tt = zeros(nz)
        for ff in eachindex(θ) # ff is the cartesian index (face, depth level)
            depth_level = ff[2]            
            temp_grid[ff] .= cₚ .* ρ .* θ[ff] .* masked_volume[ff]
            H_tt[depth_level] += sum(temp_grid[ff])
        end
        H = hcat(H, H_tt)
    end
    return H
end
"""
    function extract_3d_var(γ, nc, diagpath, expname, datafilelist, nz)  
    extract_variable 
# Arguments
- `γ`: ECCO Grid 
- `nc`: number of variable stored in fname
- `diagpath`: variable of interest
- `expname`: experiment of interest 
- `datafilelist`: list of files pertaining to variable of interest
- `nz` = number of depth levels in γ

# Output
- `vart`: variable "varname" at each time (4-D object)
"""
function extract_3d_var(γ, nc, diagpath, expname, datafilelist, nz)  
    tt = 0
    if nz > 1
        vart = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, 0)
    else
        vart = Vector{MeshArrays.gcmarray{Float32, 1, Matrix{Float32}}}(undef, 0)
    end
    for fname in datafilelist
        tt += 1
        lbound = (nc-1)*nz
        rbound = (nc)*nz
        x = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,nc*nz))[:,lbound+1:rbound]
        push!(vart, x) # load potential temperature
    end
    return vart
end

"""
    function calc_θ_bar(diagpath, expname, datafilelist,γ,masked_volume,
        nz)  
    get weight for rectangle region.
# Arguments
- `diagpath`: variable of interest
- `expname`: experiment of interest 
- `datafilelist`: lat
- `γ`: ECCO Grid 
- `masked_volume`: cell-volume (masked components should be zero) 
- `nz` = number of depth levels in γ
- `level_volume` = volume of each depth level 

# Output
- `θ̄`: volume-weighted ocean temperature
"""

function calc_θ_bar(diagpath, expname, datafilelist,γ,masked_volume,
    nz, level_volume)  
    
    tt = 0
    θ̄ = Array{Float64,2}(undef, nz, 0)
    nc = 2

    for fname in datafilelist
        tt += 1
        # println("time index ",tt)

        x = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,nc*nz))
        θ = x[:, 1:nz] # load potential temperature
        θ̄_tt = zeros(nz)
        for ff in eachindex(θ) # ff is the cartesian index (face, depth level)
            depth_level = ff[2]            
            θ̄_tt[depth_level] += sum(θ[ff] .* masked_volume[ff]) / level_volume[depth_level]
        end
        θ̄ = hcat(θ̄, θ̄_tt)
    end
    return θ̄
end

"""
    function get_cell_volumes(cell_area, cell_depths)
    get the volume in area of interst 
# Arguments
- `cell_area`: s
- `cell_depths`: experiment of interest 
# Output
- `basin_volume`: volume-weighted ocean temperature
"""

function get_cell_volumes(cell_area, cell_depths)

    basin_volume = similar(cell_depths)

    # ff is the cartesian index (LLC90_face, depth level)
    for ff in eachindex(basin_volume)
        f_i = ff[1]
        for k in eachindex(basin_volume[ff])
          basin_volume[ff][k] = cell_area[f_i][k] * cell_depths[ff][k]
        end
    end

    return basin_volume
end

"""
    function get_cell_depths(basin_mask, Δz, hFacC)
    get true depth for each grid square in area of interest.
# Arguments
- `basin_mask`: mask with entries of interest equal to 1.0 else 0.0
- `Δz`: possible thickness of the calls at each level
- 'hFacC': percentage of possible thickness of each cell 
# Output
- `basin_volume`: volume-weighted ocean temperature
"""

function get_cell_depths(basin_mask, Δz, hFacC)

    cell_depths = similar(hFacC) .* 0.0f0
    cell_depths .= 0.0f0
    for ff in eachindex(cell_depths)
        #generate the depths of each tracer cell, done in main function to avoid copying 
        cell_depths.f[ff] .= basin_mask.f[ff[1]] .* (Δz[ff[2]] .* hFacC.f[ff])  #
    end

    return cell_depths
end

function OHC_outputpath(file_name, plotsdir, use_constant_density)
            
    final_plots_dir = joinpath(plotsdir, "OHC")
    if use_constant_density
        density_folder = "const_density"
    else 
        density_folder = "variable_density"
    end
    final_plots_dir = joinpath(final_plots_dir, density_folder)
    isdir(final_plots_dir) ? nothing : mkpath(final_plots_dir)
    println("Saving plots to: \n", final_plots_dir)
    return joinpath(final_plots_dir, file_name), final_plots_dir
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

"""
    function plot_patch(Γ, ocean_mask, ocean_name)
    plot a single patch of ocean 
# Arguments
- `Γ`: grid
- `ocean_mask`: boolean mask of patch 
- 'ocean_name': name of patch 
# Output
- `H`: volume-integrated ocean heat content 
"""

function plot_patch(Γ, patch_mask, patch_name)
    pth = MeshArrays.GRID_LLC90
    γ = GridSpec("LatLonCap",pth)
    basins=read(joinpath(pth,"v4_basin.bin"),MeshArray(γ,Float32))

    μ =Γ.hFacC[:,1]
    μ[findall(μ.>0.0)].=1.0
    μ[findall(μ.==0.0)].=NaN

    fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
	ax = Mkie.Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=patch_name* " (shown in red)")

	# basinID=findall(basin_list.==basin_name)[1]
	
	mx=maximum(basins)
	for ff in 1:length(Γ.RAC)
		col=μ[ff][:].*(patch_mask[ff][:])
		kk=findall(col.>0.0)
		!isempty(kk) ? Mkie.scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=:red,markersize=2.0) : nothing
		kk=findall((col.==0.0).*(!isnan).(μ[ff][:]))
        colors = (basins[ff][kk].*0.0) .+ 1.0 
		!isempty(kk) ? Mkie.scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=colors,
			colorrange=(0.0,mx),markersize=2.0,colormap=:lisbon) : nothing
	end
	Mkie.Colorbar(fig[1,2], colormap=:lisbon, colorrange=(0.0, mx), height = Mkie.Relative(0.65))
    save(patch_name*"_patch.png",fig)
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

function compute_β(var_ts, F::Matrix{Float32}, γ::gcmgrid; no_szn = false)
    nz = 50
    β = MeshArray(γ,Float32,nz)
    fill!(β,0.0f0) 
    if no_szn 
        tecco = collect(1992+1/24:1/12:2018)
        fcycle = 1 # units: yr^{-1}
        overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
        Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)
        var_ts_noszn = remove_seasonal(var_ts,Ecycle,Fcycle,γ) 
        for tt in 1:length(var_ts_noszn)
            var = var_ts_noszn[tt]
            β .+= F[2,tt] .* var
        end
        return β
    else 
        for tt in 1:length(var_ts)
            var = var_ts[tt]
            β .+= F[2,tt] .* var
        end
        return β
    end
end

"""
    function calc_OHC(diagpath, expname, datafilelist,γ,
        masked_volume, nz)  
    get weight for rectangle region.
# Arguments
- `diagpath`: variable of interest
- `expname`: experiment of interest 
- `datafilelist`: lat
- `γ`: ECCO Grid 
- `masked_volume`: cell-volume (masked components should be zero) 
- `nz` = number of depth levels in γ
# Output
- `H`: volume-integrated ocean heat content 
"""

function calc_OHC(θ, expname, datafilelist,γ,
    masked_volume, temp_grid, nz)  
    
    tt = 0
    nc = 2
    # c_p = 3991.86795711963
    cₚ = 3850.0
    ρ = 1029
    H = Array{Float64,2}(undef, nz, 0)

    for tt in 1:length(nz)
        tt += 1
        H_tt = zeros(nz)
        for ff in eachindex(θ) # ff is the cartesian index (face, depth level)
            depth_level = ff[2]            
            temp_grid[ff] .= cₚ .* ρ .* θ[ff] .* masked_volume[ff]
            H_tt[depth_level] += sum(temp_grid[ff])
        end
        H = hcat(H, H_tt)
    end
    return H
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


# function plot_ocn_field_histogram!(ma_field_dict, ma_weights_vec, keys, lvls, ax; bins = 10)
#     i = 0 
#     for key in keys
#             i += 1 
#             ma_field_vec = vec(ma_field_dict[key][:, lvls])
#             ax[i].hist(ma_field_vec, bins = bins, density = true, 
#             weights = ma_weights_vec)
#             ax[i].set_title(key)
#     end
# end



function vec(ma::MeshArrays.gcmarray)
    ma_vec = Array{Float32}[] 
    for fz in eachindex(ma)
        ma_vec = cat(ma_vec, vec(ma[fz]), dims = 1)
    end
    return ma_vec
end

function wet_pts(Γ)
    ocean_mask = similar(Γ.hFacC[:,1] )
    ocean_mask[findall(Γ.hFacC[:,1] .>0.0)].=1.0
    ocean_mask[findall(Γ.hFacC[:,1] .<=0.0)].=0.0
    return Float32.(ocean_mask)
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
        bounds[1] = 60.0
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
 
function plot_div_0bf!(var_dict, tecco, shortnames, ignore_list, ax; 
    ylabel = "", linestyle = "-", baseline = var_dict["iter0_bulkformula"])
    for (keys,values) in shortnames
        if values ∉ ignore_list
            println(values)
            diff = var_dict[keys] .- baseline
            ax.plot(tecco, diff, label = values*" - 0bf", ls = linestyle)
            MSD_str = @sprintf "MRD: %.2E" mean(RelDiff.(var_dict[keys], baseline))
            RMSE_str = @sprintf "RMSE: %.2E" RMSE(diff)
            println(RMSE_str)
            ax.set_ylabel(ylabel)
            ax.set_xlabel("Time")
            ax.legend()
        end
    end
end

function plot_ts!(var_dict, tecco, shortnames, ignore_list, ax; 
    ylabel = "", linestyle = "-", baseline = var_dict["iter0_bulkformula"][1],
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"])
    i = 0
    for (keys,values) in shortnames
        i+=1
        c = colors[i]
        if values ∉ ignore_list
            ax.plot(tecco, var_dict[keys] .- baseline, label = values, ls = linestyle,
                    color = c)
            ax.set_ylabel(ylabel)
            ax.set_xlabel("Time")
            ax.legend()
        end
    end
end 

function smush(ma, γ)
    sum_temp = MeshArray(γ,Float32)
    fill!(sum_temp, 0.0f0)
    for ijk in eachindex(ma)
        sum_temp.f[ijk[1]] .+= ma.f[ijk] #ijk[1] is the face
    end
    return sum_temp 
end

function sum_vertical(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, γ) where T<:Real
    sum_temp = MeshArray(γ,Float32)
    fill!(sum_temp, 0.0f0)
    for ijk in eachindex(ma)
        sum_temp.f[ijk[1]] .+= ma.f[ijk] #ijk[1] is the face
    end
    return sum_temp 
end

function depth_sum(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    temp = 0.0f0 .* similar(ma[:, 1])
    temp .= 0.0f0
    for ff in eachindex(ma)
        temp[ff[1]] .+= ma[ff]
    end
    return temp 
end


function levels_2_string(uplvl, botlvl)
    return string(abs(uplvl)) * " to " * string(abs(botlvl))
end
function level_mean(ma, depths)
    weight_ma = 0.0 .* similar(ma[:, 1])
    var_sum_ma = 0.0 .* similar(ma[:, 1])
    depth_sum_ma = 0.0 .* similar(ma[:, 1])
    weight_ma .= 0.0; var_sum_ma  .= 0.0;depth_sum_ma .= 0.0
    for ff in eachindex(weight_ma)
        for ij in eachindex(weight_ma[ff[1]])
            var_sum_ma[ff[1]][ij] = ma[ff][ij] .* depths[ff][ij]
            depth_sum_ma[ff[1]][ij] = depths[ff][ij]
        end
    end
    var_mean = var_sum_ma./depth_sum_ma
    return var_mean
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
function reduce_dict(names, ignore_list)
    temp = Dict() 
    for (keys,values) in names
        if values ∉ ignore_list
            temp[keys] = values
        end
    end
    return temp
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

function get_geothermalheating(Γ, γ)
    rho = 1029; cp = 3994; 
    ocean_mask = wet_pts(Γ);
    cell_depths = get_cell_depths(ocean_mask, Γ.DRF, Γ.hFacC); 

    bathy_mask = deepcopy(Γ.hFacC) #mask indicating bottom points
    for a in eachindex(bathy_mask)
        bathy_mask.f[a][bathy_mask.f[a] .> 0.0] .= 1.0
    end

    for ff in 1:5, k = 1:49
        bathy_mask.f[ff, k] .= bathy_mask.f[ff, k] .- bathy_mask.f[ff, k+1]
    end

    fname = "/batou/eccodrive/files/Version4/Release4/input_init/geothermalFlux.bin";
    GTF = read_bin(fname,Float32,γ)
    
    GTF3d = MeshArray(γ, Float32, 50)
    for a in eachindex(GTF3d)
        ρcpdz = cell_depths.f[a] .* (rho * cp); 
        ρcpdz[ρcpdz .== 0] .= Inf #keep empty cell empty
        GTF3d.f[a] .= (GTF.f[a[1]] .* bathy_mask.f[a]) ./ (ρcpdz)
    end

    #returns in terms of θ/s 
    return GTF3d
end

function get_theta_init(γ)
    fname = "/batou/eccodrive/files/Version4/Release4/input_init/xx_theta.0000000129.data";
    θadjust = read_bin(fname,Float32,γ)
    
    #returns in terms of θ/s 
    return θadjust
end

function calc_UV_conv(fldU::Vector{MeshArrays.gcmarray{T, 2, Matrix{T}}}, 
    fldV::Vector{MeshArrays.gcmarray{T, 2, Matrix{T}}}) where T<:Real
    temp = Vector{MeshArrays.gcmarray{T, 2, Matrix{T}}}(undef, length(fldU))
    for tt in 1:length(fldU)
        temp[tt] = calc_UV_conv(fldU[tt], fldV[tt])
    end
    return temp
end

function calc_UV_conv(fldU::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
                        fldV::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
                        
    nlvls = size(fldU, 2)
    temp = similar(fldU)

    for lvl in 1:nlvls
        temp.f[:, lvl] .= MeshArrays.convergence(fldU[:, lvl], fldV[:, lvl]).f
    end
    return temp
end


function prune(x)
    idx = sortperm(x)
    remove_amt = Int(round(length(idx) * 0.025))
    y = copy(x)
    y[idx[end-remove_amt:end]] .= 0.0
    y[idx[1:remove_amt]] .= 0.0 
    return y
end

function slice_mavec(mavec::Vector{MeshArrays.gcmarray{T, 2, Matrix{T}}}, 
                    lvls::Vector{Int64}, γ::gcmgrid) where T <: Real
    nlevels = length(lvls)
    nlevels > 1 ? numdims = 2 : numdims = 1
    nt = length(mavec)
    temp_vec = Vector{MeshArrays.gcmarray{Float32, numdims, Matrix{Float32}}}(undef, nt)
    maxlevels=50
    nlevels = length(lvls)
    nt = length(mavec)
    @inbounds for tt in 1:nt
        if lvls[end]< (maxlevels+1)
            temp_vec[tt] = mavec[tt][:, lvls]
        else
            tmp = MeshArray(γ,Float32,nlevels)
            clip_ma = mavec[tt].f[:, lvls[1:end-1]]
            for ff in 1:5
                for k in 1:nlevels-1
                    tmp.f[ff, k] .=  clip_ma[ff, k]
                end
                tmp[ff, nlevels] .= 0
            end
            temp_vec[tt] = tmp
        end 
    end
    return temp_vec
end

function get_trend(var,tecco,F)
    β = [0.0]
    for tt in 1:length(tecco)
        β .+= F[2,tt] * var[tt]
    end
    return β[1]
end

function ma_horiz_avg(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    weights::MeshArrays.gcmarray{S, 2, Matrix{S}}) where {T<:Real, S<:Real}
    nlev = size(ma, 2)
    numerator = @views [sum(ma[:, lvl] .* weights[:, lvl]) for lvl in 1:nlev]
    denominator = @views [sum(weights[:, lvl]) for lvl in 1:nlev]
    return Float32.(numerator) ./ Float32.(denominator)
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

function load_object_compress(filename::String)
    jldopen(filename, "r"; compress= true) do file
        all_keys = keys(file)
        length(all_keys) == 0 && throw(ArgumentError("File $filename does not contain any object"))
        length(all_keys) > 1 && throw(ArgumentError("File $filename contains more than one object. Use `load` or `@load` instead"))
        file[all_keys[1]] #Uses HDF5 functionality of treating the file like a dict
    end
end

function load_and_calc_UV_conv(xpath::String, ypath::String, lvls::Vector{Int64}, γ::gcmgrid)
    var_exp = load_object_compress(xpath);
    x = slice_mavec(var_exp, lvls, γ)
    var_exp = load_object_compress(ypath);
    y = slice_mavec(var_exp, lvls, γ)
    println("Computing horizontal convergences...")
    convH = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, length(x))
    for tt in 1:length(fldU)
        convH[tt] = calc_UV_conv(x[tt], y[tt])
    end
    return convH
end

function load_and_calc_UV_conv(xpath::String, ypath::String)
    x = load_object_compress(xpath);
    y = load_object_compress(ypath);
    println("Computing horizontal convergences...")
    nt = length(x)
    convH = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, nt)
    for tt in 1:nt
        convH[tt] = calc_UV_conv(x[tt], y[tt])
    end
    return convH
end

function diff_ma_vec(var_ts::Vector{MeshArrays.gcmarray{T, 2, Matrix{T}}}, 
    lvls::Vector{Int64}, γ::gcmgrid, nt::Int64) where T<:Real
    temp_vec = Vector{MeshArrays.gcmarray{T, 2, Matrix{T}}}(undef, nt)
    ma_vec = slice_mavec(var_ts, lvls, γ); 
    ma_vecp1 = slice_mavec(var_ts, lvls.+1, γ);
    temp_vec = @views ma_vec .- ma_vecp1
    return temp_vec
end

function diff_ma_vec(var_ts::Vector{MeshArrays.gcmarray{T, 2, Matrix{T}}}, 
    γ::gcmgrid) where T<:Real
    ma_vec = var_ts
    ma_vecp1 = slice_mavec(var_ts, collect(2:51), γ);
    temp_vec = @views ma_vec .- ma_vecp1
    return temp_vec
end

function xz_avg(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    weights::MeshArrays.gcmarray{S, 2, Matrix{S}}) where {T<:Real,S<:Real} 
    numer = ma_zonal_sum(depth_sum(ma .* weights))
    denom = ma_zonal_sum(depth_sum(weights))
    returval = numer ./ denom
    return returval
end
function integratedθ(dθ::Vector{T}, θ₀) where T<:Real
    avgdθ = fwd_mean(dθ) #put dθ onto beginning of the month 
    dt = 2.628f6 # one month time step 
    return cumsum(vcat(θ₀, avgdθ .* dt))
end


function level_timeseries(θC, θT, crop_vols, msk, lvls)
    level_reconstruction = Dict()
    nt = length(θC["AdvH"])
    nlvls = length(lvls)
    depths_arrays = zeros(nlvls, nt)

    for key in collect(keys(θC))
        level_reconstruction[key] = zeros(nlvls, nt)
        for tt in 1:nt, i in 1:nlvls
            depths_arrays[i, tt] = sum(θC[key][tt][:, i] .* msk) / sum(crop_vols[:, i])
            if (key == "AdvR") || (key == "DiffZ")
                depths_arrays[i, tt] = -depths_arrays[i, tt]
            end
        end
        θ₀ = [OHC_helper.volume_mean(θT[1][:, lvl]; weights =  crop_vols[:, lvl]) for lvl in 1:nlvls]
        for i in 1:nlvls
            level_reconstruction[key][i, :] .= integratedθ(depths_arrays[i, :], θ₀[i])
        end
    end
    return level_reconstruction
end

function level_timeseries!(level_reconstruction, θC, crop_vols, lvls, msk, tt)
    nlvls = length(lvls)

    for key in collect(keys(level_reconstruction))
        for i in 1:nlvls
            level_reconstruction[key][i, tt] = sum(θC[key][:, i] .* msk) / sum(crop_vols[:, i])
            if (key == "AdvR") ||  (key == "DiffR")
                level_reconstruction[key][i, tt] = -level_reconstruction[key][i, tt]
            end
        end

    end
end


function total_level_change(level_reconstruction::Dict, tecco, F)
    total_change = Dict()
    nlvls = size(level_reconstruction["AdvH"], 1)
    level_reconstruction["total"] = level_reconstruction["AdvH"] .+ level_reconstruction["DiffH"]+ 
                                level_reconstruction["DiffZ"] .+ level_reconstruction["AdvR"]
    for key in keys(level_reconstruction)
        total_change[key] = zeros(nlvls)
        for lvl in 1:nlvls
            diff_ = level_reconstruction[key][lvl, end] - level_reconstruction[key][lvl,1]
            total_change[key][lvl] = diff_ * 26
        end

    end
    return total_change
end

function total_level_change(level_reconstruction::Dict, θ1, crop_vols, lvls)
    nlev = length(lvls)
    θ₀ = [volume_mean(θ1[:, i]; weights =  crop_vols[:, i]) for i in 1:nlev]
    for key in collect(keys(level_reconstruction)), i in 1:nlev
        level_reconstruction[key][i, :] .= integratedθ(level_reconstruction[key][i, :], θ₀[i])
    end

    total_change = Dict()
    level_reconstruction["total"] = level_reconstruction["AdvH"] .+ level_reconstruction["DiffH"]+ 
                                level_reconstruction["DiffR"] .+ level_reconstruction["AdvR"]
    for key in keys(level_reconstruction)
        total_change[key] = zeros(nlev)
        for lvl in 1:nlev
            diff_ = level_reconstruction[key][lvl, end] - level_reconstruction[key][lvl,1]
            total_change[key][lvl] = diff_ * 26
        end

    end
    return total_change
end

function exch_UV_cs3D(fldU::MeshArrays.gcmarray{T, 2, Matrix{T}},
    fldV::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    fillval=0f0
    #step 1

    s=size.(fldU[:, 1].f)
    nz = size(fldU, 2)
    nf=fldU.grid.nFaces
    s=vcat(s,s[3]) #always has 5 faces in LLC90
    tp=fldU.grid.class

    FLDU=similar(fldU)
    FLDV=similar(fldV)
    (ovfW,ovfE,ovfS,ovfN,evfW,evfE,evfS,evfN)=MeshArrays.exch_cs_viewfunctions();
    for lvl=1:nz, a=1:nf
        FLDU.f[a, lvl] = fill(fillval,s[a][1]+1,s[a][2]);
        FLDV.f[a, lvl] = fill(fillval,s[a][1],s[a][2]+1);
        FLDU.f[a, lvl][1:s[a][1],1:s[a][2]] .= fldU.f[a, lvl];
        FLDV.f[a, lvl][1:s[a][1],1:s[a][2]] .= fldV.f[a, lvl];

        (jW, jE, jS, jN)=MeshArrays.exch_cs_target(s[a],1)
        (aW,aE,aS,aN,iW,iE,iS,iN)=MeshArrays.exch_cs_sources(a,s,1)

        if (!iseven)(a)
            (aE <= nf) && (FLDU.f[a, lvl][jE[1].-1,jE[2].-1].=ovfE(fldU.f[aE, lvl],iE[1],iE[2]))
            (aN <= nf) && (FLDV.f[a, lvl][jN[1].-1,jN[2].-1].=ovfN(fldU.f[aN, lvl],iN[1],iN[2]))
        else
            (aE <= nf) && (FLDU.f[a, lvl][jE[1].-1,jE[2].-1].=evfE(fldV.f[aE, lvl],iE[1],iE[2]))
            (aN <= nf) && (FLDV.f[a, lvl][jN[1].-1,jN[2].-1].=evfN(fldV.f[aN, lvl],iN[1],iN[2]))
        end
    end
    return FLDU,FLDV

end

function calc_Wconv3D!(wFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, CONV::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    
    for ff=1:5, k=1:49
        @inbounds CONV.f[ff, k] .= wFLD.f[ff, k+1] .- wFLD.f[ff, k] #in - out 
    end

    for ff=1:5
        @inbounds CONV.f[ff, 50] .= -wFLD.f[ff, 50] #in - out 
    end
end

function calc_UV_conv3D!(uFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
vFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, CONV::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    tmpU, tmpV = exch_UV_cs3D(uFLD,vFLD)
    tmpU = T.(tmpU); tmpV = T.(tmpV); 
    for a in eachindex(uFLD.f)
        (s1,s2)=size(uFLD.f[a])
        @inbounds tmpU1=view(tmpU.f[a],1:s1,1:s2)
        @inbounds tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
        @inbounds tmpV1=view(tmpV.f[a],1:s1,1:s2)
        @inbounds tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
        @inbounds CONV.f[a] .= tmpU1-tmpU2+tmpV1-tmpV2
    end
end

"""
function extract_sst34
extract by reading multiple files
"""
function extract_θHbudget(expname::String,diagpath::Dict{String, String}, 
γ::gcmgrid, fnameH::String)
    dθλ = γ.read(diagpath[expname]*fnameH,MeshArray(γ,Float32,200))
    κθx, κθy = MeshArray(γ,Float32,50),MeshArray(γ,Float32,50) 
    uθx, uθy = MeshArray(γ,Float32,50),MeshArray(γ,Float32,50)
    d = Dict(name => MeshArray(γ,Float32,50) for name in ["AdvH", "DiffH"])
    @inbounds κθx.f .= dθλ.f[ :, 1:50]; @inbounds κθy.f .= dθλ.f[:, 51:100]
    @inbounds uθx.f .= dθλ.f[:, 101:150]; @inbounds uθy.f .= dθλ.f[:, 151:200]
    calc_UV_conv3D!(κθx, κθy, d["DiffH"]); 
    calc_UV_conv3D!(uθx, uθy, d["AdvH"]);
    return d
end

function extract_θRbudget(expname::String,diagpath::Dict{String, String}, 
γ::gcmgrid, fnameR::String)
    d = Dict(name => MeshArray(γ,Float32,50) for name in ["AdvR", "DiffR"])
    dθr = γ.read(diagpath[expname]*fnameR,MeshArray(γ,Float32,150))
    @inbounds d["AdvR"].f .= dθr.f[:, 1:50];
    @inbounds d["DiffR"].f .= @views dθr.f[ :, 51:100] .+ dθr.f[ :, 101:150];
    @inbounds d["AdvR"].f[ :, 1:49] .= @views d["AdvR"].f[:, 1:49] .- d["AdvR"].f[:, 2:50]
    @inbounds d["DiffR"].f[ :, 1:49] .= @views d["DiffR"].f[:, 1:49] .- d["DiffR"].f[:, 2:50]
    return d
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
    ThroughFlow(VectorField,IntegralPath,Γ::NamedTuple)

Compute transport through an integration path. Taken from Forget JuliaClimate/transport 
notebook. Removed if-statements and related code for a speedup. 
"""
function ThroughFlowDim(U::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
                        V::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
                        IntegralPath) where {T<:Real}

    #Note: vertical intergration is not always wanted; left for user to do outside
    nd=ndims(U)
    n=fill(1,4)
    tmp=size(U)
    n[1:nd].=tmp[1:nd]
    trsp=Array{Float64}(undef,1,n[3],n[4])
    for i3=1:n[3]
        tabW=IntegralPath.tabW
        tabS=IntegralPath.tabS
        for i4=1:n[4]
            trsp[1,i3,i4]=0.0
            for k=1:size(tabW,1)
                @inbounds (a,i1,i2,w)=tabW[k,:]
                @inbounds u=U.f[a,i3,i4][i1,i2]
                @inbounds trsp[1,i3,i4]=trsp[1,i3,i4]+w*u
            end
            for k=1:size(tabS,1)
                @inbounds (a,i1,i2,w)=tabS[k,:]
                @inbounds v=V.f[a,i3,i4][i1,i2]
                @inbounds trsp[1,i3,i4]=trsp[1,i3,i4]+w*v
            end
        end
    end
    return trsp[1]
end


"""
    function extract_meridionalΨ
    extract the total stream function by reading multiple files
"""
function extract_meridionalΨ̄(expname,diagpath, Γ, γ, mask)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    #[1:36, 312-36+1:312] comparing first and last years
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)
    U=0.0*Γ.hFacW #varmeta not quite correct
    V=0.0*Γ.hFacS #varmeta not quite correct
    fill!(U, 0); fill!(V, 0)
    for tt=1:nt
        fname = datafilelist[tt]
        UV = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,100))
        u = UV[:, 1:50]
        v = UV[:, 51:100]
        for i in eachindex(U)
            U[i]=U[i]+(u[i]/nt)
            V[i]=V[i]+(v[i]/nt)
        end
    end

    (Utr,Vtr)=UVtoTransport(U,V,Γ)
    Utr = Utr .* mask; Vtr = Vtr .* mask;
    ov=Array{Float64,2}(undef,nl,nz)
    for z=1:nz
        UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
        [ov[l,z]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
    end
    ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
    ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
    return ψ̄
end

"""
    function extract_meridionalΨ
    extract the total stream function by reading multiple files
"""
function extract_meridionalΨ̄wBolus(expname,diagpath, Γ, γ, mask)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    #[1:36, 312-36+1:312] comparing first and last years
    tecco = 1992+1/24:1/12:2018
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)
    U=0.0*Γ.hFacW #varmeta not quite correct
    V=0.0*Γ.hFacS #varmeta not quite correct
    fill!(U, 0); fill!(V, 0)
    mskC, mskW, mskS = get_msk(Γ)

    for tt=1:nt
        fname = datafilelist[tt]
        u, v, w, Ub, Vb, Wb = extract_velocities_and_bolus(diagpath, expname, 
                                fname, γ, Γ, mskC, mskW, mskS)
        for i in eachindex(U)
            U.f[i] .+= ((u.f[i] .+ Ub.f[i])./nt)
            V.f[i] .+= ((v.f[i] .+ Vb.f[i])./nt)
        end
    end

    (Utr,Vtr)=UVtoTransport(U,V,Γ)
    Utr = Utr .* mask; Vtr = Vtr .* mask;
    ov=Array{Float64,2}(undef,nl,nz)
    for z=1:nz
        UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
        [ov[l,z]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
    end
    ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
    ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
    return ψ̄
end


"""
    function extract_meridionalΨ
    extract the total stream function by reading multiple files
"""
function extract_meridionalΨ̄wBolus_timeseries(expname,diagpath, Γ, γ, mask)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    #[1:36, 312-36+1:312] comparing first and last years
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)
    U=0.0*Γ.hFacW #varmeta not quite correct
    V=0.0*Γ.hFacS #varmeta not quite correct
    fill!(U, 0); fill!(V, 0)
    ψ = zeros(nl,nz, nt)
    mskC, mskW, mskS = get_msk(Γ)

    for tt=1:nt
        fname = datafilelist[tt]
        u, v, w, Ub, Vb, Wb = extract_velocities_and_bolus(diagpath, expname, 
                                fname, γ, Γ, mskC, mskW, mskS)
        u .+= Ub
        v .+= Vb
        (Utr,Vtr)=UVtoTransport(u,v,Γ)
        Utr = Utr .* mask; Vtr = Vtr .* mask;
        ov=Array{Float64,2}(undef,nl,nz)
        for z=1:nz
            UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
            [ov[l,z]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
        end
        ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
        ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
        ψ[:, :, tt] .= ψ̄
    end

    return ψ

end

"""
    function extract_meridionalΨ
    extract the stream function trend in units Sv/century
"""
function extract_meridionalΨ̄timeseries(expname,diagpath, Γ, γ, mask)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    #[1:36, 312-36+1:312] comparing first and last years
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)
    U=0.0*Γ.hFacW #varmeta not quite correct
    V=0.0*Γ.hFacS #varmeta not quite correct
    fill!(U, 0); fill!(V, 0)
    ψ = zeros(nl,nz, nt)
    for tt=1:nt
        println(tt)
        fname = datafilelist[tt]
        UV = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,100))
        u = UV[:, 1:50]
        v = UV[:, 51:100]
        (Utr,Vtr)=UVtoTransport(u,v,Γ)
        Utr = Utr .* mask; Vtr = Vtr .* mask;
        ov=Array{Float64,2}(undef,nl,nz)
        for z=1:nz
            UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
            [ov[l,z]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
        end
        ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
        ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
        ψ[:, :, tt] .= ψ̄
    end

    return ψ
end


"""
    function extract_meridionalΨ
    extract the total stream function by reading multiple files
    only look at 3 year averages at beginning and end of run 
"""
function extract_meridionalΨ̄_diff(expname,diagpath, Γ, γ, mask)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)
    Ψs=Array{Float64,3}(undef,2, nl,nz)
    # bounds_array = [1:numyrs*12, 312-(numyrs*12)+1:312]
    U=0.0*Γ.hFacW #varmeta not quite correct
    V=0.0*Γ.hFacS #varmeta not quite correct
    fill!(U, 0); fill!(V, 0)
    for tt in 1:nt
        fname = datafilelist[tt]
        UV = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,100))
        u = UV[:, 1:50]
        v = UV[:, 51:100]
        for i in eachindex(U)
            U[i]=U[i]+(u[i]/nt)
            V[i]=V[i]+(v[i]/nt)
        end
    end
    (Utr,Vtr)=UVtoTransport(U,V,Γ)
    Utr = Utr .* mask; Vtr = Vtr .* mask;
    ov=Array{Float64,2}(undef,nl,nz)
    for z=1:nz
        UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
        [ov[l,z]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
    end
    ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
    ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
    Ψs[j, :, :] = ψ̄
    return Ψs
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

function get_msk(Γ)
    mskW = 0 .*Γ.hFacW; mskS = 0 .*Γ.hFacS; 
    mskC = 0 .*Γ.hFacC;
    @inbounds for ff in eachindex(mskC)
        mskW.f[ff][Γ.hFacW.f[ff] .> 0] .= 1
        mskS.f[ff][Γ.hFacS.f[ff] .> 0] .= 1;
        mskC.f[ff][Γ.hFacC.f[ff] .> 0] .= 1;
    end
    return Float32.(mskC), Float32.(mskW), Float32.(mskS)
end

function calc_bolus(GM_PsiX::MeshArrays.gcmarray{T, 2, Matrix{T}},
    GM_PsiY::MeshArrays.gcmarray{T, 2, Matrix{T}}, Γ)  where {T<:Real}
    nr=length(Γ.RC);
    println(sum(GM_PsiX))
    mskW = 0 .*Γ.hFacW; mskS = 0 .*Γ.hFacS; 
    mskC = 0 .*Γ.hFacC;
    @inbounds for ff in eachindex(mskC)
        mskW.f[ff][Γ.hFacW.f[ff] .> 0] .= 1
        mskS.f[ff][Γ.hFacS.f[ff] .> 0] .= 1;
        mskC.f[ff][Γ.hFacC.f[ff] .> 0] .= 1;
        GM_PsiX.f[ff][(!isfinite).(GM_PsiX.f[ff])] .=0; #check if this is actually needed later
        GM_PsiY.f[ff][(!isfinite).(GM_PsiY.f[ff])] .=0;
    end

    bolusU=0*Γ.hFacW;
    bolusV=0*Γ.hFacS;
    for k=1:nr-1;
        @inbounds bolusU[:,k] =(GM_PsiX[:,k+1].-GM_PsiX[:,k]) ./Γ.DRF[k];
        @inbounds bolusV[:,k] =(GM_PsiY[:,k+1].-GM_PsiY[:,k]) ./Γ.DRF[k];
    end;
    bolusU[:, nr] = -GM_PsiX[:,nr]./Γ.DRF[nr];
    bolusV[:, nr] = -GM_PsiY[:,nr]./Γ.DRF[nr];

    #and its vertical part
    #   (seems correct, leading to 0 divergence)
    #tmpx and tmpy are the BOLUS transports
    tmp_x = T.(similar(GM_PsiX)); fill!(tmp_x, 0.0)
    tmp_y = T.(similar(GM_PsiX)); fill!(tmp_y, 0.0)
    bolusW = T.(similar(GM_PsiX)); fill!(bolusW, 0.0)

    for a in eachindex(tmp_x)
        tmp_x.f[a] .= GM_PsiX.f[a] .* Γ.DYG.f[a[1]]
        tmp_y.f[a] .= GM_PsiY.f[a] .* Γ.DXG.f[a[1]]
    end

    tmpU, tmpV = exch_UV_cs3D(tmp_x,tmp_y)
    @inbounds for a in eachindex(tmp_x)
        (s1,s2)=size(tmp_x.f[a])
        tmpU1=view(tmpU.f[a],1:s1,1:s2)
        tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
        tmpV1=view(tmpV.f[a],1:s1,1:s2)
        tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
        bolusW.f[a] .= tmpU1-tmpU2+tmpV1-tmpV2
        bolusW.f[a] .= bolusW.f[a] ./ Γ.RAC.f[a[1]]

        bolusU.f[a] .= bolusU.f[a] .*mskW.f[a];
        bolusV.f[a] .= bolusV.f[a] .*mskS.f[a];
        bolusW.f[a] .= bolusW.f[a] .*mskC.f[a];
    end

    return bolusU, bolusV, bolusW
end


function calc_bolus(GM_PsiX::MeshArrays.gcmarray{T, 2, Matrix{T}},
    GM_PsiY::MeshArrays.gcmarray{T, 2, Matrix{T}}, Γ, 
    mskC::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    mskW::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    mskS::MeshArrays.gcmarray{T, 2, Matrix{T}})  where {T<:Real}
    nr=length(Γ.RC);

    bolusU = T.(similar(Γ.hFacW));
    bolusV = T.(similar(Γ.hFacS));

    @inbounds for k=1:nr-1;
        bolusU[:,k] =(GM_PsiX[:,k+1].-GM_PsiX[:,k]) ./ Γ.DRF[k];
        bolusV[:,k] =(GM_PsiY[:,k+1].-GM_PsiY[:,k]) ./ Γ.DRF[k];
    end
    @inbounds bolusU.f[:, nr] = -GM_PsiX.f[:,nr]./Γ.DRF[nr];
    @inbounds bolusV.f[:, nr] = -GM_PsiY.f[:,nr]./Γ.DRF[nr];

    #and its vertical part
    #   (seems correct, leading to 0 divergence)
    #tmpx and tmpy are the BOLUS transports
    tmp_x = T.(similar(GM_PsiX));
    tmp_y = T.(similar(GM_PsiX));
    bolusW = T.(similar(GM_PsiX));

    @inbounds for a in eachindex(tmp_x)
        tmp_x.f[a] .= GM_PsiX.f[a] .* Γ.DYG.f[a[1]]
        tmp_y.f[a] .= GM_PsiY.f[a] .* Γ.DXG.f[a[1]]
    end

    calc_UV_conv3D!(tmp_x, tmp_y, bolusW)
    @inbounds for a in eachindex(tmp_x)
        #negative for divergence instead of convergence
        #need to rescale bolus W
        bolusW.f[a] .= -1 .* bolusW.f[a] ./ Γ.RAC.f[a[1]] 
        bolusU.f[a] .= bolusU.f[a] .* mskW.f[a];
        bolusV.f[a] .= bolusV.f[a] .* mskS.f[a];
        bolusW.f[a] .= bolusW.f[a] .* mskC.f[a];
    end
    # tmpU, tmpV = exch_UV_cs3D(tmp_x,tmp_y)
    # @inbounds for a in eachindex(tmp_x)
    #     (s1,s2)=size(tmp_x.f[a])
    #     tmpU1=view(tmpU.f[a],1:s1,1:s2)
    #     tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
    #     tmpV1=view(tmpV.f[a],1:s1,1:s2)
    #     tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
    #     bolusW.f[a] .= tmpU2 - tmpU1+tmpV2 -tmpV1
    #     bolusW.f[a] .= bolusW.f[a] ./ Γ.RAC.f[a[1]]

    #     bolusU.f[a] .= bolusU.f[a] .* mskW.f[a];
    #     bolusV.f[a] .= bolusV.f[a] .* mskS.f[a];
    #     bolusW.f[a] .= bolusW.f[a] .* mskC.f[a];
    # end

    return bolusU, bolusV, bolusW
end

nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,dims=y)

function interpolate_to_vertical_faces!(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    ma_interp::MeshArrays.gcmarray{T, 2, Matrix{T}})  where T<:Real
    for k = 1:49 
        ma_interp[:, k+1] = (ma[:, k] .+ ma[:, k+1]) ./2 #linear interpolate to faces
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

function TransportThruLat(U::MeshArrays.gcmarray{T, 2, Matrix{T}}, V::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    lat, Γ::NamedTuple) where T<:Real
    (Utr,Vtr)=UVtoTransport(U,V,Γ)
    LC=LatitudeCircles(lat,Γ)
    nz = size(U, 2) 
    ov=Array{Float64,2}(undef,1,nz)
    @time for z=1:nz
        UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
        ov[1,z]=ThroughFlow(UV,LC[1],Γ)
    end
    ψ̄ = -cumsum(reverse(ov,dims=2),dims=2)' ; ψ̄[ψ̄.==0.0].=NaN
    return ψ̄
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


function depth_average(ds::MeshArray, Δh::MeshArray, H::MeshArray, γ)
    depth_avg = MeshArray(γ,Float32)
    fill!(depth_avg, 0.0)
    nz = size(ds, 2)
    for ff=1:5, k=1:nz
        depth_avg[ff] .+= (ds[ff, k] .* Δh[ff, k]) ./ H[ff]
    end
    return depth_avg
end

function tofaces(uFLD::MeshArrays.gcmarray{T, 1, Matrix{T}},
    vFLD::MeshArrays.gcmarray{T, 1, Matrix{T}}, Γ) where T<:Real
    Ugrid = T.(similar(uFLD))
    Vgrid = T.(similar(uFLD))
    (tmpU,tmpV)=exch_UV(uFLD,vFLD)

    for a in 1:5
        (s1,s2)=size(uFLD.f[a])
        tmpU1=view(tmpU.f[a],1:s1,1:s2)
        tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
        tmpV1=view(tmpV.f[a],1:s1,1:s2)
        tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
        Ugrid.f[a] .=(tmpU1+tmpU2) ./ 2
        Vgrid.f[a] .=(tmpV1+tmpV2) ./ 2
    end

return Ugrid, Vgrid
end

function tofaces(uFLD::MeshArrays.gcmarray{T, 2, Matrix{T}},
    vFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, Γ) where T<:Real
    Ugrid = T.(similar(uFLD))
    Vgrid = T.(similar(uFLD))
    DXG = T.(Γ.DXG); DYG = T.(Γ.DXG); 
    (DXG,DYG)=exch_UV(DXG,DYG)

    (tmpU,tmpV)=exch_UV_cs3D(uFLD,vFLD)

    for a in eachindex(uFLD)
        (s1,s2)=size(uFLD.f[a])
        tmpU1=view(tmpU.f[a],1:s1,1:s2)
        tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
        tmpV1=view(tmpV.f[a],1:s1,1:s2)
        tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
        #setup weights
        tmpDXG1 =view(DXG.f[a[1]],1:s1,1:s2)
        tmpDXG2 =view(DXG.f[a[1]],2:s1+1,1:s2)
        tmpDYG1 =view(DYG.f[a[1]],1:s1,1:s2)
        tmpDYG2 =view(DYG.f[a[1]],1:s1,2:s2+1)

        Ugrid.f[a] .= ((tmpU1 .* tmpDXG1) + (tmpU2 .* tmpDXG2)) ./ (tmpDXG1 .+ tmpDXG2)
        Vgrid.f[a] .= ((tmpV1 .* tmpDYG1) + (tmpV2 .* tmpDYG2)) ./ (tmpDYG1 .+ tmpDYG2)
    end

return Ugrid, Vgrid
end

function UVtoTrsp(U::MeshArray,V::MeshArray,G::NamedTuple)
    uTr=deepcopy(U)
    vTr=deepcopy(V)
    for i in eachindex(U)
        uTr.f[i] .= G.DRF[i[2]]*U.f[i].*G.DYG.f[i[1]]
        vTr.f[i] .= G.DRF[i[2]]*V.f[i].*G.DXG.f[i[1]]
    end
    return uTr,vTr
end

nanmeans(x) = mean(filter(!isnan,x))
nanmeans(x,y) = mapslices(nanmeans,x,dims=y)

function velocity2centers(u::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    v::MeshArrays.gcmarray{T, 1, Matrix{T}}, Γ) where T<:Real

        for iF in eachindex(u.f) 
            u.f[iF][Γ.hFacW.f[iF] .== 0] .= NaN
            v.f[iF][Γ.hFacS.f[iF] .== 0] .= NaN
        end
        (tmpU,tmpV)=exch_UV(u,v); uC=T.(similar(u)); vC=T.(similar(v))

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

function velocity2centerfast(u::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    v::MeshArrays.gcmarray{T, 1, Matrix{T}}, Γ) where T<:Real

        # for iF in eachindex(u.f) 
        #     u.f[iF][Γ.hFacW.f[iF] .== 0] .= NaN
        #     v.f[iF][Γ.hFacS.f[iF] .== 0] .= NaN
        # end
        (tmpU,tmpV)=exch_UV(u,v); uC=T.(similar(u)); vC=T.(similar(v))
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

function extract_ocnTAU(diagpath, expname , τdatafilelist, tt, γ)
    if expname == "seasonalclimatology"
        @time EXF = γ.read(diagpath[expname]*τdatafilelist[tt],MeshArray(γ,Float32,10))
        τx = EXF[:, 9]; τy = EXF[:, 10]; 
        return τx, τy
    else
        @time EXF = γ.read(diagpath[expname]*τdatafilelist[tt],MeshArray(γ,Float32,15))
        τx = EXF[:, 14]; τy = EXF[:, 15]; 
        return τx, τy
    end
end

function extract_velocities(diagpath, expname , fnameuvw, γ)
    UV = γ.read(diagpath[expname]*fnameuvw,MeshArray(γ,Float32,150))
    u = UV[:, 1:50]; v = UV[:, 51:100]; w = UV[:, 101:150]
    return u, v, w
end

function extract_velocities_and_bolus(
    diagpath, expname , fnameuvw, γ, Γ, 
    mskC::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    mskW::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    mskS::MeshArrays.gcmarray{T, 2, Matrix{T}})  where T<:Real
    UV = γ.read(diagpath[expname]*fnameuvw,MeshArray(γ,Float32,250))
    u = UV[:, 1:50]; v = UV[:, 51:100]; w = UV[:, 101:150]
    GM_PsiX = UV[:, 151:200]; GM_PsiY = UV[:, 201:250]
    Ub, Vb, Wb = calc_bolus(GM_PsiX,GM_PsiY, Γ, mskC, mskW, mskS)
    return u, v, w, Ub, Vb, Wb
end

function extract_heatbudgetH(diagpath, expname, fnameHθ, γ)
    dθλ = γ.read(diagpath[expname]*fnameHθ,MeshArray(γ,Float32,200))

    @inbounds κUθ = dθλ[:, 1:50]; 
    @inbounds κVθ = dθλ[:, 51:100]; 
    @inbounds Uθ = dθλ[:, 101:150]; 
    @inbounds Vθ = dθλ[:, 151:200]; 

    return κUθ, κVθ, Uθ, Vθ
end

function extract_heatbudgetR(diagpath, expname, fnameRθ, γ)
    #vertical convergences
    dθr = γ.read(diagpath[expname]*fnameRθ,MeshArray(γ,Float32,150))
    @inbounds wθ = dθr[:, 1:50]; 
    @inbounds κzθ = dθr[ :, 51:100] .+ dθr[ :, 101:150];

    return κzθ, wθ
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

    evel = T.(similar(uvel))
    nvel = T.(similar(vvel))
    for ff in eachindex(evel)
        evel.f[ff] .= uvel.f[ff] .* cs.f[ff[1]] .- vvel.f[ff] .*sn.f[ff[1]]
        nvel.f[ff] .= uvel.f[ff] .* sn.f[ff[1]] .+ vvel.f[ff] .*cs.f[ff[1]]
    end
    return evel,nvel
end

function get_msk(Γ)
    mskW = Float32.(Γ.hFacW); mskS = Float32.(Γ.hFacW); 
    mskC = Float32.(Γ.hFacW);
    @inbounds for ff in eachindex(mskC)
        mskW.f[ff][Γ.hFacW.f[ff] .> 0] .= 1.0;
        mskS.f[ff][Γ.hFacS.f[ff] .> 0] .= 1.0;
        mskC.f[ff][Γ.hFacC.f[ff] .> 0] .= 1.0;
    end
    return mskC, mskW, mskS
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
end
# @load datadir("ECCO_vars/ADVx_TH_iter129_bulkformula.jld2")