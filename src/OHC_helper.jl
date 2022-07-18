
module OHC_helper
export patchvolume, calc_OHC, get_cell_volumes, standardize, 
       OHC_outputpath, do_FFT, calc_θ_bar, standard_error, conf_int,
       lin_reg, get_cell_depths, get_GH19, plot_patch, extract_3d_var, 
       compute_β, get_sfcfrc_fnames, vec, get_diff_arrs, wet_pts,
       rm_vals, PAC_mask, plot_div_0bf!, plot_ts!, div_ma,
       smush, levels_2_string, level_mean, plot_field!,
       reduce_dict, vert_int_ma,  get_geothermalheating, calc_UV_conv,
       plot_resids!

export RMSE, RelDiff
using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, FFTW, NetCDF,
    Printf

import CairoMakie as Mkie
import Base: vec
import Statistics: mean, std
mean(x::Vector, weights::Vector) = sum(x .* weights) / sum(weights)
std(x::Vector,  weights::Vector) = sqrt(inv(sum(weights)) * sum((x .- mean(x, weights)).^2))
RMSE(Δ) = sqrt(sum(Δ.^2))
RMSE(x, y) = sqrt(sum((x.-y).^2))
RelDiff(x, y) = (x-y)/ (abs(x) + abs(y))
notfinite(x) = isnan(x) || isinf(x)

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

function mean(x::MeshArray, weights::MeshArray)
    numerator = [0.0]
    denom =[0.0]
    if size(x) == size(weights)
        for idx in eachindex(x)
            numerator .+= sum(x[idx] .* weights[idx])
            denom .+= sum(weights[idx])
        end
        return numerator[1]/denom[1]
    else 
        error("x and weights are not the same size, cannot do weighted average")
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
        x = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,nc*nz))
        lbound = (nc-1)*nz
        rbound = (nc)*nz
        push!(vart, x[:,lbound+1:rbound]) # load potential temperature
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
        cell_depths[ff] .= basin_mask[ff[1]] .* (Δz[ff[2]] .* hFacC[ff])  #
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
    cell_depths::MeshArrays.gcmarray{Float64, 2, Matrix{Float64}})
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
    return ocean_mask
end

function rm_vals(dict, new_vals)
    new_dict = Dict()
    for (keys, values) in dict
        values ∈ new_vals && (new_dict[keys] = values)
    end
    return new_dict
end

function PAC_mask(Γ, basins, basin_list, ϕ; region = "")
    ocean_mask = wet_pts(Γ)
    basin_name="Pacific"
    basinID=findall(basin_list.==basin_name)[1]
    PAC_msk=similar(basins)
    bounds = [0.0, 0.0]
    if region == "NPAC"
        bounds[1] = 60.0
        bounds[2] = 20.0
    elseif region == "SPAC"
        bounds[1] = 20.0
        bounds[2] = -40.0
    elseif region == "EPAC"
        bounds[1] = 20.0
        bounds[2] = -20.0
    elseif region == "PAC56"
        bounds[1] = 50.0
        bounds[2] = -56.0
    else 
        bounds[1] = 50.0
        bounds[2] = -40.0
    end
    println(region)
    for ff in 1:5
        above_SO = (bounds[2] .< ϕ[ff] .< bounds[1]) #removes southern ocean 
        PAC_msk[ff] .= ocean_mask[ff].*(basins[ff].==basinID) .* above_SO
    end
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
    ylabel = "", linestyle = "-", baseline = var_dict["iter0_bulkformula"][1])
    for (keys,values) in shortnames
        if values ∉ ignore_list
            ax.plot(tecco, var_dict[keys] .- baseline, label = values, ls = linestyle)
            ax.set_ylabel(ylabel)
            ax.set_xlabel("Time")
            ax.legend()
        end
    end
end

function smush(ma)
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

function get_geothermalheating(γ)
    fname = "/batou/eccodrive/files/Version4/Release4/input_init/geothermalFlux.bin"
    GTF = γ.read(fname,MeshArray(γ,Float32))
    return GTF
end


"""
    function calc_UV_conv(fldU,fldV)
        - fldU 
        - fldV
    This function calculates the divergences within a control volume 
    for the llc90 grid. Some faces must be rotated to to fit with the others, 
    so this cannot be looped through easily using the standard ECCO output. 
    Divergences for the cap cannot be computed at this time. 
     


"""
# function calc_UV_conv(fldU::MeshArrays.gcmarray{Float32, 1, Matrix{Float32}},
#     fldV::MeshArrays.gcmarray{Float32, 1, Matrix{Float32}})
#     ∇ₕ = similar(fldV)
#     ∇ₕ .= 0.0

#     ∇ₕ[1][1:end-1, 1:end-1] .= (fldV[1][1:end-1, 2:end] -  fldV[1][1:end-1, 1:end-1]) + 
#                                  (fldU[1][2:end, 1:end-1] -  fldU[1][1:end-1, 1:end-1]);
#     ∇ₕ[1][end, :] .+= fldU[1][end, :] - fldU[2][1, :]
#     ∇ₕ[1][:, end] .+= reverse(fldV[3][1, :]) .- fldV[1][:, end]   

#     ∇ₕ[2][1:end-1, 1:end-1] .= (fldV[2][1:end-1, 2:end] -  fldV[2][1:end-1, 1:end-1]) + 
#                                  (fldU[2][2:end, 1:end-1] -  fldU[2][1:end-1, 1:end-1])
#     ∇ₕ[2][end, :] .+= reverse(fldU[4][:, 1]) .- fldU[2][end, :]
#     #2 does not connect to face 3(the poles)

#     ∇ₕ[4][2:end, 1:end-1] .= (fldV[4][1:end-1, 1:end-1] - fldV[4][2:end, 1:end-1]) + 
#                                (fldU[4][2:end, 2:end] -  fldU[4][2:end, 1:end-1]) 
#     ∇ₕ[4][:, end] .+= fldU[5][:, 1]- fldU[4][:, end]
#     ∇ₕ[4][1, :] .+= fldV[3][end, :] - fldV[4][1, :]   

#     ∇ₕ[5][2:end, 1:end-1] .= (fldV[5][1:end-1, 1:end-1] - fldV[5][2:end, 1:end-1]) + 
#                                (fldU[5][2:end, 2:end] -  fldU[5][2:end, 1:end-1]) 
#     ∇ₕ[5][:, end] .+= reverse(fldU[1][1, :]) - fldU[5][:, end]   
#     ∇ₕ[5][1, :] .+= reverse(fldV[3][:, end]) - fldV[5][1, :]   

#     return ∇ₕ 
# end

function calc_UV_conv(fldU::Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}, 
    fldV::Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}})
    temp = similar(fldU)
    for tt in 1:length(fldU)
        temp[tt] = calc_UV_conv(fldU[tt], fldV[tt])
    end
    return temp
end
function calc_UV_conv(fldU::MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}, 
                        fldV::MeshArrays.gcmarray{Float32, 2, Matrix{Float32}})
    temp = similar(fldU)
    nlvls = size(fldU, 2)

    for lvl in 1:nlvls
        temp[:, lvl] = MeshArrays.convergence(fldU[:, lvl], fldV[:, lvl])
    end
    return temp
end

function plot_resids_!(ax)
    NaN
end 

end

# @load datadir("ECCO_vars/ADVx_TH_iter129_bulkformula.jld2")