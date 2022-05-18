
module OHC_helper
export patchvolume, calc_OHC, get_basin_volumes, standardize, 
       OHC_outputpath, do_FFT, calc_θ_bar, standard_error, conf_int,
       lin_reg, get_basin_depths, get_GH19, plot_patch

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, FFTW, NetCDF

import CairoMakie as Mkie


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
    function get_basin_volumes(cell_area, cell_depths)
    get the volume in area of interst 
# Arguments
- `cell_area`: s
- `cell_depths`: experiment of interest 
# Output
- `basin_volume`: volume-weighted ocean temperature
"""

function get_basin_volumes(cell_area, cell_depths)

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
    function get_basin_depths(basin_mask, Δz, hFacC)
    get true depth for each grid square in area of interest.
# Arguments
- `basin_mask`: mask with entries of interest equal to 1.0 else 0.0
- `Δz`: possible thickness of the calls at each level
- 'hFacC': percentage of possible thickness of each cell 
# Output
- `basin_volume`: volume-weighted ocean temperature
"""

function get_basin_depths(basin_mask, Δz, hFacC)

    cell_depths = similar(hFacC) .* 0.0 
    
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
    function standard_errorlin_reg(x, y_true)
    get standard error for a set of data and predictions
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

end

