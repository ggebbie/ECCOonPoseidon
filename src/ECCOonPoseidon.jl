module ECCOonPoseidon
#
# Define functions that are specific to the ECCO runs stored on Poseidon @ WHOI.

using ECCOtour, DrWatson, GoogleDrive, DelimitedFiles, PyCall, PyPlot, MAT, MeshArrays

export fluxdir, rectangle, exprootdir, sig1dir,
    diagdir, listexperiments,
    expnames, expsymbols, regpolesdir, rundir,
    Nino34file, historicalNino34, readNino34,
    sigma1grid, basin_mask, smooth

# add a method to this function
import ECCOtour.sigma1grid

const mpl = PyNULL()
const plt = PyNULL()
const cmocean = PyNULL()
const cartopy = PyNULL()

#Initialize all Python packages - install with conda through Julia
function __init__()

    # following ClimatePlots.jl
    copy!(mpl, pyimport_conda("matplotlib", "matplotlib", "conda-forge"))
    copy!(cartopy, pyimport_conda("cartopy", "cartopy", "conda-forge"))

    #copy!(plt, pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge"))
    #copy!(cmocean, pyimport_conda("cmocean", "cmocean", "conda-forge"))

    println("Python libraries installed")
 end

""" function sigma1grid()
    Choice of sigma1 surfaces for gridding
# Arguments
- `focus`: what part of water column is the focus?
# Output
- `σ₁grid`: list (vector) of σ₁ values
"""
function sigma1grid(focus::String)

    if focus == "mixed layer"
        σ₁grida = 24:0.15:31.05
        σ₁gridb = 31.175:0.075:31.55
        σ₁gridc = 31.64:0.04:32.4
        σ₁gridd = 32.408:0.008:32.768
        σ₁grid = vcat(σ₁grida,σ₁gridb,σ₁gridc,σ₁gridd)

        return σ₁grid
        
    else
        # call the standard version
        # from ECCOtour/IsopycnalSurfaces
        # has focus = "thermocline"
        
        return sig1grid()
    end
end

fluxdir() = "/batou/eccodrive/files/Version4/Release4/other/flux-forced/forcing/"

fluxdir(expt::String) = "/batou/eccodrive/files/Version4/Release4/other/flux-forced-"*expt*"/forcing/"

"""
    function exprootdir(expt::String) 
    Root directory of the ECCOv4r4 output
"""
function exprootdir(expt::String)
    rootdir = "/batou/ECCOv4r4/exps/"*expt*"/"

    # If the experiment hasn't been copied to batou, look on poseidon.
    !isdir(rootdir) ? rootdir = "/poseidon/ECCOv4r4/exps/"*expt*"/" : nothing
    return rootdir
end

"""
    function exprootdir() 
    Root directory of the ECCOv4r4 output
"""
exprootdir() = "/batou/ECCOv4r4/exps"

rundir(expt::String) = exprootdir(expt)*"run/"
sig1dir(expt::String) = rundir(expt)*"sigma1/"
diagdir(expt::String) = rundir(expt)*"diags/"
regpolesdir(expt::String) = rundir(expt)*"regularpoles/"

function rectangle(region::String)

    if region == "southpac"
        # southpac definitions: move to `src`?
        latrect = (-90, -15) # immutable
        lonrect = [150,-67] # mutable for wraparound 
    elseif region == "test"
        latrect = (0, 10)
        lonrect = [0, 10]
        dlat = 20
        dlon = 20
    else
        error("No matching region, please define region in `src/ECCOonPoseidon.jl`")
    end
    
    if lonrect[1] > lonrect[2] # then handle wraparound
        lonrect[2] += 360  # shift one way, could've shifted the opposite way
    end
    return latrect, lonrect
end

##########################################
# list of experiments on poseidon
# runpath,diagpath = listexperiments(poseidonexpdir());

"""
    function listexperiments(exppath)
    get a dictionary of experiments and their locations
# Arguments
- `exppath`: directory of experiments
# Output
- `runpath`: dictionary with keys=experiments, values=run paths
- `diagpath`: dictionary with keys=experiments, values=diagnostic paths
"""
function listexperiments(exppath)

    # add a trailing / if needed
    exppath[end] != "/" ? exppath *= "/" : nothing
    
    explist = searchdir(exppath,"") # all files in directory
    filter!(x -> !occursin("README",x),explist) # remove README to get explist
    filter!(x -> !occursin(".out",x),explist) # remove output to get explist

    runpath = Dict(explist .=> exppath.*explist.*"/run/")
    diagpath = Dict(explist .=> exppath.*explist.*"/run/diags/")
    regpolespath = Dict(explist .=> exppath.*explist.*"/run/regularpoles/")

    # print to screen all the available
    println("Available experiments are")
    for (key,value) in runpath
        println(key)
    end
    return runpath,diagpath,regpolespath
end

"""
    function expnames()
    Abbreviations/shortnames for some experiments, useful for labels
    Hand-coded and needs manual change with new experiments.
# Output
- `shortnames`: dictionary with keys=experiments, values=abbreviations
"""
function expnames()
    shortnames = Dict("iter129_bulkformula" => "129bf")
    push!(shortnames,"iter0_bulkformula" => "0bf")
    push!(shortnames,"iter129_fluxforced" => "129ff")
    push!(shortnames,"noinitadjust" => "noINIT")
    push!(shortnames,"nosfcadjust" => "noSFC")
    push!(shortnames,"nointerannual" => "noIA")
    return shortnames
end

"""
    function expsymbols()
    List of symbols for some experiments, useful for distinguishing plots
    Hand-coded and needs manual change with new experiments.
# Output
- `marks`: dictionary with keys=experiments, values=abbreviations
"""
function expsymbols()
    marks = Dict("iter129_bulkformula" => "o")
    push!(marks,"iter0_bulkformula" => "x")
    push!(marks,"iter129_fluxforced" => "*")
    push!(marks,"noinitadjust" => ".")
    push!(marks,"nosfcadjust" => "+")
    push!(marks,"nointerannual" => "s")
    return marks
end

"""
    function Nino34file()
    Get location of the historical Nino34 Google Drive file
    Download if necessary
    FUTURE UPDATES: use GoogleDrive package
# Output
- `fileloc`: local location of gdrive file
"""
function Nino34file()
    !isdir(datadir()) && mkdir(datadir())

    # download from google drive and save location
    fileloc = datadir("nino34.hadisst1.1870-2020.txt")
    if !isfile(fileloc) 

        url = "https://docs.google.com/uc?export=download&id=1kOtOnD6B3Y9SAI5W6ezP-o_tgNkonn5n"
        google_download(url,datadir())

        # if GoogleDrive.jl doesn't work, try shell command
        #run(`wget "https://drive.google.com/file/d/1kOtOnD6B3Y9SAI5W6ezP-o_tgNkonn5n/view?usp=sharing" -O $fileloc`)
        
    end
    return fileloc
end

"""
    function historicalNino34(baselineyears)
    Get historical Nino3.4 data and SST climatology from HadISST 1
# Argument
- `baselineyears`: for computation of SST climatology, i.e., (1960,2020)
# Output
- `nino34`: historical Nino3.4 index
- `tnino34`: time in years CE corresponding to index
- `SSTclimatology`: monthly mean values in Nino3.4 patch
"""
function historicalNino34(baselineyears)
    SSTnino34 = readNino34()
    nyr,nmon = size(SSTnino34)

    SSTclimatology = zeros(12) # 12 months
    count = 0
    for i = 1:nyr
        if baselineyears[1] <= SSTnino34[i,1] <= baselineyears[2]
            SSTclimatology += SSTnino34[i,2:end]  # skip the year in column 1
            count += 1
        end
    end
    SSTclimatology = SSTclimatology/count

    # get nino3.4 timeseries
    nino34 = []
    tnino34 = []
    monlist= 1/24:1/12:1
    for i  = 1:nyr
        append!(nino34,SSTnino34[i,2:end]-SSTclimatology)
        append!(tnino34,SSTnino34[i,1] .+ collect(monlist))
    end
    return nino34,tnino34,SSTclimatology
end

"""
    function readNino34()
    Get historical Nino3.4 data from HadISST 1
# Output
- `SST_nino34`: local location of gdrive file
"""
function readNino34()
    filename = Nino34file()
    SST_nino34 = DelimitedFiles.readdlm(filename)
    return SST_nino34
end

"""
    function basin_mask(basin_name,hemisphere)
# Arguments
- `basin_name`: string options are Arctic, Atlantic, Baffin Bay, Barents Sea, Bering Sea,
East China Sea, GIN Seas, Gulf, Gulf of Mexico, Hudson Bay, indian, Japan Sea, Java Sea,
Mediterranean Sea, North Sea, Okhotsk Sea, Pacific, Red Sea, South China Sea, Timor Sea.
-'hemisphere': optional argument. 0 = North, 1 = South, 2 = both
# Output
- 'mask': space and time field of surface forcing, value of zero inside
designated lat/lon rectangle and fading to 1 outside sponge zone on each edge. This is
because this field ends up being SUBTRACTED from the total forcing
""" 

# function basin_mask(basin_name,hemisphere)

#     pth = MeshArrays.GRID_LLC90
#     γ = GridSpec("LatLonCap",pth)
#     Γ = GridLoad(γ;option="full")
#     basins=read(joinpath(pth,"v4_basin.bin"),MeshArray(γ,Float32))

#     basin_list=["Pacific","Atlantic","indian","Arctic","Bering Sea",
#             "South China Sea","Gulf of Mexico","Okhotsk Sea",
#             "Hudson Bay","Mediterranean Sea","Java Sea","North Sea",
#             "Japan Sea", "Timor Sea","East China Sea","Red Sea",
#             "Gulf","Baffin Bay","GIN Seas","Barents Sea"];

#     # Γ.hFacC[:,1] can be used as an indicator for wet points
#     # (there might be a better way to do this)
#     ocean_mask = Γ.hFacC[:,1]
#     ocean_mask[findall(ocean_mask.>0.0)].=1.0
#     if hemisphere == 0 #North
#     hemisphere_mask = Γ.YC .> 0.0;
#     elseif hemisphere == 1 #South
#         hemisphere_mask = Γ.YC < 0.0;
#     elseif hemisphere == 2 #both
#         hemisphere_mask = Γ.YC > 0.0 | Γ.YC <= 0.0; #optional argument?
#     end

#     basinID=findall(basin_list.==basin_name)[1]
#     basin_mask=similar(basins)
#     for ff in 1:5
#         basin_mask[ff] .= hemisphere_mask[ff].*ocean_mask[ff].*(basins[ff].==basinID)
#     end

#     return basin_mask 

# end

function basin_mask(latpt,lonpt,basin_name,hemisphere)
    file = matopen(datadir("basin_grids/GRID_LLC90_"*basin_name))
    for ii in 1:5
        mask[ii] = read(file,basin_name*"_mask"*string(ii))
    end 
    close(file)

    return mask
end

basinlist()=["Pacific","Atlantic","Indian","Arctic","Bering Sea",
                "South China Sea","Gulf of Mexico","Okhotsk Sea",
                "Hudson Bay","Mediterranean Sea","Java Sea","North Sea",
                "Japan Sea", "Timor Sea","East China Sea","Red Sea",
                "Gulf","Baffin Bay","GIN Seas","Barents Sea"]

function basin_mask(basin_name::String,γ)
    pth = MeshArrays.GRID_LLC90
    Γ = GridLoad(γ;option="full")
    basins=read(joinpath(pth,"v4_basin.bin"),MeshArray(γ,Float32))
    basin_list=ECCOonPoseidon.basinlist()
    basinID=findall(basin_list.==basin_name)[1]
    basinmask=similar(basins)
    for ff in 1:5
        basinmask[ff] .= (basins[ff].==basinID) 
    end
    land2nan!(basinmask,γ)
    return basinmask
end

"""
    function smooth(msk::MeshArrays.gcmarray,lengthscale)

    Smooth a gcmarray with a lengthscale of `X` points

    Based off Gael Forget, MeshArrays.jl
"""
function smooth(msk::MeshArrays.gcmarray,X,γ)
    Γ = GridLoad(γ;option="full")
    DXCsm=X*Γ.DXC; DYCsm=X*Γ.DYC;
    #apply smoother
    land2nan!(msk,γ)
    return msk_smooth=MeshArrays.smooth(msk,DXCsm,DYCsm,Γ);
end

function apply_hemisphere_mask!(mask,hemisphere,γ)
    Γ = GridLoad(γ;option="full")
    if hemisphere == :north
        hemisphere_mask = Γ.YC .> 0.0;
    elseif hemisphere == :south #South
        hemisphere_mask = Γ.YC .< 0.0;
    elseif hemisphere == :both #both
        hemisphere_mask = Γ.YC .> 0.0 || Γ.YC .≤ 0.0; #optional argument?
    else
        error("no definition for hemisphere")
    end

    # need loop to assign/mutate mask
    for ff in 1:5
        mask[ff] .*= hemisphere_mask[ff]
    end
    #land2nan!(mask,γ)
end

"""
function basin_mask(basin_names,hemisphere)
# Arguments
- `basin_name`: vector of strings. string options are Arctic, Atlantic, Baffin Bay, Barents Sea, Bering Sea,
East China Sea, GIN Seas, Gulf, Gulf of Mexico, Hudson Bay, indian, Japan Sea, Java Sea,
Mediterranean Sea, North Sea, Okhotsk Sea, Pacific, Red Sea, South China Sea, Timor Sea.
-'hemisphere': optional argument. 0 = North, 1 = South, 2 = both
# Output
- 'mask': space and time field of surface forcing, value of zero inside
designated lat/lon rectangle and fading to 1 outside sponge zone on each edge. This is
because this field ends up being SUBTRACTED from the total forcing
""" 
function basin_mask(basin_names::Vector,γ;hemisphere=nothing,Lsmooth=nothing)
    pth = MeshArrays.GRID_LLC90
    Γ = GridLoad(γ;option="full")
    basins=read(joinpath(pth,"v4_basin.bin"),MeshArray(γ,Float32))
    mask = 0 * basins # needs NaN on land
    for (ii,nn) in enumerate(basin_names)
        mask += basin_mask(nn,γ)
    end

    if !isnothing(hemisphere)
        apply_hemisphere_mask!(mask,hemisphere,γ)
    end

    if !isnothing(Lsmooth)
        mask = smooth(mask,Lsmooth,γ)
    end
    
    return mask
end

end
