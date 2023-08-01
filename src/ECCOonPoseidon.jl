module ECCOonPoseidon
#
# Define functions that are specific to the ECCO runs stored on Poseidon @ WHOI.

using ECCOtour, DrWatson, GoogleDrive, DelimitedFiles, 
PyCall, PyPlot, MeshArrays, MITgcmTools

export fluxdir, rectangle, exprootdir, sig1dir,
    diagdir, listexperiments,
    expnames, expsymbols, regpolesdir, rundir,
    Nino34file, historicalNino34, readNino34,
    sigma1grid

#statements for budgets.jl
export extract_ocnTAU, extract_eulerian_velocities, 
extract_eulerian_and_bolus_velocities,
extract_lateral_heatbudget, extract_vertical_heatbudget, 
calc_bolus, extract_sθ

#statements for grid_tools.jl
export get_msk, findlatlon, findmin, 
densityJMD95, get_cell_volumes, get_cell_thickness, 
get_geothermalheating,
calc_W_conv3D!, calc_UV_conv3D!, exch_UV_llc90, 
interpolate_to_lateral_faces, interpolate_to_vertical_faces!

#statements for transports.jl
export ThroughFlow, extract_meridional_Ψ_Mean, 
extract_meridional_Ψ, 
extract_meridionalΨ̄timeseries, LatitudeCirclesMask

#statements for mask_tools.jl
export wet_pts, region_mask

#statements for Operations.jl 
export lateral_sum, vertical_sum, zonal_sum, zonal_average

#statements for PacificOcean.jl
export PAC_mask
export get_min_lat, get_max_lat, within_lon, 
get_cs_and_sn, rotate_UV_native, get_ϕ_max_min_mask

import Base: sum
import ECCOtour.sigma1grid

include("Budgets.jl")
include("Grid_Tools.jl")
include("Transports.jl")
include("Mask_Tools.jl")
include("Operations.jl")
include("PacificOcean.jl")

# add a method to this function

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
    push!(shortnames,"seasonalclimatology" => "sznF")

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

end
