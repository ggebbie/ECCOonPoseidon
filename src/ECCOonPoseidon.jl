module ECCOonPoseidon
#
# Define functions that are specific to the ECCO runs stored on Poseidon @ WHOI.

using ECCOtour, DrWatson, GoogleDrive, DelimitedFiles

export fluxdir, rectangle, exprootdir, sig1dir,
    diagdir, listexperiments,
    expnames, expsymbols, regpolesdir, rundir,
    Nino34file, historicalNino34, readNino34,
    sigma1grid

""" function sigma1grid()
    Standard (from Susan Wijffels, WHOI) choice of sigma1 surfaces
# Arguments
- `version`: to distinguish from IsopycnalSurfaces.sigma1grid
# Output
- `σ₁grid`: list (vector) of σ₁ values
"""
function sigma1grid()
   return [24.0000
   24.1500
   24.3000
   24.4500
   24.6000
   24.7500
   24.9000
   25.0500
   25.2000
   25.3500
   25.5000
   25.6500
   25.8000
   25.9500
   26.1000
   26.2500
   26.4000
   26.5500
   26.7000
   26.8500
   27.0000
   27.1500
   27.3000
   27.4500
   27.6000
   27.7500
   27.9000
   28.0500
   28.2000
   28.3500
   28.5000
   28.6500
   28.8000
   28.9500
   29.1000
   29.2500
   29.4000
   29.5500
   29.7000
   29.8500
   30.0000
   30.1500
   30.3000
   30.4500
   30.6000
   30.7500
   30.9000
   31.0500
   31.1750
   31.2500
   31.3250
   31.4000
   31.4750
   31.5500
   31.6400
   31.6800
   31.7200
   31.7600
   31.8000
   31.8400
   31.8800
   31.9200
   31.9600
   32.0000
   32.0400
   32.0800
   32.1200
   32.1600
   32.2000
   32.2400
   32.2800
   32.3200
   32.3600
   32.4000
   32.4080
   32.4160
   32.4240
   32.4320
   32.4400
   32.4480
   32.4560
   32.4640
   32.4720
   32.4800
   32.4880
   32.4960
   32.5040
   32.5120
   32.5200
   32.5280
   32.5360
   32.5440
   32.5520
   32.5600
   32.5680
   32.5760
   32.5840
   32.5920
   32.6000
   32.6080
   32.6160
   32.6240
   32.6320
   32.6400
   32.6480
   32.6560
   32.6640
   32.6720
   32.6800
   32.6880
   32.6960
   32.7040
   32.7120
   32.7200
   32.7280
   32.7360
   32.7440
   32.7520
   32.7600
   32.7680]

    # σ₁grida = 24:0.05:31
    # σ₁gridb = 31.02:0.02:33
    # σ₁grid = vcat(σ₁grida,σ₁gridb)
    # σ₁grid = σ₁grid[1:3:end]
   # return σ₁grid
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

end
