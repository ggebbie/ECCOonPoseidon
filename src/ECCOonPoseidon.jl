module ECCOonPoseidon
#
# Define functions that are specific to the ECCO runs stored on Poseidon @ WHOI.

using ECCOtour

export fluxdir, rectangle, exprootdir, sig1dir, diagdir, listexperiments
export expnames, expsymbols, regpolesdir, rundir

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

end
