println("Configure the experiment variables and directories")

# Check that the output has been copied to batou from poseidon
!isdir(exprootdir(expt)) && error("Experiment missing")

if isdir(rundir(expt)) #only look at these paths if the experiment has already been run
    diagpath = diagdir(expt)
    # print output here
    path_out = sig1dir(expt)
    !isdir(path_out) ? mkdir(path_out) : nothing;

    pathout = regpolesdir(expt)
    !isdir(pathout) ? mkdir(pathout) : nothing;
end

################################################################
# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
pth = MeshArrays.GRID_LLC90
γ = GridSpec("LatLonCap",pth)
Γ = GridLoad(γ;option="full")

# no longer needed?
#γ = setupLLCgrid(datadir("grid/"))

nf = length(γ.fSize)

# get standard levels of MITgcm
p₀ = 1000.0; # dbar
z = depthlevels(γ)
pstdz = pressurelevels(z);
