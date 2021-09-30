# Check that the output has been copied to batou from poseidon
!isdir(exprootdir(expt)) ? error("Experiment missing") : nothing

diagpath = diagdir(expt)

# print output here
path_out = sig1dir(expt)
!isdir(path_out) ? mkdir(path_out) : nothing;

pathout = regpolesdir(expt)
!isdir(pathout) ? mkdir(pathout) : nothing;

################################################################
# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
γ = setupLLCgrid(datadir("grid/"))
nf = length(γ.fSize)

# get standard levels of MITgcm
z = depthlevels(γ)
pstdz = pressurelevels(z)
p₀ = 1000.0 ; # dbar
