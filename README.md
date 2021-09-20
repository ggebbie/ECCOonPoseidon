# ECCOonPoseidon

Analysis of ECCO version 4 release 4 using output from Poseidon @ WHOI

# Getting started

* from Emacs editor (one possible method)

Install julia-mode, julia-repl, and magit \
Skip the next 5 steps if you have already cloned the repository \
`M-x magit-clone` \
Select `u` to clone from url\
Enter ` https://github.com/ggebbie/ECCOonPoseidon` as url to clone \
Select `y` in response to `remote.pushDefault' to "origin"?` \
Clone to your favorite location and rename project if necessary \
Go to any directory in the project: e.g., `C-x C-f ECCOonPoseidon/scripts`\
Then activate the project and initialize a julia session: `C-c C-a`

* from the Julia REPL

`;`\
`git clone https://github.com/ggebbie/ECCOonPoseidon # only do this the first time on each machine`\
`cd ECCOonPoseidon`\
`]`\
`activate .`\
`instantiate # only do this the first time on each machine`\
To verify you are in the project environment, `]` should return `(ECCOonPoseidon) pkg>`\
Type backspace to return to command mode.

* Using an editor like Atom/Juno or Visual Studio Code, activate the environment on one of the frame panels. The default environment is @v1.x and should be changed.

# Running a script

Add DrWatson to your default enviroment.
One way is to open the REPL and type:\
`]` \
`add DrWatson`\
`cd BASEDIR/ECCOonPoseidon`, where BASEDIR is where you cloned the repository.\
`activate .`\
`instantiate`\
`exit()`

Alternatively, see the instructions under Reproducibility.

From a shell, use the following commands to run a script:\
`cd BASEDIR/ECCOonPoseidon/scripts`, where `BASEDIR` is where you cloned the repository.\
`julia filter_interannual.jl southpac`, one example where the script takes one argument.

# Directory structure
- `scripts`: production-ready scripts
- `src`: definitions of poseidon-specific variables and functions

# Reproducibility

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ECCOonPoseidon

It is authored by G Jake Gebbie.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia REPL and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
