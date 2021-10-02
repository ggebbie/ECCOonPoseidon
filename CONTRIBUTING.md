All users are invited and encouraged to make their own branches to the repository and to submit pull requests to bring their changes into the main branch. Please do submit Issues on GitHub for any bugs or inconveniences, no matter how big or small. 

A recommended workflow is the following.
1. Clone this repository on to a machine with access to poseidon or batou. 
2. Make your own branch (or multiple branches).
3. Push your commits to GitHub often (or whenever it is convenient).
4. Development of new code will likely involve making changes to ECCOtour.jl. ECCOonPoseidon contains code specific to WHOI; ECCOtour.jl is more general code that can be used by any user of ECCO v4r4. 
5. I recommend using `develop --local ECCOtour` so that you can modify ECCOtour without needing to push your commits to GitHub. ECCOtour.jl will be cloned to the `ECCOonPoseidon/dev` directory.
6. When you are satisfied with your edits to ECCOtour.jl, make a branch for that repository as well, and push your commits to GitHub. Then you can point ECCOonPoseidon to use your branch of ECCOtour via something like these commands in the Julia package manager:  `rm ECCOtour`, `add https://github.com/ggebbie/ECCOtour.jl#yourbranch`. 
7. Alternatively, you can clone ECCOtour.jl in your favorite location on your machine and point the Julia package manager to that location via: `dev /location/of/ECCOtour.jl`. Then the ECCOonPoseidon environment will immediately know about your revisions without needing to push to GitHub.
8. Happy hacking! We look forward to collaborating with you.
