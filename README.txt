README for R Code STIR Simulations

There are 4 folders:
- functions 
- simulations
- results
- Rmkd

The 'functions' folder includes R source code for used own functions, e.g. STIR(), LSIR().

In 'simulations' there are the R scripts performing the simulations and saving it
(into 'results' folder) in a .rds object containing all relevant results.
Each of these files has the following structure:
- read in needed functions
- define function for data generation
- define function for performing simulations
- define list of parameter settings
- for each parameter setting perform simulation and save results
These files should be run from this folder (set the PATH).

The folder 'results' contains the saved simulation results in .rds format. The folder includes one example, all other ones can be generated from the simulation files or downloaded in the .zip-file in releases.

In 'Rmkd' there are .rmd files reading in the results (.rds files) and creating the Tables and Figures.
