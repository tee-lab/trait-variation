# Trait Variation in a savanna-woodland bistable system

This repository contains codes used to simulate our extension of the savanna-woodland model given by Staver & Levin, 2011, to include trait variations. We investigate how trait variations affect the dyanmics of this bistable system. The codes are part of the manuscript: ___________

The repository contains two folders:
1. [simulations](https://github.com/tee-lab/trait-variation/tree/8aebff27515786bcf7eb8262f35760350e353f1b/simulations) which contains **python code** to numerically simulate the model.
2. [figures](https://github.com/tee-lab/trait-variation/tree/8aebff27515786bcf7eb8262f35760350e353f1b/figures) which contains **R code** to plot the data obtained from running simulations using the codes in *simulations* folder.

Packages required for each of the codes are mentioned in the files. 

## Ecosystem-level Dynamics

### No variation model

To simulate the savanna-woodland model with no variation, use the file [GST_no_var.py](https://github.com/tee-lab/trait-variation/blob/33576342042848f3cf356af291941adec28216ca/simulations/GST_no_var.py) in the *simulations* folder. This file can be run through the Terminal or using an IDE. To run the file through the terminal, follow these steps:
1. Download the file to your system. 
2. Go to the folder containing this file and open terminal in that folder OR Open terminal, type `cd /location-of-the-file/`.
3. Type `python GST_no_var.py` in the terminal and press Enter.

This would generate two .csv files with the steady-state grass cover and tree cover, respectively.
Check the section [Figure 2: Bifurcation Diagram](#figure-2\:-bifurcation-diagram) for instructions on how to plot the bifurcation diagram for this case.

### Trait variation model

To understand the affect of trait variations on the ecosystem-level dynamics of the savanna-woodland model, we use the file [all_dist_trait.py](https://github.com/tee-lab/trait-variation/blob/8aebff27515786bcf7eb8262f35760350e353f1b/simulations/all_dist_trait.py) in the [simulations](https://github.com/tee-lab/trait-variation/tree/8aebff27515786bcf7eb8262f35760350e353f1b/simulations) folder. In this file, one can select 
1. the trait to be varied: sapling death rate ("u", as referred to in the code), tree death rate ("v") or sapling resistance to fire ("th"),
2. initial distribution of traits: uniform ("unif", as referred to in the code), unimodal beta ("beta") or bimodal beta distribution ("bimod"), and
3. the level of variation in the trait: "high" or "low".

After downloading the all_dist_trait.py file, it can either be run directly through the terminal or using an IDE. To run the file directly through the terminal, go to the folder containing the python script, open the terminal and type:
`python all_dist_trait.py -t (insert trait here) -d (insert distribution here) -v (insert level of variation here)`
For example, you want to run the code for variation in sapling death rate (u), which has a bimodal distribution and high level of variation, the command will be:
`python all_dist_trait.py -t u -d bimod -v high`

To run the file in an IDE (such as Spyder), comment out Code Segment 1, as mentioned in the code, and uncomment Code Segment 2 to set the parameters of your choice.

This python script runs the model with variation in the specified trait for different values of sapling birth rate (denoted by "b") and different initial values of Grass cover.
This generates two .csv files with the steady-state grass cover and tree cover, respectively.

### Figure 2: Bifurcation Diagram
To generate Figure 2 mentioned in the manuscript, we plot the data in these .csv files: specifically the files with high, low and no variation in the trait of interest while keeping other traits constant. We use the R script file [Figure2.R](https://github.com/tee-lab/trait-variation/blob/504ac4bc1de92da170a0100e91266c7a6ceed034/figures/Figure2.R) present in the [figures](https://github.com/tee-lab/trait-variation/tree/8aebff27515786bcf7eb8262f35760350e353f1b/figures) folder. This yields a stability/bifurcation diagram, with either the steady-state Grass cover or Tree cover on the y-axis and sapling birth rate on the x-axis. 
To obtain the figure, open Rstudio. In RStudio, set the working directory as the folder containing all these .csv files and run the code.

## Population-level Dynamics

To look at the population-level trait distribution at steady state, use [indi_prop_all_dist_trait.py](https://github.com/tee-lab/trait-variation/blob/a8c5d451990c6ef6574ca172dc08f14fed9f0b8a/simulations/indi_prop_all_dist_trait.py) in the *simulations folder*. Similar to [all_dist_trait.py](https://github.com/tee-lab/trait-variation/blob/8aebff27515786bcf7eb8262f35760350e353f1b/simulations/all_dist_trait.py), users can set the varying trait, initial trait distribution and the level of variation in this file. For a specific value of sapling birth rate (b), this file generates two .csv files with steady state Tree and Sapling cover respectively corresponding to each trait value present in the population. For instance, the varying trait is *sapling death rate* and the level of variation is *high*, the python script will give the tree (and sapling) cover of the 10 tree (and sapling) types with trait value: 0.05, 0.15,..., 0.95. 

### Figure 3: Steady-state Trait Distribution

Plot the data in these .csv files to obtain the population-level trait distribution at steady state using [Figure3.R](https://github.com/tee-lab/trait-variation/blob/a8c5d451990c6ef6574ca172dc08f14fed9f0b8a/figures/Figure3.R). This generate Figure 3 mentioned in the manuscript.

## Timeseries

To see how grass, tree and sapling covers change with time for different initial conditions (i.e., different initial values of G,S and T), sapling birth rate (b), varying trait, trait distribution and level of variation, use the python script [TimeSeries_all_dist_trait.py](https://github.com/tee-lab/trait-variation/blob/a8c5d451990c6ef6574ca172dc08f14fed9f0b8a/simulations/TimeSeries_all_dist_trait.py) in the *simulations* folder.

### Figure 4: Timeseries

Plot the data in the .csv file obtained from TimeSeries_all_dist_trait.py using the R script [Figure4.R](https://github.com/tee-lab/trait-variation/blob/a8c5d451990c6ef6574ca172dc08f14fed9f0b8a/figures/Figure4.R) in the *figures* folder. This generates Figure 4 mentioned in the manuscript.

7. To run the model for 100 or 1000 sapling or tree types, use the file [INSERT FILE NAME]. This generates Figure S1 in Supplementary Information.
