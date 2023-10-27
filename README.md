# trait-variation

# trait-variation
We investigate how trait variations affect the dyanmics of savanna-woodland bistable system. 

1. To run the model with variation in the following traits: sapling death rate (u, as referred to in the codes), tree death rate (v) and sapling resistance to fire (th), use the file all_dist_trait.py. Users can also set the initial trait distribution and level of variation using this file. This file gives the steady state grass and tree cover for a specific case of trait variation, trait distribution and level of variation for values of sapling birth rate (b) ranging from 0 to 2. Data is generated in the form of 2 .csv files, one for Grass cover and the other for Tree cover. 
2. Plot the data in these .csv files to obtain the stability diagram and study the ecosystem-level properties. For this, use the plot_bifurcation_diag.R file. This will generate Figure 2 mentioned in the manuscript.

3. To look at the population-level trait distribution at steady state, use indi_prop_all_dist_trait.py. Similar to case 1, users can set the varying trait, initial trait distribution and the level of variation in this file. At a particular value of sapling birth rate (b), this file generates two .csv files with steady state Tree and Sapling cover, respectively corresponding to each trait value present in the population.
4. Plot the data in these .csv files to obtain the population-level trait distribution at steady state using plot_indi_prop.R. This will generate Figure 3 mentioned in the manuscript.

5. To see how grass, tree and sapling covers change with time for different cases of initial conditions (G,S and T), sapling birth rate (b), varying trait, trait distribution and level of variation, use the TimeSeries_all_dist_trait.py. 
6. Plot the data in these .csv files using plot_timeseries.R. This generates Figure 4 mentioned in the manuscript.

7. To run the model for 100 or 1000 sapling or tree types, use the file [INSERT FILE NAME]. This generates Figure S1 in Supplementary Information.
