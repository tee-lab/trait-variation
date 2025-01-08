# Using the csv file generated from indi_prop_all_dist_trait.py, we plot the initial and final distribution of traits
library(ggplot2)
library(ggthemes) # For theme_clean() in ggplot

d = read.csv('/home/T_indi_prop_10_types_beta_dist_varying_th_low_var_beta_0.8.csv')
# Replace the '/home/' with the folder address where the files are saved
# In dataframe 'd', each column refers to a particular initial condition (or initial G) and the rows refer to each trait type
# column initial_prop refers to the relative proportion of trait values, and not the actual proportion cover 
# of tree types.
d = d[,-1]

Gseq = seq(0,0.9,0.1) # Same as Gseq in indi_prop_all_dist_trait.py

name = c() 
for (i in Gseq) {
  name=append(name, paste('G0_',i, sep = ''))
}
colnames(d)= c('trait_value','initial_prop',name) 

# Plot for initial G = 0.1 using the column 'G0_0.1' in the dataframe d. Change 'G0_0.1' to other values
# to see the distribution of T for other values of G0.
G0 = 0.1
T0 = (1 - G0)/2
ggplot(data = d, aes(x=trait_value)) + 
  geom_bar(aes(y=G0_0.1), stat="identity", position="identity",fill='mediumpurple3', width=0.01)+
  # plotting initial distribution:
  geom_bar(aes(y=(initial_prop*T0)), stat="identity", position="identity",  alpha=0.5, fill='grey', width=0.03)+
  scale_y_continuous(name="Proportion cover of trees") +
  scale_x_continuous(name="Trait value",breaks=seq(0.05,0.95,0.1), limits=c(0,1)) + 
  # Plots initial G value (useful when looking at variation in theta):
  geom_vline(aes(xintercept = G0), col="darkgreen", linetype="dashed", size=0.7) + 
  # You can also plot the steady-state G value (from the bifurcation diagram) using the above line of code
  # by changing the xintercept value
  theme_clean(base_size = 15)+
  ggtitle("Initial G = 0.1")
