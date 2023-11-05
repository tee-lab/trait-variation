# Using the csv file generated from all_dist_trait.py, we plot bifurcation diagrams

library(reshape2) # For melt()
library(ggplot2)
library(ggthemes) # For theme_clean() in ggplot

Bseq = seq(0,2,0.1) # Same as the Bseq in the python script. List of beta values
data_bif = data.frame('b'=Bseq)

dataf = data.frame('b'=Bseq)

# First, we'll read the csv file with high var
d =  read.csv('/home/G_10_types_unif_dist_varying_u_high_var.csv') # If you have set the working directory as the folder containing these .csv files,
# remove '/home/', otherwise replace the '/home/' with the folder address where the files are saved
d = d[,-1]
d = t(d)
d = d[-1,]

name = c() 
for (i in 1:ncol(d)) {
  name=append(name, paste('high_',i, sep = ''))
}
colnames(d)= name
dataf = cbind(dataf,d)
d_main = melt(dataf, id.vars='b')
d_main$variable = rep('high')

# Second, we'll read the csv file with low var
d =  read.csv('/home/G_10_types_unif_dist_varying_u_low_var.csv') # Replace the '/home/' with the folder address where the files are saved
d = d[,-1]
d = t(d)
d = d[-1,]

name = c() 
for (i in 1:ncol(d)) {
  name=append(name, paste('low_',i, sep = ''))
}
colnames(d)= name
dataf = data.frame('b'=Bseq)
dataf = cbind(dataf,d)
d_temp = melt(dataf, id.vars='b')
d_temp$variable = rep('low')

d_main = rbind(d_main, d_temp)

# For the no variation case, we had used Wolfram Mathematica, 
# but if you want to simulate that numerically, you can use GST_no_var.py
# We'll read the csv file with no variation
d =  read.csv('/home/G_u_0.1_v_0.1_th_0.5_no_var.csv') # Replace the '/home/' with the folder address where the files are saved
d = d[,-1]
d = t(d)
d = d[-1,]

name = c() 
for (i in 1:ncol(d)) {
  name=append(name, paste('no_var_',i, sep = ''))
}
colnames(d)= name
dataf = data.frame('b'=Bseq)
dataf = cbind(dataf,d)
d_temp = melt(dataf, id.vars='b')
d_temp$variable = rep('no_var')

d_main = rbind(d_main, d_temp)

# Plotting the data
ggplot(d_main, aes(x=b, y=value, col=variable))+
  geom_point()+
  xlab("Birth Rate of Saplings")+
  ylab("G*")+
  theme_clean(base_size=12)+
  scale_x_continuous(breaks = seq(0,2,0.2))+
  scale_y_continuous(breaks = seq(0,1,0.1))
