# Using the csv file generated from TimeSeries_all_dist_trait.py, we see how proportion covers of different variables
# change with time
library(ggplot2)
library(ggthemes) # For theme_clean() in ggplot

d = read.csv('/home/TimeSeries_10_types_beta_dist_varying_th_low_var_beta_0.8.csv')
# Replace the '/home/' with the folder address where the files are saved
d = d[,-1]
d_ts = melt(d, id.vars='Time')

ggplot(d_ts, aes(x=Time, y=value, col=variable))+
  geom_line() +
  xlab("Time")+
  ylab("Proportion Cover")+
  theme_clean(base_size=12)+
  #scale_x_continuous(breaks = seq(0,2,0.2))+
  scale_y_continuous(breaks = seq(0,1,0.1))
