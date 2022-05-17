#100 types Individual proportion
library(reshape2)
library(ggrepel)
library(rBeta2009) #for rbeta
library(ggplot2)
library(ggthemes)

dataf = read.csv('filename.csv') #Load the data 
dataf = dataf[,-1]# Remove the first column if it contains the indices

distr = c('unif', 'beta_uni','beta_bimod') # Set of distributions 
dis = distr[1] # Choose the distribution here
bins = 10 #No. of Trait Types
var_range = c('high','low') # Level of trait variation
var0 = var_range[1]
if (var0=='low'){
  min_var = 0.3
  max_var = 0.7
  diff = max_var - min_var
  binsize = (diff)/(bins)
  q = seq(min_var+(binsize/2),max_var-(binsize/2),binsize) # Set of trait values
}
if (var0=='high'){
  q= seq((1/(2*bins)),1,(1/bins))
}

p=numeric()
s=0
# p[i] tells the proportion of each trait value out of 1
if (dis == "unif"){
  for(i in q){
    s=s+dunif(i,0,1)
  }
  for (i in 1:bins) {
    p[i] = dunif(q[i],0,1)/s
  }
} else if(dis=="beta_uni") {
  for(i in q){
    s=s + dbeta(i,4,4)
  }
  for (i in 1:bins) {
    p[i] = dbeta(q[i],4,4)/s
  }
} else if(dis=="beta_bimod") {
  for(i in q){
    s= s + dbeta(i,2,8)/2+dbeta(i,8,2)/2
  }
  for (i in 1:bins) {
    p[i] = (dbeta(q[i],2,8)/2+dbeta(q[i],8,2)/2)/s
  }
}
p=p*0.25 # Multiplying p with Total Tree/Sapling cover (=0.25) to get the proportion cover of each type
dataf = cbind(dataf,p) # Adding a column with proprtion cover values
dataf=cbind(dataf,q) # Adding a column with the trait values

# The following code plots the steady-state distribution for different initial values of G (G0)
# For example, to obtain the plot for G0=0.2, change 'y=X0.2' and 'xintercept = 0.2'
# Uncomment (by removing the '#') and Comment (by adding '#') the line in the code to plot for various values of G0
ggplot(dataf,aes(x=q)) +
  geom_bar(aes(y=X0.1),stat = 'identity',fill='#480089',alpha=0.7,col='#480089',width=0.03)+ geom_vline(xintercept = 0.1, col='#2BDA35',size=1,linetype='dashed')+
  #geom_bar(aes(y=X0.6),stat = 'identity',fill='#00b19a',alpha=0.7,col='#00b19a',width=0.03)+ geom_vline(xintercept = 0.6, col='#2BDA35',size=1,linetype='dashed')+
  #geom_bar(aes(y=X0.8),stat = 'identity',fill='#d95272',alpha=0.7,col='#d95272',width=0.03)+ geom_vline(xintercept = 0.8, col='#2BDA35',size=1,linetype='dashed')+
  #geom_area(aes(y=X0.7),fill='#f58231',alpha=0.7,col='#f58231')+
  #geom_bar(aes(y=p),stat = 'identity',fill='#575861',alpha=0.7,col='#575861',width=0.03)+ # This is the initial distribution of traits
  theme_clean(base_size = 17)+
  xlab("Trait Value")+ ylab("Tree Cover") +
  scale_x_continuous(breaks = seq(0.05,1,0.1),limits = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,0.2))+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
