#100 types Individual proportion
library(reshape2)
library(ggrepel)
library(rBeta2009) #for rbeta
library(ggplot2)
library(ggthemes)

{dataf = read.csv('T_Indi_prop_th_unif_MODIFIED.csv') #Load the data 
dataf = dataf[,-1]# Remove the first column if it contains the indices

# for (i in 2:ncol(dataf)) {
#   #s = sum(dataf[i,2:(ncol(dataf)-1)])
#   #for (j in 2:(ncol(dataf)-1)) {
#   s = sum(dataf[1:nrow(dataf),i])
#   for (j in 1:nrow(dataf)) {
#     if(s!=0){
#       dataf[j,i] = dataf[j,i]/s
#     }
#   }
# }
# dataf=t(dataf)
# dataf= dataf[-1,]
#Gi=seq(0,0.9,0.1)
#dataf=cbind(dataf,Gi)

distr = c('unif', 'betahi','betalo','bimod') # Set of distributions 
dis = distr[1] # Choose the distribution here
bins = 10 #No. of Types
var_range = c('high','low')
var0 = var_range[1]
if (var0=='low'){
  min_var = 0.3
  max_var = 0.7
  diff = max_var - min_var
  binsize = (diff)/(bins)
  q = seq(min_var+(binsize/2),max_var-(binsize/2),binsize)
}
if (var0=='high'){
  q= seq((1/(2*bins)),1,(1/bins))
}

p=numeric()
s=0
if (dis == "unif"){
  for(i in q){
    s=s+dunif(i,0,1)
  }
  for (i in 1:bins) {
    p[i] = dunif(q[i],0,1)/s
  }
} else if(dis=="betalo") {
  for(i in q){
    s=s + dbeta(i,4,4)
  }
  for (i in 1:bins) {
    p[i] = dbeta(q[i],4,4)/s
  }
} else if(dis=="betahi") {
  for(i in q){
    s=s + dbeta(i,2,2)
  }
  for (i in 1:bins) {
    p[i] = dbeta(q[i],2,2)/s
  }
} else if(dis=="bimod") {
  for(i in q){
    s= s + dbeta(i,2,8)/2+dbeta(i,8,2)/2
  }
  for (i in 1:bins) {
    p[i] = (dbeta(q[i],2,8)/2+dbeta(q[i],8,2)/2)/s
  }
}
p=p*0.25
#p[length(p)+1]=999.999
dataf = cbind(dataf,p)
#q=seq(0,0.9,0.1)
dataf=cbind(dataf,q)}
#q2= round(q,3)
#q2= c(q2,'Gi')
#colnames(dataf) = q2
d2 = melt(dataf,id.vars = 'q')
d2$variable = rep(c(seq(0,0.9,0.1),"At t=0"), each=length(q))
colnames(d2)[colnames(d2)=="q"] = "Tree th value"
colnames(d2)[colnames(d2)=="variable"] = "Initial G"
d3 = dataf[,c(2,5,8,11,12)]
d4=melt(d3,id.vars = 'q')
ggplot(d4, aes(x=q,fill=variable,y=value))+
  geom_bar(color='#e9ecef', alpha=0.6, position='identity')
  
ggplot(dataf,aes(x=q)) +
  geom_bar(aes(y=X0.1),stat = 'identity',fill='#480089',alpha=0.7,col='#480089',width=0.03)+ geom_vline(xintercept = 0.1, col='#2BDA35',size=1,linetype='dashed')+
  #geom_bar(aes(y=X0.6),stat = 'identity',fill='#00b19a',alpha=0.7,col='#00b19a',width=0.03)+ geom_vline(xintercept = 0.6, col='#2BDA35',size=1,linetype='dashed')+
  #geom_bar(aes(y=X0.8),stat = 'identity',fill='#d95272',alpha=0.7,col='#d95272',width=0.03)+ geom_vline(xintercept = 0.8, col='#2BDA35',size=1,linetype='dashed')+
  #geom_area(aes(y=X0.7),fill='#f58231',alpha=0.7,col='#f58231')+
  #geom_bar(aes(y=p),stat = 'identity',fill='#575861',alpha=0.7,col='#575861',width=0.03)+
  theme_clean(base_size = 17)+
  xlab("v")+ ylab("Tree Cover") +
  #geom_vline(xintercept = th_beta[46,2], col='#05954a',size=1)+
  #ggtitle("Woodland Regime,10 types,Unif Dist,Range:0.35-0.65,Relative Proportions") +
  ggtitle("G0<0.5")+
  #ggtitle("G0 = 0.8")+
  #ggtitle("5 Types, Low Var (0.4-0.6), Relative proportions, Unif Dist")+ , High Var
  labs(fill = "Tree th
       value") +
  scale_fill_viridis_d(3)+
  scale_x_continuous(breaks = seq(0.05,1,0.1),limits = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,0.2))+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

colnames(d2)[colnames(d2)=="q"] = "Sapling th value"
colnames(d2)[colnames(d2)=="variable"] = "Initial G"}
ggplot(d2,aes(x=`Initial G`, y=value, fill=factor(`Sapling th value`))) +
  geom_bar(position = 'stack', stat = 'identity')+
  theme_clean(base_size = 15)+
  xlab("Initial G")+ ylab("Sapling Cover") +
  #ggtitle("Woodland Regime,10 types,Unif Dist,Range:0.35-0.65,Relative Proportions") +
  ggtitle("10 Types, Unif Dist, Range: 0-1, Grassland Regime")+
  #ggtitle("5 Types, Low Var (0.4-0.6), Relative proportions, Unif Dist")+ , High Var
  labs(fill = "Sapling th 
       value") +
  scale_fill_viridis_d()
#scale_x_discrete(breaks = c(seq(0,1,0.1)),'At t=0')
ggplot(d2,aes(x=variable, y=value, fill=factor(Types))) +
  geom_bar(position = 'stack', stat = 'identity')+
  theme_clean(base_size = 10)+
  xlab("Initial G")+ ylab("Tree Cover") +
  #ggtitle("Woodland Regime,10 types,Unif Dist,Range:0.35-0.65,Relative Proportions") +
  ggtitle("10 Types, Uniform Distribution, Range: 0 - 1")+
  #ggtitle("5 Types, Low Var (0.4-0.6), Relative proportions, Unif Dist")+ , High Var
  labs(fill = "Tree Types") +
  scale_fill_viridis_d()
  scale_x_discrete(breaks = c(seq(0,1,0.1)),'At t=0')
load('theta_bif.RData')
d=theta_bif_T[,c(1,2,3)]
d[seq(1,15),2]=NA
d[seq(17,nrow(d)),3]=NA
ggplot(d,aes(x=b))+
  geom_point(aes(y=high_1),size=4,col='#3cb44b')+
  geom_point(aes(y=high_2),size=4,col='#3cb44b')+
  geom_point(aes(x=1.4,y=0.82),col='black',size=6)+
  scale_x_continuous(breaks = seq(0,2,0.4),limits = c(0,2))+ 
  scale_y_continuous(breaks = seq(0.1,1,0.2),limits = c(0,1))+
  xlab("b")+
  ylab("G*")+
  theme_classic(base_size=20)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in th")
