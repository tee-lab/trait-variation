library(reshape2)
library(ggrepel)
library(ggplot2)
library(ggthemes)

dataf = read.csv("filename.csv")
dataf=dataf[,-1] # Remove the column with row indices

# For plotting Bif Diagrams
library(reshape2)
library(ggrepel)
library(rBeta2009) #for rbeta
library(ggplot2)
library(ggthemes)
library(plotly)
load("types_u_T.RData")
Gi = seq(0,0.9,0.1)
dataf  = data.frame('Gi'=Gi)
Bseq = seq(0,2,0.01)
dataf = data.frame()
x=seq(1,8)
x=c('1','2','_new','_old')
x=c(0,0.3,0.6,1)
#no_b = seq(0.025,1,0.05)
for (i in 1:(length(x))) {
  #if (i %in% no_b){
    #next
 # }
  #else {
  
    d = read.csv(paste('Si_10ty_Si_unif_Ti_unif_vary_uandth_b_0_2_full',x[i],'.csv', sep = ''))
    #d = read.csv(paste('T_1000ty_Si_unif_Ti_unif_fullrange_vary_u_btill',x[i],'.csv', sep = ''))
    #dataf=cbind(dataf,d$Gf)
    dataf=cbind(dataf,d[,seq(3,ncol(d))])
    #colnames(dataf)[ncol(dataf)] = paste(i)
  #}
}
dataf1 = dataf
dataf=cbind(dataf,dataf1[,seq(2,ncol(dataf1))])

{dataf = read.csv("TimeSeries_full_woodland_th_unif_MODIFIED.csv")
dataf=dataf[,-1]}
dataf = cbind(dataf, theta_bif[,c(6,7,8)])
u_beta_T = dataf
save(u_beta_T, file = "u_beta_T.RData")}

t_dataf= t_dataf
t_dataf = dataf[,c(1,2)]
t_dataf["S"] = rowSums(dataf[,seq(3,12)])
t_dataf["T"] = rowSums(dataf[,seq(13,22)])
dataf = dataf[1:70000,]
load('types_th.RData')
uandth_narrow=cbind(uandth_narrow,dataf)
uandth_full=cbind(uandth_full,dataf)

dat = dataf
for( i in 1:nrow(dat)){
  for (j in 1:ncol(dat)) {
    dat[i,j]= as.numeric(round(as.numeric(dat[i,j]),4))
  }
}
dataf = dat
d2 = dat[c(1,10),-1]}
#d2 = unique(d2)}
#dat = dat[,-1]
# colnames(d_unif_100) = c("Gi",Bseq)
d2 = dat[1:200000,]
dataf['Types']=types
{d2 = melt(t_dataf, id.vars='Time')
#d2$variable = rep(seq(0,0.9,0.1), each=10)
data_label <- d2                            # Modify data
data_label$label <- NA
data_label$label[which(data_label$Time == max(data_label$Time))] <- data_label$variable[which(data_label$Time == max(data_label$Time))]
#data_label$label[which(data_label$Time == mean(data_label$Time))] <- data_label$variable[which(data_label$Time == mean(data_label$Time))]
#niter=50000
niter=nrow(dataf)
l=seq(niter,21*niter,niter)
#l_name = c("G","S1","S2",'S3','S4','S5','S6','S7','S8','S9','S10',"T1","T2",'T3','T4','T5','T6','T7','T8','T9','T10')
l_name = c("G","0.05","0.15",'0.25','0.35','0.45','0.55','0.65','0.75','0.85','0.95',"0.05","0.15",'0.25','0.35','0.45','0.55','0.65','0.75','0.85','0.95')
#l_name = c("G","0","0.1",'0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9',"0","0.1",'0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9')
l_name = c("G","S","T")
l=seq(niter,3*niter,niter)
#l_name = c("G","S","T")
# 
data_label$label= replace(data_label$label,l,l_name)}
palette_blues <- colorRampPalette(colors = c("#5cc2af", "#0a6265"))(10)
scales::show_col(palette_blues)
#palette_blues <- colorRampPalette(colors = c("#70b837", "#004b88"))(10)
#scales::show_col(palette_blues)
palette_reds <- colorRampPalette(colors = c("#eca14d", "#da4224"))(10)
scales::show_col(palette_reds)
#
colour_pal = c( "#440154FF", viridis::viridis(10), viridis::magma(10))
ggplot(data_label, aes(x=Time,y=value, col=variable)) +    # Draw ggplot2 plot with labels
  xlim(0,niter+10000) + 
  #ylim(0,0.1)+
  scale_color_manual(values = c("#575A59","#04689C", "#CA6F05"))+
  #scale_color_manual(values = c("#575A59",palette_blues,palette_reds))+
  #scale_color_viridis_d(21) +
  geom_line(size=1) + 
  geom_label_repel(aes(label = label, size=27, fontface="bold"), nudge_x = 10000, nudge_y = 0.03,na.rm = TRUE, max.overlaps=Inf) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))+
  #scale_x_continuous(breaks = seq(0,60000,2000), limits = c(0,niter+10000))+
  #scale_x_continuous(breaks = seq(0,50000,10000))+
  #scale_y_continuous(breaks = seq(0,1,0.1))+
  labs(x="Time",y="Proportion cover",
       title="Variation in o, Woodland Regime")+
  #geom_vline(xintercept = 10500, colour='black',linetype="dashed")+
  theme_clean(base_size = 20)+theme(legend.position = "none")

dataf = read.csv("Si_10ty_Si_unif_Ti_unif_vary_uandth_b_0_2_NAR.csv")
dataf = dataf[,-1]
dat = dataf
for( i in 1:nrow(dat)){
  for (j in 1:ncol(dat)) {
    dat[i,j]= as.numeric(round(as.numeric(dat[i,j]),3))
  }
}
dataf = dat
d2 = melt(uandth_full, id.vars = 'b')
d2$variable = rep(c('1','2'), each=201,times=499)
ggplot(T_u_th, aes(x=u, y=th,fill=RelCover))+
  geom_point( shape=1,alpha=0.9)+
  #stat_density_2d(aes(fill = RelCover), geom = "polygon")
  #geom_bar(stat = 'identity')
  #scale_x_continuous(breaks = seq(0,2,0.5))+
  #scale_y_continuous(breaks = seq(0,1,0.050))+
  xlab("Death Rate of trees")+
  ylab("G*")

ggplot(d2,aes(x=b,y=value))+
  geom_point(alpha=0.1,col='darkblue')+
  #ylim(c(0,0.05))+xlim(c(0,2000))+
  #scale_color_manual(values = c('#453064','#453064','#5f64c0','#5f64c0','#3995fc','#3995fc','lightgreen','lightgreen'))+
  #scale_size_manual(values = c(5,5,4,4,3,3,1,1))+
  #scale_x_continuous(breaks = seq(0,1,0.1))+
  #scale_y_continuous(breaks = seq(0,1,0.1))+
  #scale_color_viridis_b(21)+
  xlab("B")+ scale_color_discrete(name="")+
  ylab("G*")+
  theme_clean(base_size=15)+
  ggtitle("10 Types, Varying u and th, Range: 0-1")
  #ggtitle("10 types, Bimodal Distribution, Range:0.3-0.7, 0-1")
  #ggtitle("10 types, Beta Distribution (High Variation), Range: 0-1")

types_th=theta_bif[,c(1,2,3)]
colnames(types_th)=c('b','10ty1','10ty2')
dataf = read.csv("100ty_Si_unif_Ti_unif_vary_th_b_0_2.csv")
dataf = dataf[,-1]
types_th[,'100ty1']=dataf[,2]
types_th[,'100ty2']=dataf[,3]
x=c(0,1,2)
i=1
dataf = read.csv(paste('1000ty_Si_unif_Ti_unif_vary_th_b_',x[i],'_',x[i+1],'.csv', sep = ''))
i=2
d = read.csv(paste('1000ty_Si_unif_Ti_unif_vary_th_b_',x[i],'_',x[i+1],'.csv', sep = ''))
dataf=rbind(dataf,d)
dataf = dataf[,-1]
types_th[,'1000ty1']=dataf[,2]
types_th[,'1000ty2']=dataf[,3]

{i=1
while(i<ncol(d2)) {
  if(d2[1,i]!=d2[2,i]){
    a1=i-1
    break;
  }
  i = i + 1
}
i = ncol(d2)
while(i>0) {
  if(d2[1,i]!=d2[2,i]){
    a2=i+1
    break;
  }
  i = i-1
}
d2[1,1:a1] = NA
d2[2,a2:ncol(d2)] = NA
}

d2 = t(d2)
v0.5_bif_T[,'S1_1000_unif'] = d2[,'10']
v0.5_bif_T[,'S2_1000_unif'] = d2[,'1']

d2 = dat[c(5,7),2:ncol(dat)]
d2[2,1:13] = NA
d2[1,21:ncol(d2)] = NA
d2[1,1:27] = NA
d2[2,31:ncol(d2)] = NA
colnames(d2)=c("S2_u_100","S1_u_100")
bnew = as.numeric(rownames(d2))
d2 = t(d2)
u_bif[,'S1_no_var'] = d2[,'7']
u_bif[,'S2_no_var'] = d2[,'5']
d2 = cbind(d2,bnew)
d2 = d2[,c(1,3,2)]
colnames(d2)=c("b","S1_u_100","S2_u_100")
dmerge2 = melt(d2, id.vars="b")
dmerge2 =dmerge2[-seq(1,60),]
colnames(dmerge2)=c("b","variable","value")
d.merged = rbind(dmerge,dmerge2)
d3 = melt(newfb2, id.vars = 'b')
d.merged = rbind(d.merged,d3[seq(649,810),])
newfb = data.frame(b=seq(0,1,0.05))
newfb4[,'S1_b_100']=d2[,'10']
newfb4[,'S2_b_100']=d2[,'1']
finbif_b[,'S1_b_bi'] = d2[,'10']
finbif_b[,'S2_b_bi'] = d2[,'79']
finbif_b=finbif_b[,c(1,6,7,4,5,2,3,8,9)]
d3 = melt(finbif_b, id.vars=finbif_b$b)
ggplot(d.merged, aes(x=b, y= value,size=variable,colour=variable))+
  geom_line()+
  scale_color_manual(values = c('#453064','#453064','#3995fc','#3995fc','black','black'))+
  scale_size_manual(values = c(5,5,3,3,1,1))+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Birth Rate")+
  ylab("G*")+
  theme_clean(base_size=15)+
  #ggtitle("10 types, Uniform Distribution , Range: 0 - 1")
  #ggtitle("10 types, Bimodal Distribution, Range: 0-1")
  ggtitle("10 types, Beta Distribution (High Variation), Range: 0-1")
d3 = melt(newfb, id.vars = 'b')
ggplot(d3, aes(x=b, y= value, colour=variable))+
  geom_line()+
  scale_color_manual(values = c('#453064','#453064','#5f64c0','#5f64c0','#3995fc','#3995fc','lightgreen','lightgreen'))+
  scale_size_manual(values = c(5,5,4,4,3,3,1,1))+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Birth Rate")+
  ylab("G*")+
  theme_clean(base_size=15)+
  #ggtitle("10 types, Uniform Distribution , Range: 0 - 1")
  #ggtitle("10 types, Bimodal Distribution, Range: 0-1")
  ggtitle("10 types, Beta Distribution (High Variation), Range: 0-1")
#save(finbif_b, file='FinalBifBeta.RData')

d_unif[sapply(d_unif, is.character)]<-lapply(d_unif[sapply(d_unif,is.character)], as.numeric)
d2[sapply(d2, is.integer)]<-lapply(d2[sapply(d2,is.integer)], as.numeric)
typeof(d2[1,3])
class(d2$variable) = 'Numeric'

newfb2 = data.frame('b'=Bseq, 'S1_u'=finbif[,'S1_u'],'S2_u'=finbif[,'S2_u'],'S1_b_bi'=finbif_b[,'S1_b_bi'],
                    'S2_b_bi'=finbif_b[,'S2_b_bi'],'S1_b_h'=finbif_b[,'S1_b_h'],'S2_b_h'=finbif_b[,'S2_b_h'],
                    'S1_nv'=finbif_b[,'S1_nv'],'S2_nv'=finbif_b[,'S2_nv'])
d3 = melt(newfb2, id.vars = 'b')
ggplot(d3, aes(x=b, y= value, size=variable,colour=variable))+
  geom_line()+
  scale_color_manual(values = c('#453064','#453064','#5f64c0','#5f64c0','#3995fc','#3995fc','lightgreen','lightgreen','black','black'))+
  scale_size_manual(values = c(5,5,4,4,3,3,2,2,1,1))+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Birth Rate")+
  ylab("G*")+
  theme_clean(base_size=15)+
  #ggtitle("10 types, Uniform Distribution , Range: 0 - 1")
  #ggtitle("10 types, Bimodal Distribution, Range: 0-1")
  ggtitle("10 types, Beta Distribution (High Variation), Range: 0-1")

Bseq= seq(0,1,0.0125)
newfb3=data.frame('b'=Bseq, 'S1_u'=finbif[,2],
                  'S2_u'=finbif[,3])
dnew = 
dmerge = melt(newfb3,id.vars = 'b')
#'S1_few'=newfb[,'S1_few'],'S2_few'=newfb[,'S2_few'])
newfb3=newfb3[-15,]
newfb3 = cbind(newfb3,'S2_u_100'=d2[,'9'])
d3 = melt(newfb4, id.vars = 'b')
ggplot(d3, aes(x=b, y= value, size=variable,colour=variable))+
  geom_line()+
  scale_color_manual(values = c('#453064','#453064','#5f64c0','#5f64c0','#3995fc','#3995fc','lightgreen','lightgreen'))+
  scale_size_manual(values = c(5,5,4,4,3,3,2,2,1,1))+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Birth Rate")+
  ylab("G*")+
  theme_clean(base_size=15)+
  ggtitle("10 vs 100 types, Beta Distribution , Range: 0 - 1")
  
newfb4=data.frame('b'=Bseq, finbif_b[,c(4,5,6,7,8,9)])
newfb3 = cbind(newfb3,'S2_u_100'=d2[,'9'])
newfb4=newfb4[,c(1,2,3,4,5,8,9,6,7)]
newfb4 = newfb4[,c(1,2,3,6,7,4,5,8,9)]
load("FinalBifUnif.RData")
d1 = melt(finbif, id.vars = 'b')
finbif[,'S1_unif_1_indi'] = d2[,'10']
finbif[,'S2_unif_1_indi'] = d2[,'8']
finbif = finbif[,-1]

ggplot(finbif_b, aes(x=b))+
  geom_line(aes(y=S1_b_1000-0.02), col='#453064', size=0.9)+
  geom_line(aes(y=S2_b_1000-0.02), col='#453064', size=0.9)+
  geom_line(aes(y=S1_b_100), col='#5f64c0', size=0.9)+
  geom_line(aes(y=S2_b_100), col='#5f64c0', size=0.9)+
  geom_line(aes(y=S1_b_h+0.02), col='#3995fc', size=0.9)+
  geom_line(aes(y=S2_b_h+0.02), col='#3995fc', size=0.9)+
  geom_line(aes(y=S1_nv+0.04), col='lightgreen', size=0.9)+
  geom_line(aes(y=S2_nv+0.04), col='lightgreen', size=0.9)+
  # geom_line(aes(y=S1_u_1000-0.02), col='#453064', size=0.9)+
  # geom_line(aes(y=S2_u_1000-0.02), col='#453064', size=0.9)+
  # geom_line(aes(y=S1_u_100), col='#5f64c0', size=0.9)+
  # geom_line(aes(y=S2_u_100), col='#5f64c0', size=0.9)+
  # geom_line(aes(y=S1_u+0.02), col='#3995fc', size=0.9)+
  # geom_line(aes(y=S2_u+0.02), col='#3995fc', size=0.9)+
  # geom_line(aes(y=S1_nv+0.04), col='lightgreen', size=0.9)+
  # geom_line(aes(y=S2_nv+0.04), col='lightgreen', size=0.9)+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Birth Rate")+
  ylab("G*")+
  theme_clean(base_size=16)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Various types, Beta Distribution")
ggtitle("10 types, Analytical Solutions")

save(theta_bif_T,file = 'theta_bif_T.RData')

v0.5_bif = data.frame('b'=Bseq)
##############################################################

ggplot(u_bif, aes(x=b))+
  geom_line(aes(y=Si_beta_Ti_beta_1+0.02), col='#453064', size=0.9)+
  geom_line(aes(y=Si_beta_Ti_beta_2+0.02), col='#453064', size=0.9)+
  geom_line(aes(y=Si_unif_Ti_beta_1-0.02), col='#5f64c0', size=0.9)+
  geom_line(aes(y=Si_unif_Ti_beta_2-0.02), col='#5f64c0', size=0.9)+
  geom_line(aes(y=Si_beta_Ti_unif_1), col='#3995fc', size=0.9)+
  geom_line(aes(y=Si_beta_Ti_unif_2), col='#3995fc', size=0.9)+
  geom_line(aes(y=S1_no_var+0.04), col='lightgreen', size=0.9)+
  geom_line(aes(y=S2_no_var+0.04), col='lightgreen', size=0.9)+
  #geom_line(aes(y=u_bif$`Si_beta_Ti_beta_1,4_1`-0.06), col='orange', size=0.9)+
  #geom_line(aes(y=u_bif$`Si_beta_Ti_beta_1,4_2`-0.06), col='orange', size=0.9)+
  geom_line(aes(y=Si_unif_Ti_unif_1-0.04), col='black', size=0.9)+
  geom_line(aes(y=Si_unif_Ti_unif_2-0.04), col='black', size=0.9)+
  # geom_line(aes(y=S1_u_1000-0.02), col='#453064', size=0.9)+
  # geom_line(aes(y=S2_u_1000-0.02), col='#453064', size=0.9)+
  # geom_line(aes(y=S1_u_100), col='#5f64c0', size=0.9)+
  # geom_line(aes(y=S2_u_100), col='#5f64c0', size=0.9)+
  # geom_line(aes(y=S1_u+0.02), col='#3995fc', size=0.9)+
  # geom_line(aes(y=S2_u+0.02), col='#3995fc', size=0.9)+
  # geom_line(aes(y=S1_nv+0.04), col='lightgreen', size=0.9)+
  # geom_line(aes(y=S2_nv+0.04), col='lightgreen', size=0.9)+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Birth rate of saplings")+
  ylab("G*")+
  theme_clean(base_size=16)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in Death Rate of Saplings")

ggplot(v_bif, aes(x=b))+
  # geom_line(aes(y=S1_u_1000-0.04), col='#453064', size=0.9)+
  # geom_line(aes(y=S2_u_1000-0.04), col='#453064', size=0.9)+
  # geom_line(aes(y=S1_u_100-0.02), col='#5f64c0', size=0.9)+
  # geom_line(aes(y=S2_u_100-0.02), col='#5f64c0', size=0.9)+
  # geom_line(aes(y=S1_unif), col='#3995fc', size=0.9)+
  # geom_line(aes(y=S2_unif), col='#3995fc', size=0.9)+
  geom_line(aes(y=S1_no_var+0.02), col='lightgreen', size=0.9)+
  geom_line(aes(y=S2_no_var+0.02), col='lightgreen', size=0.9)+
  # geom_line(aes(y=S1_b_h), col='orange', size=0.9)+
  # geom_line(aes(y=S2_b_h), col='orange', size=0.9)+
  # geom_line(aes(y=S1_b_l+0.04), col='black', size=0.9)+
  # geom_line(aes(y=S2_b_l+0.04), col='black', size=0.9)+
  # geom_line(aes(y=S1_u_1000-0.02), col='#453064', size=0.9)+
  # geom_line(aes(y=S2_u_1000-0.02), col='#453064', size=0.9)+
  # geom_line(aes(y=S1_u_100), col='#5f64c0', size=0.9)+
  # geom_line(aes(y=S2_u_100), col='#5f64c0', size=0.9)+
  # geom_line(aes(y=S1_u+0.02), col='#3995fc', size=0.9)+
  # geom_line(aes(y=S2_u+0.02), col='#3995fc', size=0.9)+
  # geom_line(aes(y=S1_nv+0.04), col='lightgreen', size=0.9)+
  # geom_line(aes(y=S2_nv+0.04), col='lightgreen', size=0.9)+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Death rate of trees")+
  ylab("G*")+
  theme_clean(base_size=16)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in Birth rate, Uniform Distribution, Various types")
l = c()
ggplot(types_v2[seq(1,54,2),], aes(x=b))+
  geom_point(aes(y=`1000ty2`), col='#262626', size=6,shape=1,stroke=1.5)+
  geom_point(aes(y=`1000ty1`), col='#262626', size=6,shape=1,stroke=1.5)+
  geom_point(aes(y=`100ty1`), col='#F03762', size=3.5,shape=2,stroke=2)+
  geom_point(aes(y=`100ty1`), col='#F03762', size=3.5,shape=2,stroke=2)+
  geom_point(aes(y=`10ty1`), col='#FFC107', size=3,shape=3,stroke=1.5)+
  geom_point(aes(y=`10ty1`), col='#FFC107', size=3,shape=3,stroke=1.5)+
  # geom_line(aes(y=S1_no_var), col='#000004', size=1, alpha=1)+
  # geom_line(aes(y=S2_no_var), col='#000004', size=1, alpha=1)+
  # geom_line(aes(y=US1_no_var), col='#000004', size=1, alpha=1, linetype="dashed")+
  scale_x_continuous(breaks = seq(0,2,0.1),limits = c(0,0.5))+
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
  xlab("Sapling birth rate")+
  ylab("G*")+
  theme_clean(base_size=16)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in v")

ggplot(types_th2, aes(x=b))+
   geom_point(aes(y=`1000ty2`), col='#32008c', size=5.5,shape=16)+
   geom_point(aes(y=`1000ty1`), col='#32008c', size=5.5,shape=16)+
   geom_point(aes(y=`100ty1`), col='#ff6d00', size=3.5,shape=17)+
   geom_point(aes(y=`100ty1`), col='#ff6d00', size=3.5,shape=17)+
   geom_point(aes(y=`10ty1`), col='#bc2f31', size=3,shape='+',stroke=5.5)+
   geom_point(aes(y=`10ty1`), col='#bc2f31', size=3,shape='+',stroke=5.5)+
   geom_line(aes(y=S1_no_var), col='#f0f921', size=0.8, alpha=1)+
   geom_line(aes(y=S2_no_var), col='#f0f921', size=0.8, alpha=1)+
   geom_line(aes(y=US1_no_var), col='#f0f921', size=0.8, alpha=1, linetype="dashed")+
   scale_x_continuous(breaks = seq(0,2,0.1),limits = c(0,0.5))+
   scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
   xlab("Sapling birth rate")+
   ylab("G*")+
   theme_clean(base_size=16)+
   #ggtitle("Various types, Uniform Distribution")
   ggtitle("Variation in th")

ggplot(types_u2_T, aes(x=b))+
  geom_point(aes(y=`1000ty2`), col='#a4a4a4', size=6.5,shape=16)+
  geom_point(aes(y=`1000ty1`), col='#a4a4a4', size=6.5,shape=16)+
  geom_point(aes(y=`100ty1`), col='#6e009a', size=4.5,shape=2)+
  geom_point(aes(y=`100ty1`), col='#6e009a', size=4.5,shape=2)+
  geom_point(aes(y=`10ty1`), col='#0264b2', size=4,shape='+',stroke=5.5)+
  geom_point(aes(y=`10ty1`), col='#0264b2', size=4,shape='+',stroke=5.5)+
  geom_line(aes(y=S1_no_var), col='#f7dc15', size=1, alpha=1)+
  geom_line(aes(y=S2_no_var), col='#f7dc15', size=1, alpha=1)+
  geom_line(aes(y=US1_no_var), col='#f7dc15', size=1, alpha=1, linetype="dashed")+
  scale_x_continuous(breaks = seq(0,2,0.2),limits = c(0,2))+
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
  xlab("Sapling birth rate")+
  ylab("T*")+
  theme_clean(base_size=16)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in u")

d3= melt(types_th2,id.vars = 'b')
ggplot(d3,aes(x=b,y=value,col=variable,shape=variable,size=variable))+
  geom_point()+
  scale_color_manual(values = c('#000000','#49acec','#e15400','#00838e','#00838e','#00838e','#00838e','#00838e','#00838e'))+
  scale_shape_manual(values = c(1,2,3,2,2,2,2,2,2,2))+
  scale_size_manual(values=c(4.5,2.5,1.5,3,3,3,3,3,3,3,3))
  theme_clean()

d1 = data.frame('S1_1000'=c(beta_bif[1:17,'S_1000'], rep(NA, times=(nrow(beta_bif)-17))))
d1 = data.frame('S1_1000'=c(rep(NA, times=17),beta_bif[18:nrow(beta_bif),'S_1000']))
beta_bif[,'S2_1000']= d1[,'S1_1000']
  
ggplot(u0.5_bif,aes(x=b))+
  geom_line(aes(y=S1_10_beta_low), col='#3995fc', size=0.9)+
  geom_line(aes(y=S2_10_beta_low), col='#3995fc', size=0.9)+
  geom_line(aes(y=S1_100_beta_low-0.02), col='#5f64c0', size=0.9)+
  geom_line(aes(y=S2_100_beta_low-0.02), col='#5f64c0', size=0.9)+
  geom_line(aes(y=S1_1000_beta_low-0.04), col='#453064', size=0.9)+
  geom_line(aes(y=S2_1000_beta_low-0.04), col='#453064', size=0.9)+
  geom_line(aes(y=S1_no_var+0.02), col='lightgreen', size=0.9)+
  geom_line(aes(y=S2_no_var+0.02), col='lightgreen', size=0.9)+
  # geom_line(aes(y=S1_10_beta_high-0.02), col='#5f64c0', size=0.9)+
  # geom_line(aes(y=S2_10_beta_high-0.02), col='#5f64c0', size=0.9)+
  # geom_line(aes(y=S1_10_beta_low-0.04), col='#453064', size=0.9)+
  # geom_line(aes(y=S2_10_beta_low-0.04), col='#453064', size=0.9)+
  # geom_line(aes(y=S1_no_var+0.04), col='lightgreen', size=0.9)+
  # geom_line(aes(y=S2_no_var+0.04), col='lightgreen', size=0.9)+
  # geom_line(aes(y=S1_10_bimod+0.02), col='black', size=0.9)+
  # geom_line(aes(y=S2_10_bimod+0.02), col='black', size=0.9)+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Birth Rate of Saplings")+
  ylab("G*")+
  theme_clean(base_size=16)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in u, Beta Dist (Low Var)")

ggplot(ranges_v0.5,aes(x=b))+
  geom_line(aes(y=S_no_var), col='#00b4e1', size=2, alpha=0.7)+
  #geom_line(aes(y=S2_no_var), col='#00b4e1', size=2, alpha=0.7)+
  #geom_point(aes(y=S1_10_unif_narrow), fill='#5f64c0', size=7, shape="_")+
  geom_point(aes(y=S_10_unif_narrow), fill='#5f64c0', size=7, shape="_")+
  geom_point(aes(y=S1_10_unif_modrange), fill='#3995fc', size=5, shape="/")+
  geom_point(aes(y=S2_10_unif_modrange), fill='#3995fc', size=5, shape="/")+
  geom_point(aes(y=S1_10_unif), fill='#453064', size=6.5, shape="*")+
  geom_point(aes(y=S2_10_unif), fill='#453064', size=6.5, shape="*")+
  #geom_point(aes(y=S1_10_bimod), fill='black', size=2, shape=4)+
  #geom_point(aes(y=S2_10_bimod), fill='black', size=2, shape=4)+
  scale_x_continuous(breaks = seq(0,1,0.1),limits = c(0,1))+ 
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0.1,1))+
  xlab("Birth Rate of Saplings")+
  ylab("G*")+
  theme_clean(base_size=12)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in v")

ggplot(dataf,aes(x=b))+
  geom_line(aes(y=S1_no_var), col='#9df100', size=1.5, alpha=1)+
  geom_line(aes(y=S2_no_var), col='#9df100', size=1.5, alpha=1)+
  geom_line(aes(y=US1_no_var), col='#9df100', size=1.5, alpha=1, linetype="dashed")+
  #geom_point(aes(y=S1_10_unif_narrow), fill='#5f64c0', size=7, shape="_")+
  #geom_point(aes(y=S_10_unif_narrow), fill='#5f64c0', size=7, shape="_")+
  geom_point(aes(y=high_1), col='#55007a', size=4, shape=1)+
  geom_point(aes(y=high_2), col='#55007a', size=4, shape=1)+
  geom_point(aes(y=low_1), col='#00838e', size=2.5, shape=2,alpha=0.7)+
  geom_point(aes(y=low_2), col='#00838e', size=2.5, shape=2, alpha=0.7)+
  #geom_point(aes(y=S1_10_bimod), fill='black', size=2, shape=4)+
  #geom_point(aes(y=S2_10_bimod), fill='black', size=2, shape=4)+
  scale_x_continuous(breaks = seq(0,2,0.4),limits = c(0,2))+ 
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
  xlab("b")+
  ylab("T*")+
  theme_clean(base_size=13)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in th")

ggplot(th_beta, aes(x=b))+  
  geom_line(aes(y=S1_no_var), col='#9df100', size=1.5, alpha=1)+
  geom_line(aes(y=S2_no_var), col='#9df100', size=1.5, alpha=1)+
  geom_line(aes(y=US1_no_var), col='#9df100', size=1.5, alpha=1, linetype="dashed")+
  #geom_point(aes(y=S1_10_unif_narrow), fill='#5f64c0', size=7, shape="_")+
  #geom_point(aes(y=S_10_unif_narrow), fill='#5f64c0', size=7, shape="_")+
  geom_point(aes(y=high_1), col='#55007a', size=3, shape=1)+
  geom_point(aes(y=high_2), col='#55007a', size=3, shape=1)+
  geom_point(aes(y=low_1), col='#00838e', size=1.5, shape=2,alpha=0.7)+
  geom_point(aes(y=low_2), col='#00838e', size=1.5, shape=2, alpha=0.7)+
  #geom_point(aes(y=S1_10_bimod), fill='black', size=2, shape=4)+
  #geom_point(aes(y=S2_10_bimod), fill='black', size=2, shape=4)+
  scale_x_continuous(breaks = seq(0,2,0.2),limits = c(0,2))+ 
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
  xlab("b")+
  ylab("G*")+
  theme_clean(base_size=13)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in th")
types_theta=theta_bif[,c(1,8,9,10,11,16,17,18,19)]
d3=melt(types_u0.5_2, id.vars = 'b')
d3= melt(ranges_theta[,c(1,seq(5,9))], id.vars='b')

ggplot(d3, aes(x=b,y=value, group=variable))+
  geom_point(aes(shape=variable, col=variable),size=3)+
  scale_shape_manual(values = c(15,16,16,17,17))+
  #scale_color_brewer(palette="Dark2")+
  #scale_fill_brewer(palette="Dark2")+
#(values = c('#69c856','#69c856','orange','orange','#1d8ebc','#1d8ebc',
 #                              '#430052','#430052','#f0f921','#f0f921'))+
    #values = c('#69c856','#69c856','orange','orange','#1d8ebc','#1d8ebc',
     #                           '#430052','#430052','#f0f921','#f0f921'),alpha=0.4)+
  scale_color_manual(values = c("#0D0887FF", "#7E03A8FF", "#CC4678FF" ,"#F89441FF", "#F0F921FF"))+
  #scale_color_viridis_d(option = "plasma")+
  # scale_shape_discrete(name="Type of Distribution", 
  #                      labels=c("No Variation, u=0.5","No Variation, u=0.5","Uniform","Uniform",
  #   "Beta high var","Beta high var",
  #                               "Beta low var","Beta low var","Bimodal","Bimodal" ))+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Birth Rate of Saplings")+
  ylab("G*")+
  theme_clean(base_size=16)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("10 types, Variation in u, Various Distributions")
  theme(legend.position = "none")
t1=seq(0,1,0.01)
T_types_u0.5_2=T_types_u0.5[1:length(t1),]
T_types_u0.5_2[,1]=t1
for (i in 1:length(t1)) {
  T_types_u0.5_2[i,2:ncol(T_types_u0.5_2)]= T_types_u0.5[which(T_types_u0.5$b==paste(t1[i])),2:ncol(T_types_u0.5)]
}
T_types_u0.5_2[1:20,7]=NA
T_types_u0.5_2[21:nrow(T_types_u0.5_2),6]=NA

save(T_types_u0.5, file = "T_types_u0.5.RData")

ggplot(types_u0.5,aes(x=b))+
  #geom_line(aes(y=S1_no_var), col='#00b4e1', size=1)+
  #geom_line(aes(y=S2_no_var), col='#00b4e1', size=1)+
  # geom_point(aes(y=S1_10_unif), col='#e8a300', fill='#efc35e', size=4.5, shape=2, alpha=1,stroke=1.2)+
  # geom_point(aes(y=S2_10_unif), col='#e8a300', fill='#efc35e', size=4.5, shape=2, alpha=1,stroke=1.2)+
  # geom_point(aes(y=S1_100_unif), col='#2d5cae',fill='#88c4f0', size=4, shape=1, alpha=0.9,stroke=0.8)+
  # geom_point(aes(y=S2_100_unif), col='#2d5cae',fill='#88c4f0', size=4, shape=1, alpha=0.9,stroke=0.8)+
  # geom_point(aes(y=S1_1000_unif), col='#5ac060', size=4, shape=4)+
  # geom_point(aes(y=S2_1000_unif), col='#5ac060', size=4, shape=4)+
  #geom_point(aes(y=S1_10_bimod), fill='black', size=2, shape=4)+
  # #geom_point(aes(y=S2_10_bimod), fill='black', size=2, shape=4)+
  # geom_point(aes(y=S1_10_unif), col='#55007a', fill='#efc35e', size=2, shape=18, alpha=0.5,stroke=1)+
  # geom_point(aes(y=S2_10_unif), col='#55007a', fill='#efc35e', size=2, shape=18, alpha=0.5,stroke=1)+
  # geom_point(aes(y=S1_100_unif), col='#00838e',fill='#88c4f0', size=2, shape=2, alpha=0.4,stroke=1)+
  # geom_point(aes(y=S2_100_unif), col='#00838e',fill='#88c4f0', size=2, shape=2, alpha=0.4,stroke=1)+
  #geom_point(aes(y=S1_1000_unif), col='#f6df23', size=2, shape=20)+
  #geom_point(aes(y=S_1000_unif), col='#9df100',fill='#9df100', size=2, shape=1)+
  geom_point(aes(y=S_1000_unif), size=2, shape=1)+
  geom_point(aes(y=S1_10_unif), size=2, shape=18, alpha=0.5,stroke=1)+
  geom_point(aes(y=S2_10_unif),  size=2, shape=18, alpha=0.5,stroke=1)+
  geom_point(aes(y=S1_100_unif), size=2, shape=2, alpha=0.4,stroke=1)+
  geom_point(aes(y=S2_100_unif), size=2, shape=2, alpha=0.4,stroke=1)+
  scale_color_viridis_b(option = "magma")+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Birth Rate of Saplings")+
  ylab("G*")+
  theme_clean(base_size=12)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in theta, Different no. of types")

ggplot(types_u0.5_2,aes(x=b))+
  #geom_line(aes(y=S1_no_var), col='#00b4e1', size=1)+
  #geom_line(aes(y=S2_no_var), col='#00b4e1', size=1)+
  # geom_point(aes(y=S1_10_unif), col='#e8a300', fill='#efc35e', size=4.5, shape=2, alpha=1,stroke=1.2)+
  # geom_point(aes(y=S2_10_unif), col='#e8a300', fill='#efc35e', size=4.5, shape=2, alpha=1,stroke=1.2)+
  # geom_point(aes(y=S1_100_unif), col='#2d5cae',fill='#88c4f0', size=4, shape=1, alpha=0.9,stroke=0.8)+
  # geom_point(aes(y=S2_100_unif), col='#2d5cae',fill='#88c4f0', size=4, shape=1, alpha=0.9,stroke=0.8)+
  # geom_point(aes(y=S1_1000_unif), col='#5ac060', size=4, shape=4)+
  # geom_point(aes(y=S2_1000_unif), col='#5ac060', size=4, shape=4)+
  #geom_point(aes(y=S1_10_bimod), fill='black', size=2, shape=4)+
  #geom_point(aes(y=S2_10_bimod), fill='black', size=2, shape=4)+
  geom_point(aes(y=S1_10_unif), size=6, shape="/")+
  geom_point(aes(y=S2_10_unif), size=6, shape="/")+
  geom_point(aes(y=S1_100_unif), size=2, shape=1)+
  geom_point(aes(y=S2_100_unif), size=2, shape=1)+
  geom_point(aes(y=S_1000_unif),size=2, shape=2)+
  #geom_point(aes(y=S2_1000_unif), size=6, shape="/")+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  xlab("Birth Rate of Saplings")+
  ylab("G*")+
  theme_clean(base_size=12)+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Variation in theta, Different no. of types")

stable1_u_T = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0., 0.0015015, 0.00297619, 0.00442478, 
                0.00584795, 0.00724638, 0.00862069, 0.00997151, 0.0112994, 0.012605, 
                0.0138889, 0.0151515, 0.0163934, 0.0176152, 0.0188172, 0.02, 
                0.021164, 0.0223097, 0.0234375, 0.0245478, 0.025641, 0.0267176, 
                0.0277778, 0.0288221, 0.0298507, 0.0308642, 0.0318627, 0.0328467, 
                0.0338164, 0.0347722, 0.0357143, 0.036643, 0.0375587, 0.0384615, 
                0.0393519, 0.0402299, 0.0410959, 0.0419501, 0.0427928, 0.0436242, 
                0.0444444, 0.0452539, 0.0460526, 0.046841, 0.047619, 0.0483871, 
                0.0491453, 0.0498938, 0.0506329, 0.0513627, 0.0520834, 0.0527951, 
                0.053498, 0.0541923, 0.0548781, 0.0555557, 0.0562251, 0.0568865, 
                0.0575402, 0.0581861, 0.0588246, 0.0594557, 0.0600797, 0.0606968, 
                0.0613072, 0.0619113, 0.0625093, 0.0631018, 0.0636893, 0.0642726, 
                0.0648526, 0.0654308, 0.066009, 0.0665899, 0.0671774, 0.0677771, 
                0.068398, 0.0690552, 0.0697768, 0.0706281, 0.07186, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA)

stable2_u_T = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                0.420577, 0.433934, 0.44421, 0.453512, 0.462219, 0.470449, 0.478258, 
                0.485683, 0.492753, 0.499494, 0.505929, 0.512077, 0.517958, 0.523589, 
                0.528986, 0.534161, 0.53913, 0.543905, 0.548495, 0.552912, 0.557166, 
                0.561265, 0.565217, 0.569031, 0.572714, 0.576271, 0.57971, 0.583036, 
                0.586255, 0.589372, 0.592391, 0.595318, 0.598155, 0.600909, 0.603581, 
                0.606175, 0.608696, 0.611145, 0.613527, 0.615843, 0.618096, 0.62029, 
                0.622426, 0.624506, 0.626533, 0.628509, 0.630435, 0.632313, 0.634146, 
                0.635935, 0.637681, 0.639386, 0.641052, 0.642679, 0.644269, 0.645823, 
                0.647343, 0.648829, 0.650284, 0.651706, 0.653099, 0.654462, 0.655797, 
                0.657104, 0.658385, 0.65964, 0.66087, 0.662075, 0.663257, 0.664415, 
                0.665552, 0.666667, 0.66776, 0.668834, 0.669887, 0.670921, 0.671937,
                0.672934, 0.673913, 0.674875, 0.67582, 0.676749, 0.677661, 0.678558, 
                0.67944, 0.680307, 0.681159, 0.681998, 0.682823, 0.683634, 0.684432, 
                0.685217, 0.68599, 0.686751, 0.6875, 0.688237, 0.688963, 0.689678, 
                0.690382, 0.691076, 0.691759, 0.692432, 0.693095, 0.693748, 0.694392, 
                0.695027, 0.695652, 0.696269, 0.696877, 0.697476, 0.698068, 0.698651, 
                0.699226, 0.699793, 0.700353, 0.700905, 0.701449, 0.701987, 0.702517, 
                0.703041, 0.703557, 0.704067, 0.704571, 0.705068, 0.705559, 0.706043, 
                0.706522, 0.706994, 0.707461, 0.707922, 0.708378, 0.708827, 0.709272, 
                0.709711, 0.710145, 0.710574, 0.710997, 0.711416, 0.71183, 0.712239, 
                0.712644, 0.713043, 0.713439, 0.71383, 0.714216, 0.714598, 0.714976, 
                0.71535, 0.715719, 0.716085, 0.716446, 0.716804, 0.717158, 0.717508, 
                0.717854, 0.718196, 0.718535, 0.718871, 0.719203, 0.719531, 0.719857, 
                0.720178, 0.720497, 0.720812, 0.721124, 0.721433, 0.721739)

unstable_u_T = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 0.3933, 0.378455, 0.367585, 0.358243, 0.349818, 0.342046, 0.334782, 
                 0.327936, 0.321447, 0.31527, 0.309371, 0.303723, 0.298303, 0.293094, 
                 0.28808, 0.283247, 0.278582, 0.274077, 0.269721, 0.265505, 0.261422, 
                 0.257466, 0.253628, 0.249905, 0.24629, 0.242778, 0.239364, 0.236044, 
                 0.232814, 0.229671, 0.226609, 0.223627, 0.220721, 0.217887, 0.215124, 
                 0.212428, 0.209797, 0.207228, 0.20472, 0.202269, 0.199875, 0.197534, 
                 0.195245, 0.193007, 0.190818, 0.188675, 0.186578, 0.184525, 0.182515, 
                 0.180545, 0.178616, 0.176726, 0.174873, 0.173056, 0.171275, 0.169527, 
                 0.167814, 0.166132, 0.164482, 0.162862, 0.161272, 0.15971, 0.158177, 
                 0.15667, 0.15519, 0.153736, 0.152307, 0.150902, 0.149521, 0.148163, 
                 0.146827, 0.145514, 0.144221, 0.14295, 0.141699, 0.140467, 0.139255, 
                 0.138062, 0.136887, 0.13573, 0.13459, 0.133468, 0.132362, 0.131272, 
                 0.130198, 0.129139, 0.128096, 0.127068, 0.126053, 0.125053, 0.124067, 
                 0.123094, 0.122134, 0.121187, 0.120253, 0.11933, 0.11842, 0.117521, 
                 0.116634, 0.115758, 0.114893, 0.114038, 0.113194, 0.11236, 0.111536, 
                 0.110722, 0.109917, 0.109121, 0.108334, 0.107556, 0.106787, 0.106026, 
                 0.105273, 0.104529, 0.103792, 0.103062, 0.10234, 0.101625, 0.100917, 
                 0.100216, 0.0995207, 0.0988323, 0.0981499, 0.0974735, 0.0968028, 
                 0.0961375, 0.0954775, 0.0948225, 0.0941723, 0.0935266, 0.0928852, 
                 0.0922478, 0.0916142, 0.0909839, 0.0903568, 0.0897324, 0.0891103, 
                 0.0884902, 0.0878716, 0.087254, 0.0866368, 0.0860192, 0.0854007, 
                 0.0847803, 0.0841569, 0.0835293, 0.082896, 0.082255, 0.0816039, 
                 0.0809394, 0.0802572, 0.0795509, 0.0788112, 0.0780225, 0.0771562, 
                 0.0761468, 0.0747431, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA)
u_bif_T[,"S1_no_var"] = stable1_u_T
u_bif_T[,"S2_no_var"] = stable2_u_T
u_bif_T[,"US1_no_var"]=unstable_u_T
ranges_u0.5[,c(3,10)]=u0.5_bif[,c(3,23)]

stable1_theta_T = c(0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0., 0.00653595, 0.0128205, 0.0188679, 0.0246914, 
                    0.030303, 0.0357143, 0.0409357, 0.045977, 0.0508475, 0.0555556, 
                    0.0601093, 0.0645161, 0.0687831, 0.0729167, 0.0769231, 0.0808081, 
                    0.0845771, 0.0882353, 0.0917874, 0.0952381, 0.0985916, 0.101852, 
                    0.105023, 0.108108, 0.111111, 0.114036, 0.116885, 0.119661, 0.12237, 
                    0.125016, 0.127606, 0.130155, 0.132686, 0.13525, 0.137965, 0.141214, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA)
unstable_theta_T = c(NA,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, 0.43366, 0.418536, 0.406006, 0.39487, 
                     0.384684, 0.375223, 0.366351, 0.357978, 0.350041, 0.342489, 0.335284, 
                     0.328395, 0.321795, 0.315461, 0.309374, 0.303517, 0.297875, 0.292433, 
                     0.287181, 0.282106, 0.277198, 0.272448, 0.267848, 0.26339, 0.259066, 
                     0.254869, 0.250794, 0.246833, 0.242982, 0.239235, 0.235588, 0.232035, 
                     0.228572, 0.225194, 0.221899, 0.218681, 0.215537, 0.212464, 0.209457, 
                     0.206514, 0.203631, 0.200805, 0.198032, 0.19531, 0.192633, 0.19, 
                     0.187406, 0.184847, 0.182319, 0.179816, 0.177333, 0.174862, 0.172396, 
                     0.169922, 0.167426, 0.164885, 0.162264, 0.159501, 0.156468, 0.152773, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                     NA, NA, NA, NA)
stable2_theta_T = c(NA,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, 0.49149, 0.506939, 0.520645, 0.533323, 
                    0.545158, 0.556249, 0.566666, 0.576471, 0.585714, 0.594444, 0.602703, 
                    0.610526, 0.617949, 0.625, 0.631707, 0.638095, 0.644186, 0.65, 
                    0.655556, 0.66087, 0.665957, 0.670833, 0.67551, 0.68, 0.684314, 
                    0.688462, 0.692453, 0.696296, 0.7, 0.703571, 0.707018, 0.710345, 
                    0.713559, 0.716667, 0.719672, 0.722581, 0.725397, 0.728125, 0.730769, 
                    0.733333, 0.735821, 0.738235, 0.74058, 0.742857, 0.74507, 0.747222, 
                    0.749315, 0.751351, 0.753333, 0.755263, 0.757143, 0.758974, 0.760759, 
                    0.7625, 0.764198, 0.765854, 0.76747, 0.769048, 0.770588, 0.772093, 
                    0.773563, 0.775, 0.776404, 0.777778, 0.779121, 0.780435, 0.78172, 
                    0.782979, 0.784211, 0.785417, 0.786598, 0.787755, 0.788889, 0.79, 
                    0.791089, 0.792157, 0.793204, 0.794231, 0.795238, 0.796226, 0.797196, 
                    0.798148, 0.799083, 0.8, 0.800901, 0.801786, 0.802655, 0.803509, 
                    0.804348, 0.805172, 0.805983, 0.80678, 0.807563, 0.808333, 0.809091, 
                    0.809836, 0.810569, 0.81129, 0.812, 0.812698, 0.813386, 0.814063, 
                    0.814729, 0.815385, 0.816031, 0.816667, 0.817293, 0.81791, 0.818519, 
                    0.819118, 0.819708, 0.82029, 0.820863, 0.821429, 0.821986, 0.822535, 
                    0.823077, 0.823611, 0.824138, 0.824658, 0.82517, 0.825676, 0.826174, 
                    0.826667, 0.827152, 0.827632, 0.828105, 0.828571, 0.829032, 0.829487, 
                    0.829936, 0.83038, 0.830818, 0.83125, 0.831677, 0.832099, 0.832515, 
                    0.832927, 0.833333, 0.833735, 0.834132, 0.834524, 0.834911, 0.835294, 
                    0.835673, 0.836047, 0.836416, 0.836782, 0.837143, 0.8375, 0.837853, 
                    0.838202, 0.838547, 0.838889, 0.839227, 0.83956, 0.839891, 0.840217, 
                    0.840541, 0.84086, 0.841176, 0.841489, 0.841799, 0.842105, 0.842408, 
                    0.842708, 0.843005, 0.843299, 0.84359, 0.843878, 0.844162, 0.844444, 
                    0.844724, 0.845)
theta_bif_T[,"S1_no_var"] = stable1_theta_T
theta_bif_T[,"S2_no_var"] = stable2_theta_T
theta_bif_T[,"US1_no_var"]=unstable_theta_T
ranges_theta[,c(2,3,10)]=theta_bif[,c(8,9,20)]

stable1_v_T = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0., 0.00180018, 0.00356506, 0.00529568, 
                0.00699301, 0.00865801, 0.0102916, 0.0118946, 0.013468, 0.0150125, 
                0.0165289, 0.018018, 0.0194805, 0.0209171, 0.0223285, 0.0237154, 
                0.0250784, 0.026418, 0.027735, 0.0290298, 0.030303, 0.0315552, 
                0.0327869, 0.0339985, 0.0351906, 0.0363636, 0.037518, 0.0386543, 
                0.0397727, 0.0408739, 0.041958, 0.0430257, 0.0440771, 0.0451128, 
                0.046133, 0.047138, 0.0481283, 0.0491042, 0.0500659, 0.0510137, 
                0.0519481, 0.0528691, 0.0537772, 0.0546726, 0.0555556, 0.0564263, 
                0.0572852, 0.0581324, 0.0589681, 0.0597926, 0.0606061, 0.0614089, 
                0.0622012, 0.0629831, 0.0637549, 0.0645169, 0.0652692, 0.0660121, 
                0.0667459, 0.0674709, 0.0681874, 0.0688959, 0.0695968, 0.0702909, 
                0.0709792, 0.0716627, 0.0723433, 0.0730232, 0.073706, 0.0743967, 
                0.0751031, 0.0758375, 0.0766222, 0.0775014, 0.0785937, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA)
stable2_v_T = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.162521, 
                0.165498, 0.167512, 0.169198, 0.170714, 0.172127, 0.17347, 0.174763, 
                0.176016, 0.177237, 0.178431, 0.1796, 0.180746, 0.181872, 0.182978, 
                0.184065, 0.185134, 0.186186, 0.187221, 0.18824, 0.189243, 0.190231, 
                0.191203, 0.192161, 0.193105, 0.194034, 0.19495, 0.195853, 0.196742, 
                0.197619, 0.198483, 0.199335, 0.200174, 0.201003, 0.201819, 0.202624, 
                0.203419, 0.204202, 0.204975, 0.205737, 0.20649, 0.207232, 0.207965, 
                0.208688, 0.209401, 0.210106, 0.210801, 0.211488, 0.212165, 0.212835, 
                0.213496, 0.214148, 0.214793, 0.21543, 0.216059, 0.21668, 0.217294, 
                0.2179, 0.2185, 0.219092, 0.219677, 0.220256, 0.220827, 0.221392, 
                0.221951, 0.222503, 0.223049, 0.223589, 0.224123, 0.22465, 0.225172, 
                0.225688, 0.226199, 0.226703, 0.227203, 0.227697, 0.228185, 0.228669, 
                0.229147, 0.22962, 0.230088, 0.230552, 0.23101, 0.231464, 0.231913, 
                0.232358, 0.232798, 0.233233, 0.233664, 0.234091, 0.234513)
unstable_v_T = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.155356, 
                 0.151738, 0.149152, 0.146958, 0.14499, 0.143174, 0.141472, 0.13986, 
                 0.13832, 0.136842, 0.135417, 0.134038, 0.132702, 0.131403, 0.130138, 
                 0.128904, 0.1277, 0.126522, 0.12537, 0.124242, 0.123135, 0.12205, 
                 0.120984, 0.119938, 0.118909, 0.117897, 0.116902, 0.115922, 0.114956, 
                 0.114005, 0.113068, 0.112143, 0.111231, 0.110331, 0.109442, 0.108564, 
                 0.107696, 0.106837, 0.105989, 0.105148, 0.104317, 0.103493, 0.102676, 
                 0.101866, 0.101062, 0.100264, 0.099471, 0.0986822, 0.097897, 
                 0.0971147, 0.0963342, 0.0955546, 0.0947748, 0.0939933, 0.0932085, 
                 0.0924184, 0.0916206, 0.0908117, 0.0899875, 0.0891422, 0.0882671, 
                 0.0873491, 0.0863654, 0.0852711, 0.0839472, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
v_bif_T[,"S1_no_var"] = stable1_v_T
v_bif_T[,"S2_no_var"] = stable2_v_T
v_bif_T[,"US1_no_var"]=unstable_v_T
x=seq(1,ncol(dataf)-1)
colnames(dataf) = rep(c('u','th','Cover'), times=(ncol(dataf)-1)/3)
S_u_th_full = data.frame('u'=0,'th'=0,'Cover'=0)
i=1
while(i < length(x) ) {
  S_u_th_full=rbind(S_u_th_full,dataf[,c(x[i],x[i+1],x[i+2])])
  i = i+3
}
S_u_th_full = S_u_th_full[-1,]
RelCover = numeric()
i=1
while(i < nrow(S_u_th_full) ) {
  s = sum(S_u_th_full[seq(i,i+9),3])
  RelCover = append(RelCover, (S_u_th_full[seq(i,i+9),3])/s)
  i=i+10
}
S_u_th_full= cbind(S_u_th_full,RelCover)
ggplot(T_u_th, aes(x=u, y=th,fill=RelCover) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) 
  theme(
    legend.position='none'
  )
ggplot(T_u_th,aes(x=u,y=th,col=RelCover))+
  geom_point(size=10)+
  scale_color_continuous(type = "viridis") +
  theme_bw()
ggplot(S_u_th_full)+
  geom_line(aes(x=u,y=RelCover),size=0.1)+
  scale_color_continuous(type = "viridis") +
  theme_classic(base_size=13)+
  scale_x_continuous(breaks = seq(0,2,0.1),limits = c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
  ylab("Relative Cover of Sapling Type")+
  xlab("u")+
  #ggtitle("Various types, Uniform Distribution")
  ggtitle("Proportion of Tree Types, Var in u and th, 10 types, Range:0-1")

save(S_u_th_full, file="Si_u_th_f.RData")

write.csv(d2,'d2.csv', row.names = FALSE)

load("types_th_T.RData")

types_th_S = types_th
for (i in 1:nrow(types_th)) {
 for (j in 2:ncol(types_th)) {
  types_th_S[i,j] = 1- types_th[i,j] - types_th_T[i,j]   
 }  
}
save(types_th_S,file = "types_th_S.RData")

TbyS_th_ty = types_th
for (i in 1:nrow(types_th)) {
  for (j in 2:ncol(types_th)) {
    if (types_th_S[i,j]==0){
      TbyS_th_ty[i,j] = 100
    } else TbyS_th_ty[i,j] = types_th_T[i,j]/types_th_S[i,j]   
  }  
}

load("u_bif.RData")

v_bif_S = v_bif
for (i in 1:nrow(v_bif)) {
  for (j in 2:ncol(v_bif)) {
    v_bif_S[i,j] = 1- v_bif[i,j] - v_bif_T[i,j]   
  }  
}
save(v_bif_S,file = "v_bif_S.RData")

TbyS_v = v_bif
for (j in 2:ncol(v_bif)) {
  for (i in 1:nrow(v_bif)) {
   if (!is.na(v_bif_S[i,j])) {
     if (v_bif_S[i,j]==0){
        TbyS_v[i,j] = 100
     } else {
        TbyS_v[i,j] = v_bif_T[i,j]/v_bif_S[i,j]
     } 
   }
  }  
}
{load("th_types_T_upd.RData")
dataf = read.csv("1000ty_Si_unif_Ti_unif_vary_v_b_0_2.csv")
dataf=dataf[,-1]
#types_v2_T = v_bif_T[,1:3]
types_v2_T=cbind(types_v2_T, v_bif_T[6:8])
types_v2=cbind(types_v2, v_bif[6:8])
colnames(types_th2_T)=c('b','10ty1','10ty2', '100ty1','100ty2', '1000ty1','1000ty2','S1_no_var','S2_no_var','US1_no_var')}
save(types_v2_T,file = "v_types_T_upd.RData")
