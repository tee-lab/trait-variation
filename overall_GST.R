# Program for generating the total G, S and T cover for various cases: 3,5,10 and 100 types of trees/saplings from two levels of variation: High (0.3-0.7) and 
# Low (0.4-0.6)

library(reshape2)
library(ggrepel)
library(rBeta2009) #for rbeta
library(ggplot2)

omega =  function(G,th){
  0.9 + (0.05-0.9)/(1+exp((th-G)/0.005))
}

Gseq = seq(0,1,0.1)
cat = c("Low", "High")
binseq = c(3,5,10,100)
len =length(Gseq)*length(cat)*length(binseq)
df_sc = data.frame(rep(Gseq, times=(length(cat)*length(binseq))),rep(0,times=len),rep(0,times=len),
                   rep(0,times=len),rep(binseq,each=(length(Gseq)*length(cat))),rep(rep(cat,each=(length(Gseq)),times=length(binseq))))
colnames(df_sc)= c('Gi','FinS','FinT','FinG','Types','Var')

for (e in binseq) {
  print("bins=")
  print(e)
  for (ty in cat) {
    for(r in Gseq) {
      print(r)
      bins = e
      if(ty=='High'){
        min_var = 0.25
        max_var = 0.75
        x = max_var - min_var
        binsize = (x-0.1)/(bins-1)
        q= seq(0.3,0.7,binsize) 
        th= seq(min_var, max_var,0.0001) 
      } else {
        min_var = 0.35
        max_var = 0.65
        x = max_var - min_var
        binsize = (x-0.1)/(bins-1)
        q= seq(0.4,0.6,binsize)
        th= seq(min_var, max_var,0.0001) 
      }
      l= dunif(th,min_var,max_var)
      #l= dbeta(th,4,4)
      a =  numeric(bins)
      for (i in 1:bins) {
        a[i]=mean(l[((i-1)*((length(th)-1)/bins) +1):((i)*((length(th)-1)/bins)+1)])
        a[i] = (a[i]*dunif(q[i],min_var,max_var))
        #a[i] = (a[i]*dbeta(q[i],4,4))
      }
      plot(a)
      p = a/sum(a)
      x = 10000
      pr = p*x
      a_beta1 = numeric(0)
      for (i in 1:length(pr)) {
        a_beta1 = append(a_beta1, rep(q[i], times=pr[i] ))
      }
      
      d=length(a_beta1)
      hist(a_beta1)
      #mean(a_beta1)
      
      Gi = G0 = r
      S0 = T0 = (1-G0)/2
      b = 0.35; u <- 0.1; v <- 0.1
      (G0 + S0 + T0)  
      S_t = Si = array(S0/d,dim = d)
      T_t = Ti = array(T0/d,dim = d)
      G = array(G0,dim = d)
      
      tim = 200
      dt = 0.001
      timesteps = tim/dt
      df = data.frame(a_beta1,Si,Ti)
      colnames(df)= c('th','FinS','FinT')
      start_time <- Sys.time()
      ####   Running the loop  ###########
      for(t in 0:(timesteps)){
        G_t<-G0+(u*S0 + v*T0 - b*G0*T0)*dt;
        for (i in q) {
          S_i = df$FinS[which(df$th==i)][1] + (b*G0*df$FinT[which(df$th==i)][1] - omega(G0,i)*df$FinS[which(df$th==i)][1] - u*df$FinS[which(df$th==i)][1])*dt;
          T_i = df$FinT[which(df$th==i)][1]+ (omega(G0,i)*df$FinS[which(df$th==i)][1] - v*df$FinT[which(df$th==i)][1])*dt; 
          
          df$FinS[which(df$th==i)] = S_i
          df$FinT[which(df$th==i)] = T_i
          }
       
        G0<-G_t
        S0 <- sum(df$FinS)
        T0 <- sum(df$FinT)
      }
      end_time <- Sys.time()
      print(end_time -start_time)
      df_sc$FinG[which(df_sc$Gi==r & df_sc$Types==e & df_sc$Var==ty)] = G0
      df_sc$FinS[which(df_sc$Gi==r & df_sc$Types==e & df_sc$Var==ty)] = S0
      df_sc$FinT[which(df_sc$Gi==r & df_sc$Types==e & df_sc$Var==ty)] = T0
    }
  }
}  
# For no variation case #
Gseq = seq(0,1,0.05)
df_sc2 = data.frame(Gseq,rep(0,times=length(Gseq)),rep(0,times=length(Gseq)),
                    rep(0,times=length(Gseq)),rep(1,times=length(Gseq)),rep('No Var',times=length(Gseq)))
colnames(df_sc2)= c('Gi','FinS','FinT','FinG','Types','Var')
### For n types (n = bins) ###
for(r in Gseq) {
  Gi = G0 = r
  S0 = T0 = (1-G0)/2
  b = 0.35; u <- 0.1; v <- 0.1
  (G0 + S0 + T0)  
  
  tim = 100
  dt = 0.001
  timesteps = tim/dt
  #### Running the loop  ###########
  for(t in 0:(timesteps)){
    G_t<-G0+(u*S0 + v*T0 - b*G0*T0)*dt;
    S_t<-S0+(b*G0*T0 - omega(G0,0.5)*S0 - u*S0)*dt;
    T_t<-T0+ (omega(G0,0.5)*S0 - v*T0)*dt;
    G0<-G_t
    S0 <- S_t
    T0 <- T_t
  }
  df_sc2$FinG[which(df_sc2$Gi==r & df_sc2$Types==1 & df_sc2$Var=='No Var')] = G0
  df_sc2$FinS[which(df_sc2$Gi==r & df_sc2$Types==1 & df_sc2$Var=='No Var')] = S0
  df_sc2$FinT[which(df_sc2$Gi==r & df_sc2$Types==1 & df_sc2$Var=='No Var')] = T0
}

df_s <- rbind(df_sc,df_sc2)
dff=df_s[which(df_s$Types==1),]
dff = df_s[which(df_s$Var=='Low'),]

ggplot(dff,aes(x=Gi,y=FinT))+
  #geom_point(size=2.5)+
  geom_point(size=3,aes(color = factor(Var)))+
  ggtitle('For 1 type')+
  #ggtitle('For No variation (1 type)')+
  #ggtitle('For Low Variation (0.4-0.6)') +
  #ggtitle('For High Variation (0.3-0.7)') +
  xlab("Initial G") + ylab("Final T")+
  scale_color_discrete(name = "Level of Variation")+
  #scale_color_discrete(name = "No. of Types")+
  theme_minimal(base_size = 15)+
  scale_x_continuous(breaks = seq(0,1,0.1))
  
  ggplot(df_s[c(1,4,7)],aes(x=Gi,y=FinG))+
  geom_line(size=1,aes(color = factor(Type)))+
  xlab("Initial G") + ylab("Final G")+
  scale_color_discrete(name = "")+
  theme_minimal(base_size = 15)+
  scale_x_continuous(breaks = seq(0,1,0.1))
