#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import math
import pandas as pd
import time
from scipy.stats import beta, uniform

import getopt
import sys

''' Code Segment 1 (CS1)
If you want to set the varying trait, trait distribution, level
of variation and sapling birth rate (beta) while running the script in terminal, 
keep the following lines of code till the end of Code Segment 1 (CS1).
For running the code through terminal, give the following command:
python Timeseries_all_dist_trait.py -t (insert trait here) -d (insert distribution here) -v (insert level of variation here) -b (insert sapling birth rate here)
Example: you want to run the code for variation in sapling death rate (u), 
bimodal distribution, high level of variation and for beta = 0.3, use this:
python Timeseries_all_dist_trait.py -t u -d bimod -v high -b 0.3
Check CS2 for acceptable values for each of these.

Remove or comment the segment you don't want to use'''
 
# Remove the first argument( the filename)
all_args = sys.argv[1:]

try:
   # Gather the arguments
   opts, args = getopt.getopt(all_args, 't:d:v:b:')
except:
        print("Error")
        
print("opts:", opts)
print("args:", args) 
   
for opt, arg in opts:
    if opt in ['-t']:
        trait = str(arg)
    elif opt in ['-d']:
        dist = str(arg)
    elif opt in ['-v']:
        var = str(arg)
    elif opt in ['-b']:
        b = float(arg)

'''CS1 ends here. If you don't wish to run the code like this, go to CS2'''


def cust_range(*args, rtol=1e-05, atol=1e-08, include=[True, False]):
    # process arguments
    if len(args) == 1:
        start = 0
        stop = args[0]
        step = 1
    elif len(args) == 2:
        start, stop = args
        step = 1
    else:
        assert len(args) == 3
        start, stop, step = tuple(args)

    # determine number of segments
    n = (stop-start)/step + 1

    # do rounding for n
    if np.isclose(n, np.round(n), rtol=rtol, atol=atol):
        n = np.round(n)

    # correct for start/end is exluded
    if not include[0]:
        n -= 1
        start += step
    if not include[1]:
        n -= 1
        stop -= step

    return np.linspace(start, stop, int(n))

def crange(*args, **kwargs):
    return cust_range(*args, **kwargs, include=[True, True])

def orange(*args, **kwargs):
    return cust_range(*args, **kwargs, include=[True, False])

print('crange(1, 1.3, 0.1) >>>', crange(1, 1.3, 0.1))

# Omega function with G and theta (sapling resistance to fire) as variables
def omega(G,th):  
    return 0.9 + (0.05-0.9)/(1+math.exp((th-G)/0.005))

''' Code Segment 2 (CS2)
Here you can set the values in the code itself'''

# # Based on the bifurcation diagram, choose the value of b. This would determine
# # in which regime do you want to see the timeseries
# b = 0.45 # Value of sapling birth rate (beta)

# var_list = ["high","low"] # Levels of variation
# dist_list = ["unif","beta","bimod"] # Distribution of traits
# trait_list = ["u","v","th"] 
# '''Trait being varied, u = sapling death rate (mu),
# v = tree death rate (nu), th = sapling resistance to fire (theta)'''

# #set the choice of trait to vary, distribution of traits and level of variation
# trait = trait_list[0] # 0 refers to the first entry in the list
# dist = dist_list[0]
# var = var_list[0] 

''' End of CS2 
Comment or remove this part if using CS1'''


# Set the initial value of Grass cover here:
G0 = 0.5
S0 = T0 = (1-G0)/2

tim = 1 # Total time over which the simulations will run
dt = 0.1 # step size of time 
timesteps = int(tim/dt) # total no. of timesteps

# Create a dataframe to save the proportion of each variable (G, all trees and all saplings)
# at every time point
df_time = pd.DataFrame({'Time':range(0,timesteps),'G':[0.0]*timesteps,'S1': [0.0]*timesteps, 'S2':[0.0]*timesteps,
                      'S3':[0.0]*timesteps,'S4':[0.0]*timesteps,'S5':[0.0]*timesteps,'S6': [0.0]*timesteps, 
                      'S7':[0.0]*timesteps,'S8':[0.0]*timesteps,'S9':[0.0]*timesteps,'S10':[0.0]*timesteps,
                      'T1': [0.0]*timesteps, 'T2':[0.0]*timesteps,'T3':[0.0]*timesteps,'T4':[0.0]*timesteps,
                      'T5':[0.0]*timesteps,'T6': [0.0]*timesteps,'T7':[0.0]*timesteps,'T8':[0.0]*timesteps,
                      'T9':[0.0]*timesteps,'T10':[0.0]*timesteps})

n_types = 10 #Number of sapling (or tree) types, referred to as 'k' in Methods
print('Number of types: ',n_types)

if var == "high": # Range of values for the trait
    min_var = 0 
    max_var = 1
elif var == "low": # Range of values for the trait
    min_var = 0.3
    max_var = 0.7
print("Level of variation:",var)

# for choosing equidistant trait values in the given range:
diff = max_var - min_var
binsize = (diff)/(n_types)

p =  np.empty(n_types)  # list with proportion of each type
sump=0

# Values of the traits depending on the Level of variation and No. of types
tr = crange(min_var+(binsize/2),max_var-(binsize/2),binsize)

#Distribution tells the proportion (p) of each trait type wrt other types 

if dist=="bimod":
    print('Trait values have a bimodal beta distribution')
    for i in tr:
        sump += beta.pdf(i,2,8)/2 + beta.pdf(i,8,2)/2
    for i in range(0,n_types,1):
        p[i] = (beta.pdf(tr[i],2,8)/2 + beta.pdf(tr[i],8,2)/2)/sump
if dist=="beta":
    print('Trait values have a unimodal beta distribution')
    for i in tr:
        sump += beta.pdf(i,4,4)
    for i in range(0,n_types,1):
        p[i] = beta.pdf(tr[i],4,4)/sump
if dist=="unif":
    print('Trait values have a uniform distribution')
    for i in tr:
        sump += uniform.pdf(i,min_var,max_var)
    for i in range(0,n_types,1):
        p[i] = uniform.pdf(tr[i],min_var,max_var)/sump

# Si (Ti) gives a list of the proportion of each sapling (tree) type wrt other variables  
Si = S0*p
Ti = T0*p 

if trait=="v":
    print('Varying trait is tree death rate')
    v = tr
    v = np.round(v,3)
    u = 0.05
    th = 0.5
    
    print("Start of simulation")
    start = time.time()
    for t in range(0,timesteps):
        if t!=0:
          G_t= G0 + (u*np.sum(Si) + np.sum(v*Ti) - G0*b*np.sum(Ti))*dt
          for i in range(0,n_types):
              s = Si[i]
              Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - u*s)*dt;
              Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - v[i]*Ti[i])*dt; 
          G0 = G_t
          S0 = np.sum(Si)
          T0 = np.sum(Ti)
                
        z=df_time[df_time['Time']==t].index.tolist()
        df_time.loc[z[0],'G'] = G0
        df_time.iloc[z[0],2:12] = Si
        df_time.iloc[z[0],12:22] = Ti
        
    print("Complete")
    end = time.time()
    print('total time (s)= ' + str(end-start))
        
if trait=="u":
     print('Varying trait is sapling death rate')
     u = tr
     u = np.round(u,3)
     v = 0.1
     th = 0.5
     
     print("Start of simulation")
     start = time.time()
     for t in range(0,timesteps):
         if t!=0:
            death = u*Si
            G_t= G0 + (np.sum(death) + v*np.sum(Ti) - G0*b*np.sum(Ti))*dt
            for i in range(0,n_types):
                s = Si[i]
                Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - u[i]*s)*dt;
                Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - v*Ti[i])*dt; 
            G0 = G_t
            S0 = np.sum(Si)
            T0 = np.sum(Ti)
                 
         z=df_time[df_time['Time']==t].index.tolist()
         df_time.loc[z[0],'G'] = G0
         df_time.iloc[z[0],2:12] = Si
         df_time.iloc[z[0],12:22] = Ti
         
     print("Complete")
     end = time.time()
     print('total time (s)= ' + str(end-start))
        
     
if trait=="th":
     print('Varying trait is sapling resistance to fire')
     th = tr
     th = np.round(th,3)
     u = 0.2
     v = 0.1
     
     print("Start of simulation")
     start = time.time()
     for t in range(0,timesteps):
         if t!=0:
           G_t= G0 + (u*np.sum(Si) + v*np.sum(Ti) - G0*b*np.sum(Ti))*dt
           for i in range(0,n_types):
               s = Si[i]
               Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th[i]-G0)/0.005)))*s - u*s)*dt;
               Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th[i]-G0)/0.005)))*s - v*Ti[i])*dt; 
           G0 = G_t
           S0 = np.sum(Si)
           T0 = np.sum(Ti)
                 
         z=df_time[df_time['Time']==t].index.tolist()
         df_time.loc[z[0],'G'] = G0
         df_time.iloc[z[0],2:12] = Si
         df_time.iloc[z[0],12:22] = Ti
         
     print("Complete")
     end = time.time()
     print('total time (s)= ' + str(end-start))

df_time.to_csv('TimeSeries_'+str(n_types)+'_types_'+str(dist)+'_dist_varying_'+str(trait)+'_'+str(var)+'_var_beta_'+str(b)+'.csv')
