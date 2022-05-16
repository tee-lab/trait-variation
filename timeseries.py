#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import math
import pandas as pd
import time
from scipy.stats import beta, uniform
from numba import njit, prange


# If one wishes to run this file from the terminal, one can specify the varying trait 
# and the distribution of the trait as arguments. For example: python filename -t u -d unif
# Here -t refers to the trait and -d refers to the distribution of the trait
import getopt
import sys

# Remove the first argument( the filename)
all_args = sys.argv[1:]

try:
   # Gather the arguments
   opts, args = getopt.getopt(all_args, 't:d:')
except:
        print("Error")
        
print("opts:", opts)
print("args:", args) 
   
for opt, arg in opts:
    if opt in ['-t']:
        trait = str(arg)
    elif opt in ['-d']:
        dist = str(arg)

# The following code defines a function crange, that generates a sequence 
# while making sure that the start and end points are included in it
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

d= 10 # No. of trait types 
# Range of trait values: (0,1)
min_var = 0
max_var = 1
diff = max_var - min_var
binsize = (diff)/(d)
tr = crange(min_var+(binsize/2),max_var-(binsize/2),binsize) # Set of trait values

b= 0.45 # Setting sapling birth rate = 0.45 corresponds to the woodland regime
G0 = 0.5 # Intial G value
S0 = T0 = (1-G0)/2 # Initial total S and T value

tim = 200
dt = 0.001
timesteps = int(tim/dt)

# Dataframe to store grass cover, individual tree cover and sapling cover of each trait type at every timestep
df_time = pd.DataFrame({'Time':range(0,timesteps),'G':[0.0]*timesteps,'S1': [0.0]*timesteps, 'S2':[0.0]*timesteps,
                      'S3':[0.0]*timesteps,'S4':[0.0]*timesteps,'S5':[0.0]*timesteps,'S6': [0.0]*timesteps, 
                      'S7':[0.0]*timesteps,'S8':[0.0]*timesteps,'S9':[0.0]*timesteps,'S10':[0.0]*timesteps,
                      'T1': [0.0]*timesteps, 'T2':[0.0]*timesteps,'T3':[0.0]*timesteps,'T4':[0.0]*timesteps,
                      'T5':[0.0]*timesteps,'T6': [0.0]*timesteps,'T7':[0.0]*timesteps,'T8':[0.0]*timesteps,
                      'T9':[0.0]*timesteps,'T10':[0.0]*timesteps})

# Set the Distribution for the trait with 'dist=', if not running from the terminal directly
# p[i] tells the proportion of each trait value out of 1
p =  np.empty(d)
sump=0
if dist=="bimod": # This refers to bimodal beta distribution
    for i in tr:
        sump += beta.pdf(i,2,8)/2 + beta.pdf(i,8,2)/2
    for i in range(0,d,1):
        p[i] = (beta.pdf(tr[i],2,8)/2 + beta.pdf(tr[i],8,2)/2)/sump  
if dist=="beta": # This refers to unimodal beta distribution
    for i in tr:
        sump += beta.pdf(i,4,4)
    for i in range(0,d,1):
        p[i] = beta.pdf(tr[i],4,4)/sump
if dist=="unif": # This refers to uniform distribution
    for i in tr:
        sump += uniform.pdf(i,min_var,max_var)
    for i in range(0,d,1):
        p[i] = uniform.pdf(tr[i],min_var,max_var)/sump

# Multiplying with p gives the proportion cover of each trait type
Si = S0*p
Ti = T0*p 

if trait=="v":
    v = crange(min_var+(binsize/2),max_var-(binsize/2),binsize) # Set of values of v
    v = np.round(v,3)
    u = 0.05
    th = 0.5
    
    start = time.time()
    for t in range(0,timesteps):
        if t!=0:
          G_t= G0 + (u*np.sum(Si) + np.sum(v*Ti) - G0*b*np.sum(Ti))*dt
          for i in range(0,d):
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
     u = crange(min_var+(binsize/2),max_var-(binsize/2),binsize) # Set of values of u
     u = np.round(u,3)
     v = 0.1
     th = 0.5
     
     start = time.time()
     for t in range(0,timesteps):
         if t!=0:
            death = u*Si
            G_t= G0 + (np.sum(death) + v*np.sum(Ti) - G0*b*np.sum(Ti))*dt
            for i in range(0,d):
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
     th = crange(min_var+(binsize/2),max_var-(binsize/2),binsize) # Set of values of theta
     th = np.round(th,3)
     u = 0.2
     v = 0.1
     start = time.time()
     for t in range(0,timesteps):
         if t!=0:
           G_t= G0 + (u*np.sum(Si) + v*np.sum(Ti) - G0*b*np.sum(Ti))*dt
           for i in range(0,d):
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

# Saving the dataframes to a CSV for further analysis:
df_time.to_csv('TimeSeries_full_woodland_'+trait+'_'+dist+'.csv')
