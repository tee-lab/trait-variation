#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import time
from numba import njit
from multiprocessing import Pool
from scipy.stats import beta, uniform

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

# The following functions encode for the ODEs using Euler method for different
# varying traits: sapling death rate (u), tree death rate (v) and sapling resistance to fire (th), respectively.

@njit()
def iterations_u(b):
  global Gseq, d, timesteps, u, v, dt, th, p
  l=np.zeros((len(Gseq),2))
  for g0 in range(len(Gseq)): 

          G0 = Gseq[g0]  # Intial G value
          S0 = T0 = (1-G0)/2 # Initial total S and T value
          # Multiplying with p gives the proportion cover of each trait type
          Si = S0*p 
          Ti = T0*p 
          
          # How the proportion covers evolve for 10^7 timesteps: 
          for t in range(0,timesteps):
                death = u*Si
                G_t= G0 + (np.sum(death) + v*np.sum(Ti) - G0*b*np.sum(Ti))*dt
                for i in range(0,d):
                    s = Si[i]
                    Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - u[i]*s)*dt;
                    Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - v*Ti[i])*dt; 
                G0 = G_t
                T0 = np.sum(Ti)
          l[g0]=G0, T0
  return l 

@njit()
def iterations_v(b):
  global Gseq, d, timesteps, u, v, dt, th, p
  l=np.zeros((len(Gseq),2))
  for g0 in range(len(Gseq)): 

          G0 = Gseq[g0]
          S0 = T0 = (1-G0)/2
          Si = S0*p
          Ti = T0*p 
          
          for t in range(0,timesteps):
                G_t= G0 + (u*np.sum(Si) + np.sum(v*Ti) - G0*b*np.sum(Ti))*dt
                for i in range(0,d):
                    s = Si[i]
                    Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - u*s)*dt;
                    Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - v[i]*Ti[i])*dt; 
                G0 = G_t
                T0 = np.sum(Ti)
          l[g0]=G0, T0
  return l 

@njit()
def iterations_th(b):
  global Gseq, d, timesteps, u, v, dt, th, p
  l=np.zeros((len(Gseq),2))
  for g0 in range(len(Gseq)): 

          G0 = Gseq[g0]
          S0 = T0 = (1-G0)/2
          Si = S0*p
          Ti = T0*p 
          
          for t in range(0,timesteps):
                G_t= G0 + (u*np.sum(Si) + v*np.sum(Ti) - G0*b*np.sum(Ti))*dt
                for i in range(0,d):
                    s = Si[i]
                    Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th[i]-G0)/0.005)))*s - u*s)*dt;
                    Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th[i]-G0)/0.005)))*s - v*Ti[i])*dt; 
                G0 = G_t
                T0 = np.sum(Ti)
          l[g0]=G0, T0
  return l 

Gseq = np.arange(0,1,0.1) # Initial values of G
Gseq = np.round(Gseq,2)
Bseq = crange(0,2,0.01) # Values of sapling birth rate (beta)
Bseq = np.round(Bseq,2)

df_th = pd.DataFrame({'b':Bseq}) # Dataframe to save final grass cover (G) for all values of beta 
df_th_T = pd.DataFrame({'b':Bseq}) # Dataframe to save final tree cover (T) for all values of beta

para_range = ["high","low"] # Two levels of variation

for xx in para_range:
    
    d= 10  # No. of trait types 
    if xx == "high": # High variation corresponds to the range (0,1)
        min_var = 0
        max_var = 1
    elif xx == "low": # Low variation corresponds to the range (0.3,0.7)
        min_var = 0.3
        max_var = 0.7
    print("range:",xx)    
    diff = max_var - min_var
    binsize = (diff)/(d)
    tr = crange(min_var+(binsize/2),max_var-(binsize/2),binsize) # Set of trait values
    
    tim = 10000
    dt = 0.001  # size of each time step
    timesteps = int(tim/dt)
    
    # Create dataframes to save the values of final G and T intermediately
    ar_th = np.empty([3,len(Bseq)+1])
    ar_th_T = np.empty([3,len(Bseq)+1])
    ar_th[0][1:(len(Bseq)+2)]=Bseq
    ar_th_T[0][1:(len(Bseq)+2)]=Bseq
    # We'll collect only the values at G0 = 0 and G0 = 0.9, because 
    # if there's bistability in the system, these two initial values will capture that
    ar_th[1][0]=0
    ar_th_T[1][0]=0
    ar_th[2][0]=0.9
    ar_th_T[2][0]=0.9
    
    
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
    
    # Set the varying trait with 'trait=', if not running from the terminal directly
    if trait=="v":
        v = crange(min_var+(binsize/2),max_var-(binsize/2),binsize)
        v = np.round(v,3)
        u = 0.05
        th = 0.5
        
        #This runs iteration_v() function for different values of beta in parallel:
        if __name__ == '__main__':
            pool = Pool() 
            start = time.time()
            x = pool.map(iterations_v, Bseq)
            print("Complete")
            end = time.time()
            print('total time (s)= ' + str(end-start))
     
     
    if trait=="u":
         u = crange(min_var+(binsize/2),max_var-(binsize/2),binsize)
         u = np.round(u,3)
         v = 0.1
         th = 0.5
         
         #This runs iteration_u() function for different values of beta in parallel:
         if __name__ == '__main__':
             pool = Pool() 
             start = time.time()
             x = pool.map(iterations_u, Bseq)
             print("Complete")
             end = time.time()
             print('total time (s)= ' + str(end-start))
             
    if trait=="th":
         th = crange(min_var+(binsize/2),max_var-(binsize/2),binsize)
         th = np.round(th,3)
         u = 0.2
         v = 0.1
         
         #This runs iteration_th() function for different values of beta in parallel:
         if __name__ == '__main__':
             pool = Pool() 
             start = time.time()
             x = pool.map(iterations_th, Bseq)
             print("Complete")
             end = time.time()
             print('total time (s)= ' + str(end-start))
             
    for i in range(len(Bseq)):
        ar_th[1][i+1]=x[i][0][0] # Collecting final G values for G0 = 0
        ar_th[2][i+1]=x[i][9][0] # Collecting final G values for G0 = 0.9
        
        ar_th_T[1][i+1]=x[i][0][1] # Collecting final T values for G0 = 0
        ar_th_T[2][i+1]=x[i][9][1] # Collecting final G values for G0 = 0.9
    
    # Cleaning the dataframe
    ar_th = np.round(ar_th,3)
    ar_th = np.delete(ar_th,0,0)
    ar_th = np.delete(ar_th,0,1)
    
    ar_th_T =np.round(ar_th_T,2)
    ar_th_T = np.delete(ar_th_T,0,0)
    ar_th_T = np.delete(ar_th_T,0,1)
   
    str1 = xx +"_1"
    str2 = xx +"_2"
    
    # Storing obtained values to the main dataframe
    df_th[str1] = ar_th[0]
    df_th[str2] = ar_th[1]
    
    df_th_T[str1] = ar_th_T[0]
    df_th_T[str2] = ar_th_T[1]

# Saving the dataframes to a CSV for further analysis:
df_th.to_csv('10ty_Si_Ti_'+dist+'_vary_'+trait+'.csv') 
df_th_T.to_csv('T_10ty_Si_Ti_'+dist+'_vary_'+trait+'.csv')    
