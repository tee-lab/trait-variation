#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import time
from numba import njit
from multiprocessing import Pool
from scipy.stats import beta, uniform

import getopt
import sys

''' Code Segment 1 (CS1)
If you want to set the varying trait, trait distribution and level
of variation while running the script in terminal, keep the following
lines of code till the end of Code Segment 1 (CS1).
For running the code through terminal, give the following command:
python all_dist_trait.py -t (insert trait here) -d (insert distribution here) -v (insert level of variation here)
Example: you want to run the code for variation in sapling death rate (u), bimodal distribution and high level of variation, use this:
python all_dist_trait.py -t u -d bimod -v high
Check CS2 for acceptable values for each of these.

Remove or comment the segment you don't want to use'''
 
# Remove the first argument( the filename)
all_args = sys.argv[1:]

try:
   # Gather the arguments
   opts, args = getopt.getopt(all_args, 't:d:v:')
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

'''CS1 ends here. If you don't wish to run the code like this, go to CS2'''

# Copied cust_range() and crange() from stack overflow
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

# simulating ODEs for variation in sapling death rate (mu)
@njit(fastmath=False)
def iterations_u(b):
  global Gseq, n_types, timesteps, u, v, dt, th, p
  l=np.zeros((len(Gseq),2))
  for g0 in range(len(Gseq)): 

          G0 = Gseq[g0]
          S0 = T0 = (1-G0)/2
          Si = S0*p
          Ti = T0*p 
          
          for t in range(0,timesteps):
                death = u*Si
                G_t= G0 + (np.sum(death) + v*np.sum(Ti) - G0*b*np.sum(Ti))*dt
                for i in range(0,n_types):
                    s = Si[i]
                    Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - u[i]*s)*dt;
                    Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - v*Ti[i])*dt; 
                G0 = G_t
                T0 = np.sum(Ti)
          l[g0]=G0, T0
  return l 

# simulating ODEs for variation in tree death rate (nu)
@njit(fastmath=False)
def iterations_v(b):
  global Gseq, n_types, timesteps, u, v, dt, th, p
  l=np.zeros((len(Gseq),2))
  for g0 in range(len(Gseq)): 

          G0 = Gseq[g0]
          S0 = T0 = (1-G0)/2
          Si = S0*p
          Ti = T0*p 
          
          for t in range(0,timesteps):
                G_t= G0 + (u*np.sum(Si) + np.sum(v*Ti) - G0*b*np.sum(Ti))*dt
                for i in range(0,n_types):
                    s = Si[i]
                    Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - u*s)*dt;
                    Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - v[i]*Ti[i])*dt; 
                G0 = G_t
                T0 = np.sum(Ti)
          l[g0]=G0, T0
  return l 

# simulating ODEs for variation in sapling resistance to fire (theta)
@njit(fastmath=False)
def iterations_th(b):
  global Gseq, n_types, timesteps, u, v, dt, th, p
  l=np.zeros((len(Gseq),2))
  for g0 in range(len(Gseq)): 

          G0 = Gseq[g0]
          S0 = T0 = (1-G0)/2
          Si = S0*p
          Ti = T0*p 
          
          for t in range(0,timesteps):
                G_t= G0 + (u*np.sum(Si) + v*np.sum(Ti) - G0*b*np.sum(Ti))*dt
                for i in range(0,n_types):
                    s = Si[i]
                    Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th[i]-G0)/0.005)))*s - u*s)*dt;
                    Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th[i]-G0)/0.005)))*s - v*Ti[i])*dt; 
                G0 = G_t
                T0 = np.sum(Ti)
          l[g0]=G0, T0
  return l 

Gseq = np.arange(0,1,0.1) # Initial Grass covers
Gseq = np.round(Gseq,2)
Bseq = crange(0,2,0.1) # Values of sapling birth rate 
Bseq = np.round(Bseq,2)

# Dataframes to save the steady-state values of G and T
df_G = pd.DataFrame({'Gi':Gseq})
for i in Bseq:
    df_G[i] = 0
    
df_T = pd.DataFrame({'Gi':Gseq})
for i in Bseq:
    df_T[i] = 0

''' Code Segment 2 (CS2)
Here you can set the values in the code itself
To uncomment the code, remove the # at the start of each line of code'''

#var_list = ["high","low"] # Levels of variation
#dist_list = ["unif","beta","bimod"] # Distribution of traits
#trait_list = ["u","v","th"] 
'''Trait being varied, u = sapling death rate (mu),
v = tree death rate (nu), th = sapling resistance to fire (theta)'''

'''set the choice of trait to vary, distribution of traits and level of variation'''
#trait = trait_list[0] # 0 refers to the first entry in the list
#dist = dist_list[0]
#var = var_list[0] 

''' End of CS2 
Comment or remove this part if using CS1'''

n_types = 10 #Number of sapling (or tree) types, referred to as 'k' in Methods
print('Number of types: ',n_types)

if var == "high": # Range of values for the trait
    min_var = 0 
    max_var = 1
elif var == "low": # Range of values for the trait
    min_var = 0.3
    max_var = 0.7
print("Level of variation:",var)

# for choosing equidistant trait values in the given range
diff = max_var - min_var
binsize = (diff)/(n_types)

p =  np.empty(n_types) # list with proportion of each type
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

tim = 1 # Total time over which the simulations will run
dt = 0.1 # step size of time 
timesteps = int(tim/dt) # total no. of timesteps

if trait=="v":
    print('Varying trait is tree death rate')
    v = tr
    v = np.round(v,3)
    u = 0.05
    th = 0.5
    
    if __name__ == '__main__':
        pool = Pool() 
        start = time.time()
        # We'll run the function iterations_v() over different values of 
        # sapling birth rate in parallel using multiprocessing
        x = pool.map(iterations_v, Bseq) 
        print("Complete")
        end = time.time()
        print('total time (s)= ' + str(end-start))
        
if trait=="u":
     print('Varying trait is sapling death rate')
     u = tr
     u = np.round(u,3)
     v = 0.1
     th = 0.5
     
     if __name__ == '__main__':
         pool = Pool() 
         start = time.time()
         # We'll run the function iterations_u() over different values of 
         # sapling birth rate in parallel using multiprocessing
         x = pool.map(iterations_u, Bseq)
         print("Complete")
         end = time.time()
         print('total time (s)= ' + str(end-start))
         
if trait=="th":
     print('Varying trait is sapling resistance to fire')
     th = tr
     th = np.round(th,3)
     u = 0.2
     v = 0.1
     if __name__ == '__main__':
         pool = Pool() 
         start = time.time()
         # We'll run the function iterations_th() over different values of 
         # sapling birth rate in parallel using multiprocessing
         x = pool.map(iterations_th, Bseq)
         print("Complete")
         end = time.time()
         print('total time (s)= ' + str(end-start))
 
# Saving values for the simulation to the dataframes
for i in range(1,len(Bseq)+1):
    for j in range(0,len(Gseq)):
        df_G.iat[j,i] = x[i-1][j][0]
        df_T.iat[j,i] = x[i-1][j][1]

df_G.to_csv('G_'+str(n_types)+'_types_'+str(dist)+'_dist_varying_'+str(trait)+'_'+str(var)+'_var.csv')
df_T.to_csv('T_'+str(n_types)+'_types_'+str(dist)+'_dist_varying_'+str(trait)+'_'+str(var)+'_var.csv')
