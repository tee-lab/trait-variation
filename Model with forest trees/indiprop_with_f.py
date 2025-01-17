# -*- coding: utf-8 -*-
"""indiprop_with_F.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1nrRGmTfdnxs68HslzoJO2lcWoEMTKVhJ
"""

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
# This function returns the proportion of each tree and sapling type at steady-state
@njit(fastmath=False)
def props_u(g):
    global n_types, timesteps, u, v, dt, th, b, p
    G0 = g
    F0 = 0.0001
    S0 = T0 = (1-G0-F0)/2
    Si = S0*p
    Ti = T0*p

    for t in range(0,timesteps):
          death = u*Si
          G_t= G0 + (np.sum(death) + v*np.sum(Ti) - G0*b*np.sum(Ti)+(0.05 + (0.9-0.05)/(1+np.exp((th-G0)/0.005)))*F0-al*G0*F0)*dt
          for i in range(0,n_types):
              s = Si[i]
              Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - u[i]*s-al*s*F0)*dt;
              Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - v*Ti[i]-al*Ti[i]*F0)*dt;
              Si[i] = max(min(Si[i], 1), 0)
              Ti[i] = max(min(Ti[i], 1), 0)
          F_t = F0 + (al*(1-F0)*F0-(0.05 + (0.9-0.05)/(1+np.exp((th-G0)/0.005)))*F0)*dt
          F0 = max(min(F_t, 1), 0)
          G0 = G_t
          S0 = np.sum(Si)
          T0 = np.sum(Ti)

    return Si, Ti

# simulating ODEs for variation in tree death rate (nu)
# This function returns the proportion of each tree and sapling type at steady-state
@njit(fastmath=False)
def props_v(g):
    global n_types, timesteps, u, v, dt, th, b, p
    G0 = g
    F0 = 0.0001
    S0 = T0 = (1-G0-F0)/2
    Si = S0*p
    Ti = T0*p

    for t in range(0,timesteps):
          G_t= G0 + (u*np.sum(Si) + np.sum(v*Ti) - G0*b*np.sum(Ti)+(0.05 + (0.9-0.05)/(1+np.exp((th-G0)/0.005)))*F0-al*G0*F0)*dt
          for i in range(0,n_types):
                 s = Si[i]
                 Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - u*s-al*s*F0)*dt;
                 Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - v[i]*Ti[i]-al*Ti[i]*F0)*dt;
                 Si[i] = max(min(Si[i], 1), 0)
                 Ti[i] = max(min(Ti[i], 1), 0)
          F_t = F0 + (al*(1-F0)*F0-(0.05 + (0.9-0.05)/(1+np.exp((th-G0)/0.005)))*F0)*dt
          F0 = max(min(F_t, 1), 0)
          G0 = max(min(G_t, 1), 0)
          S0 = np.sum(Si)
          T0 = np.sum(Ti)
    return Si, Ti


# simulating ODEs for variation in sapling resistance to fire (theta)
# This function returns the proportion of each tree and sapling type at steady-state
@njit(fastmath=False)
def props_th(g):
    global n_types, timesteps, u, v, dt, th, b, p
    G0 = g
    F0 = 0.0001
    S0 = T0 = (1-G0-F0)/2
    Si = S0*p
    Ti = T0*p

    for t in range(0,timesteps):
            G_t= G0 + (u*np.sum(Si) + v*np.sum(Ti) - G0*b*np.sum(Ti)+(0.05 + (0.9-0.05)/(1+np.exp((0.5-G0)/0.005)))*F0-al*G0*F0)*dt
            for i in range(0,n_types):
                s = Si[i]
                Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th[i]-G0)/0.005)))*s - u*s-al*s*F0)*dt;
                Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th[i]-G0)/0.005)))*s - v*Ti[i]-al*Ti[i]*F0)*dt;
                Si[i] = max(min(Si[i], 1), 0)
                Ti[i] = max(min(Ti[i], 1), 0)
            F_t = F0 + (al*(1-F0)*F0-(0.05 + (0.9-0.05)/(1+np.exp((0.5-G0)/0.005)))*F0)*dt
            F0 = max(min(F_t, 1), 0)
            G0 = max(min(G_t, 1), 0)
            T0 = np.sum(Ti)
            S0 = np.sum(Si)
    return Si, Ti


Gseq = np.arange(0,1,0.01) # Initial Grass covers
Gseq = np.round(Gseq,2)

''' Code Segment 2 (CS2)
Here you can set the values in the code itself'''

# Based on the bifurcation diagram, choose the value of b. This would determine
# in which regime do you want to see the final trait distribution
b = 0.9 # Value of sapling birth rate (beta)
al = 0.03

var_list = ["high","low"] # Levels of variation
dist_list = ["unif","beta","bimod"] # Distribution of traits
trait_list = ["u","v","th"]
'''Trait being varied, u = sapling death rate (mu),
v = tree death rate (nu), th = sapling resistance to fire (theta)'''

#set the choice of trait to vary, distribution of traits and level of variation
trait = trait_list[1] # 0 refers to the first entry in the list
dist = dist_list[0]
var = var_list[1]

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

tim = 100000 # Total time over which the simulations will run
dt = 0.1 # step size of time
timesteps = int(tim/dt) # total no. of timesteps

# Dataframes to save the steady-state proportion of all Sapling and Tree types
df_S = pd.DataFrame()
df_S['trait_value'] = tr
df_S['initial_prop'] = p
df_T = pd.DataFrame()
df_T['trait_value'] = tr
df_T['initial_prop'] = p

print('Sapling birth rate = ', b)
if trait=="v":
    print('Varying trait is tree death rate')
    v = tr
    v = np.round(v,3)
    u = 0.02
    th = 0.5

    if __name__ == '__main__':
        pool = Pool()
        start = time.time()
        # We'll run the function props_v() over different
        # initial values of G in parallel using multiprocessing
        print("Start of simulation")
        x = pool.map(props_v, Gseq)
        print("End of simulation")
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
         # We'll run the function props_u() over different
         # initial values of G in parallel using multiprocessing
         print("Start of simulation")
         x = pool.map(props_u,Gseq)
         print("End of simulation")
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
         # We'll run the function props_th() over different
         # initial values of G in parallel using multiprocessing
         print("Start of simulation")
         x = pool.map(props_th, Gseq)
         print("End of simulation")
         end = time.time()
         print('total time (s)= ' + str(end-start))

# Saving values for the simulation to the dataframes
for i in range(len(Gseq)):
    df_S[i/10]=x[i][0]
    df_T[i/10]=x[i][1]

df_S.to_csv('S_indi_prop_'+str(n_types)+'_types_'+str(dist)+'_dist_varying_'+str(trait)+'_'+str(var)+'_var_beta_'+str(b)+'.csv')
df_T.to_csv('T_indi_prop_'+str(n_types)+'_types_'+str(dist)+'_dist_varying_'+str(trait)+'_'+str(var)+'_var_beta_'+str(b)+'.csv')