#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 15:08:30 2022

@author: tanveen
"""
import numpy as np
import pandas as pd

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


# simulating ODEs for no variation 
def iterations_no_var(b):
  global Gseq, timesteps, u, v, dt, th
  l=np.zeros((len(Gseq),2))
  for g0 in range(len(Gseq)): 

          G0 = Gseq[g0]
          S0 = T0 = (1-G0)/2
          
          for t in range(0,timesteps):
                G_t= G0 + (u*S0 + v*T0 - G0*b*T0)*dt
                
                S_t = S0+ (b*G0*T0 - (0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*S0 - u*S0)*dt;
                T_t = T0 + ((0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*S0 - v*T0)*dt; 
                G0 = G_t
                T0 = T_t
                S0 = S_t
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

tim = 1 # Total time over which the simulations will run
dt = 0.1 # step size of time 
timesteps = int(tim/dt) # total no. of timesteps

# Set the value of other traits here: 
u = 0.1
v = 0.2
th = 0.5

x = []
for i in Bseq:
    y = iterations_no_var(i)
    x.append(y)
 
# Saving values for the simulation to the dataframes
for i in range(1,len(Bseq)+1):
    for j in range(0,len(Gseq)):
        df_G.iat[j,i] = x[i-1][j][0]
        df_T.iat[j,i] = x[i-1][j][1]

df_G.to_csv('G_u_'+str(u)+'_v_'+str(u)+'_th_'+str(th)+'_no_var.csv')
df_T.to_csv('T_u_'+str(u)+'_v_'+str(u)+'_th_'+str(th)+'_no_var.csv')
