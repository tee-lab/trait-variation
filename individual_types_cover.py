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

@njit(fastmath=False)
def props_u(g):
    global d, timesteps, u, v, dt, th, b, p
    G0 = g # Intial G value
    S0 = T0 = (1-G0)/2  # Initial total S and T value
    # Multiplying with p gives the proportion cover of each trait type
    Si = S0*p
    Ti = T0*p 
  
    for t in range(0,timesteps):
          death = u*Si
          G_t= G0 + (np.sum(death) + v*np.sum(Ti) - G0*b*np.sum(Ti))*dt
          for i in range(0,d):
              s = Si[i]
              Si[i] += (b*G0*Ti[i] - (0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - u[i]*s)*dt;
              Ti[i] += ((0.9 + (0.05-0.9)/(1+np.exp((th-G0)/0.005)))*s - v*Ti[i])*dt; 
          G0 = G_t
          S0 = np.sum(Si)
          T0 = np.sum(Ti)

    return Si, Ti 


@njit(fastmath=False)
def props_v(g):
    global d, timesteps, u, v, dt, th, b, p
    G0 = g
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
          S0 = np.sum(Si)
          T0 = np.sum(Ti)

    return Si, Ti 

@njit(fastmath=False)
def props_th(g):
    global d, timesteps, u, v, dt, th, b, p
    G0 = g
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
          S0 = np.sum(Si)
          T0 = np.sum(Ti)

    return Si, Ti 

Gseq = np.arange(0,1,0.1) # Initial values of G
Gseq = np.round(Gseq,2)

df_S_Wo = pd.DataFrame() # Dataframe to store individual sapling cover of each trait type 
df_T_Wo = pd.DataFrame() # Dataframe to store individual tree cover of each trait type 

d = 10 # No. of trait types 
# Range of trait values: (0,1)
min_var = 0
max_var = 1
b = 0.45 # Setting sapling birth rate = 0.45 corresponds to the woodland regime
print("b:",b)    
diff = max_var - min_var
binsize = (diff)/(d)
tr = crange(min_var+(binsize/2),max_var-(binsize/2),binsize) # Set of trait values


tim = 10000
dt = 0.001 # size of each time step
timesteps = int(tim/dt)

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
    
    #This runs props_v() function for different values of G0 in parallel:
    if __name__ == '__main__':
        pool = Pool() 
        start = time.time()
        x = pool.map(props_v, Gseq)
        print("Complete")
        end = time.time()
        print('total time (s)= ' + str(end-start))
        
if trait=="u":
     u = crange(min_var+(binsize/2),max_var-(binsize/2),binsize)
     u = np.round(u,3)
     v = 0.1
     th = 0.5
     
     #This runs props_u() function for different values of G0 in parallel:
     if __name__ == '__main__':
         pool = Pool() 
         start = time.time()
         x = pool.map(props_u,Gseq)
         print("Complete")
         end = time.time()
         print('total time (s)= ' + str(end-start))
         
if trait=="th":
     th = crange(min_var+(binsize/2),max_var-(binsize/2),binsize)
     th = np.round(th,3)
     u = 0.2
     v = 0.1
    
     #This runs props_th() function for different values of G0 in parallel:
     if __name__ == '__main__':
         pool = Pool() 
         start = time.time()
         x = pool.map(props_th, Gseq)
         print("Complete")
         end = time.time()
         print('total time (s)= ' + str(end-start))

# Storing obtained individual sapling and tree covers to the dataframes
for i in range(len(Gseq)):
    df_S_Wo[i/10]=x[i][0] 
    df_T_Wo[i/10]=x[i][1] 

# Saving the dataframes to a CSV for further analysis:
df_S_Wo.to_csv('S_Indi_prop_'+trait+'_'+dist+'.csv')
df_T_Wo.to_csv('T_Indi_prop_'+trait+'_'+dist+'.csv')
