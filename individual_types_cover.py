

import numpy as np
import math
import matplotlib.pyplot as plt
import statistics as st
import pandas as pd
from scipy.stats import uniform

def cust_range(*args, rtol=1e-05, atol=1e-08, include=[True, False]):
    """
    Combines numpy.arange and numpy.isclose to mimic
    open, half-open and closed intervals.
    Avoids also floating point rounding errors as with
    >>> numpy.arange(1, 1.3, 0.1)
    array([1. , 1.1, 1.2, 1.3])

    args: [start, ]stop, [step, ]
        as in numpy.arange
    rtol, atol: floats
        floating point tolerance as in numpy.isclose
    include: boolean list-like, length 2
        if start and end point are included
    """
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

def omega(G,th):  
    return 0.9 + (0.05-0.9)/(1+math.exp((th-G)/0.005))

Gseq = np.arange(0,1.05,0.1) 
df_th = pd.DataFrame({'Gi':Gseq,'T1': [0.0]*len(Gseq), 'T2':[0.0]*len(Gseq),'T3':[0.0]*len(Gseq),'T4':[0.0]*len(Gseq),'T5':[0.0]*len(Gseq)})
    #, 'T6': [0.0]*len(Gseq), 'T7':[0.0]*len(Gseq),'T8':[0.0]*len(Gseq),
     #                 'T9':[0.0]*len(Gseq),'T10':[0.0]*len(Gseq)})

for r in Gseq: 
    bins = 5
    min_var = 0.25
    max_var = 0.75
    #min_var = 0.35
    #max_var = 0.65
    th = list(np.arange(min_var, max_var,0.0001))
    x = max_var - min_var
    binsize = (x-0.1)/(bins-1)
    a =  np.empty(bins)
    l = uniform.pdf(th)
    q= list(crange(0.3,0.7,binsize))
    q  = [round(i,1) for i in q]
    v = range(0,bins,1)
    for i in v:
        a[i]=st.mean(l[int(((i)*((len(th)-1)/bins))):int(((i+1)*((len(th)-1)/bins)+1))])
        a[i] = (a[i]*uniform.pdf(q[i],min_var,max_var))
    a
    p = a/sum(a)
    x = 10000
    pr = p*x
    a_beta1 = np.empty(0)
    for i in range(0,len(pr)):
        ab = np.array([q[i]]*int(pr[i]))
        a_beta1 = np.append(a_beta1,ab)
    d=len(a_beta1)
    Gi = G0 = r
    S0 = T0 = (1-G0)/2
    b, u, v = 0.35, 0.1, 0.1
    S_t = Si = [S0/d for i in range(0,d)]
    T_t = Ti = [T0/d for i in range(0,d)]
    G = [G0 for i in range(0,d)]
    tim = 200
    dt = 0.001
    timesteps = int(tim/dt)
    df = pd.DataFrame(a_beta1, columns=['th'])
    df['FinS'] = Si
    df['FinT'] = Ti
    tic1 = time.process_time()
    for t in range(0,timesteps):
        G_t=G0+(u*S0 + v*T0 - b*G0*T0)*dt
        for i in q:
            z1=df[df['th']==i].index.tolist()
            z1[0]
            S_i = df.loc[z1[0],'FinS'] + (b*G0*df.loc[z1[0],'FinT'] - omega(G0,i)*df.loc[z1[0],'FinS'] - u*df.loc[z1[0],'FinS'])*dt;
            T_i = df.loc[z1[0],'FinT']+ (omega(G0,i)*df.loc[z1[0],'FinS'] - v*df.loc[z1[0],'FinT'])*dt; 
        
            df.loc[z1,'FinS']= S_i
            df.loc[z1,'FinT'] = T_i
            
        G0 = G_t
        S0 = sum(df.loc[:,'FinS'])
        T0 = sum(df.loc[:,'FinT'])
    toc1 = time.process_time()
    print(toc1 - tic1)
#    z = int(r*10)
#    df_th.loc[z,'T1']= sum(df.loc[df['th']==q[0],'FinT'])
#    df_th.loc[[z],['T2']]= sum(df.loc[df['th']==q[1],'FinT'])
#    df_th.loc[[z],['T3']]= sum(df.loc[df['th']==q[2],'FinT'])
#    df_th.loc[[z],['T4']]= sum(df.loc[df['th']==q[3],'FinT'])
#    df_th.loc[[z],['T5']]= sum(df.loc[df['th']==q[4],'FinT'])
#    df_th.loc[z,'T6']= sum(df.loc[df['th']==q[5],'FinT'])
#    df_th.loc[[z],['T7']]= sum(df.loc[df['th']==q[6],'FinT'])
#    df_th.loc[[z],['T8']]= sum(df.loc[df['th']==q[7],'FinT'])
#    df_th.loc[[z],['T9']]= sum(df.loc[df['th']==q[8],'FinT'])
#    df_th.loc[[z],['T10']]= sum(df.loc[df['th']==q[9],'FinT'])
    z=df_th[df_th['Gi']==r].index.tolist()
    df_th.loc[z[0],'T1']= sum(df.loc[df['th']==q[0],'FinT'])
    df_th.loc[[z[0]],['T2']]= sum(df.loc[df['th']==q[1],'FinT'])
    df_th.loc[[z[0]],['T3']]= sum(df.loc[df['th']==q[2],'FinT'])
    df_th.loc[[z[0]],['T4']]= sum(df.loc[df['th']==q[3],'FinT'])
    df_th.loc[[z[0]],['T5']]= sum(df.loc[df['th']==q[4],'FinT'])
#    df_th.loc[z[0],'T6']= sum(df.loc[df['th']==q[5],'FinT'])
#    df_th.loc[[z[0]],['T7']]= sum(df.loc[df['th']==q[6],'FinT'])
#    df_th.loc[[z[0]],['T8']]= sum(df.loc[df['th']==q[7],'FinT'])
#    df_th.loc[[z[0]],['T9']]= sum(df.loc[df['th']==q[8],'FinT'])
#    df_th.loc[[z[0]],['T10']]= sum(df.loc[df['th']==q[9],'FinT'])

q2 = [round(x,2) for x in q if type(x)==np.float64]
q2.insert(0,'Gi')
df_th.columns = q2

f, ax = plt.subplots()    
df_th.plot(x="Gi", kind="bar", stacked=True, cmap="plasma")
