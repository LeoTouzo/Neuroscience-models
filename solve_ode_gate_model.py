# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 12:20:23 2021

@author: leoto
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import erf
import seaborn as sns

#w=(wm,gm,wc,gc)

def f(t,w,mu_m,mu_c,sigma_m,sigma_c,l):
    fw=np.zeros(4)
    fw[0]=-(mu_m**2+sigma_m**2)*w[1]**2*w[0]-mu_m*mu_c*w[1]*w[2]*w[3]+mu_m*w[1]
    fw[1]=-(mu_m**2+sigma_m**2)*w[0]**2*w[1]-mu_m*mu_c*w[0]*w[2]*w[3]+mu_m*w[0]-l*np.sign(w[1])
    fw[2]=-(mu_c**2+sigma_c**2)*w[3]**2*w[2]-mu_m*mu_c*w[0]*w[1]*w[3]+mu_c*w[3]
    fw[3]=-(mu_c**2+sigma_c**2)*w[2]**2*w[3]-mu_m*mu_c*w[0]*w[1]*w[2]+mu_c*w[2]-l*np.sign(w[3])
    return fw

def accuracy(w,mu_m,mu_c,sigma_m,sigma_c,sigma):
    sigma_tot2=sigma**2+(w[0]*w[1]*sigma_m)**2+(w[2]*w[3]*sigma_c)**2
    return (1+erf((w[0]*w[1]*mu_m+w[2]*w[3]*mu_c)/np.sqrt(2*sigma_tot2)))/2

mu_m,mu_c=1,1
sigma_m,sigma_c,sigma=5,0.5,0.5
l=0.1
w0=np.array([0.05,0.5,0.05,0.5])
t_switch=20
t_end=100
dt=0.05
n_steps=int(t_end/dt)

t=np.linspace(0,t_end,n_steps+1)
w=np.zeros((4,n_steps+1))
sol=solve_ivp(f,(0,t_switch),w0,t_eval=np.linspace(0,t_switch,int(t_switch/dt)+1),args=(mu_m,0,sigma_m,sigma_c,l))
w[:,:int(t_switch/dt)+1]=sol.y

t2=t_end-t_switch
sol=solve_ivp(f,(0,t2),w[:,int(t_switch/dt)],t_eval=np.linspace(0,t2,int(t2/dt)+1),args=(mu_m,mu_c,sigma_m,sigma_c,l))
w[:,int(t_switch/dt):]=sol.y

acc=np.zeros(n_steps)
acc[:int(t_switch/dt)]=accuracy(w[:,:int(t_switch/dt)],mu_m,0,sigma_m,sigma_c,sigma)
acc[int(t_switch/dt):]=accuracy(w[:,int(t_switch/dt):],mu_m,mu_c,sigma_m,sigma_c,sigma)

fig,ax = plt.subplots()
ax.plot(t,w[0],label='w_m')
ax.plot(t,w[1],label='g_m')
ax.plot(t,w[2],label='w_c')
ax.plot(t,w[3],label='g_c')
plt.axvline(t_switch,linestyle='dashed',color='red')
ax2=ax.twinx()
ax2.plot(t,acc,color='black',label='accuracy')
ax.set_xlabel('t')
ax.set_ylabel('w')
ax2.set_ylabel('accuracy')
ax.set_xticks([0,10,20,30,40])
ax.set_xticklabels([0,10,'t_switch=20',30,40])
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1+h2, l1+l2)