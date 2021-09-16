# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 11:22:36 2021

@author: leoto
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import erf
import seaborn as sns

def f_noise(w,mu_m,mu_c,sigma_m,sigma_c,sigma_ksi,l):
    fw=np.zeros(2)
    ksi=np.random.normal(0,sigma_ksi,2)
    fw[0]=-(mu_m**2+sigma_m**2)*w[0]**3-mu_m*mu_c*w[0]*w[1]**2+(mu_m-l)*w[0]
    fw[1]=-(mu_c**2+sigma_c**2)*w[1]**3-mu_m*mu_c*w[0]**2*w[1]+(mu_c-l)*w[1]
    return fw+ksi

def accuracy(w,mu_m,mu_c,sigma_m,sigma_c,sigma):
    sigma_tot2=sigma**2+w[0]**4*sigma_m**2+w[1]**4*sigma_c**2
    return (1+erf((w[0]**2*mu_m+w[1]**2*mu_c)/np.sqrt(2*sigma_tot2)))/2

def integrate(w0,mu_m,mu_c,sigma_m,sigma_c,sigma_ksi,l,dt,t_end):
    nsteps=int(t_end/dt)
    w=np.zeros((2,nsteps+1))
    w[:,0]=w0
    for i in range(nsteps):
        w[:,i+1]=w[:,i]+dt*f_noise(w[:,i],mu_m,mu_c,sigma_m,sigma_c,l,sigma_ksi)
    return w

mu_m,mu_c=1,1
sigma_m,sigma_c,sigma,sigma_ksi=5,0.5,0.5,0.5
l=0.1
w0=np.ones(2)*1e-1
t_switch=20
t_end=80
dt=0.05
n_steps=int(t_end/dt)

t=np.linspace(0,t_end,n_steps+1)
fig,ax = plt.subplots()
for k in range(10):
    w=np.zeros((2,n_steps+1))
    w[:,:int(t_switch/t_end*n_steps+1)]=integrate(w0,mu_m,0,sigma_m,sigma_c,sigma_ksi,l,dt,t_switch)
    w[:,int(t_switch/t_end*n_steps):]=integrate(w[:,int(t_switch/t_end*n_steps)],mu_m,mu_c,sigma_m,sigma_c,sigma_ksi,l,dt,t_end-t_switch)
    acc=np.zeros(n_steps+1)
    acc[:int(t_switch/t_end*n_steps)]=accuracy(w[:,:int(t_switch/t_end*n_steps)],mu_m,0,sigma_m,sigma_c,sigma)
    acc[int(t_switch/t_end*n_steps):]=accuracy(w[:,int(t_switch/t_end*n_steps):],mu_m,mu_c,sigma_m,sigma_c,sigma)
    ax.plot(t,acc)
plt.axvline(t_switch,linestyle='dashed',color='red')
ax.set_xlabel('t')
ax.set_ylabel('accuracy')
ax.set_xticks([0,20,40,60,80])
ax.set_xticklabels([0,'t_switch=20',40,60,80])
ax.legend()