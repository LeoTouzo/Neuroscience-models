# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 19:23:28 2021

@author: leoto
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import erf
import seaborn as sns

def f(t,w,mu_m,mu_c,sigma_m,sigma_c,l):
    fw=np.zeros(2)
    fw[0]=-(mu_m**2+sigma_m**2)*w[0]**3-mu_m*mu_c*w[0]*w[1]**2+(mu_m-l)*w[0]
    fw[1]=-(mu_c**2+sigma_c**2)*w[1]**3-mu_m*mu_c*w[0]**2*w[1]+(mu_c-l)*w[1]
    return fw

def accuracy(w,mu_m,mu_c,sigma_m,sigma_c,sigma):
    sigma_tot2=sigma**2+w[0]**4*sigma_m**2+w[1]**4*sigma_c**2
    return (1+erf((w[0]**2*mu_m+w[1]**2*mu_c)/np.sqrt(2*sigma_tot2)))/2

def m1_sol(t,w0_m,mu_m,sigma_m):
    return 1/np.sqrt((mu_m**2+sigma_m**2)/(mu_m-l)+(1/w0_m**2-(mu_m**2+sigma_m**2)/(mu_m-l))*np.exp(-2*(mu_m-l)*t))

def c1_sol(t,w0_c,sigma_c):
    return 1/np.sqrt(-sigma_c**2/l+(1/w0_c**2+sigma_c**2/l)*np.exp(2*l*t))

mu_m,mu_c=1,1
sigma_m,sigma_c,sigma=5,0.5,0.5
l=0.1
w0=np.ones(2)*1e-1
t_switch=20
t_end=40
n_steps=401
dt=t_end/(n_steps-1)

t=np.linspace(0,t_end,n_steps)
w=np.zeros((2,n_steps))
w[0,:int(t_switch/dt)+1]=m1_sol(t[:int(t_switch/dt)+1],w0[0],mu_m,sigma_m)
w[1,:int(t_switch/dt)+1]=c1_sol(t[:int(t_switch/dt)+1],w0[1],sigma_c)

t2=t_end-t_switch
sol=solve_ivp(f,(0,t2),w[:,int(t_switch/dt)],t_eval=np.linspace(0,t2,int(t2/dt)+1),args=(mu_m,mu_c,sigma_m,sigma_c,l))
w[:,int(t_switch/dt):]=sol.y

acc=np.zeros(n_steps)
acc[:int(t_switch/dt)]=accuracy(w[:,:int(t_switch/dt)],mu_m,0,sigma_m,sigma_c,sigma)
acc[int(t_switch/dt):]=accuracy(w[:,int(t_switch/dt):],mu_m,mu_c,sigma_m,sigma_c,sigma)

fig,ax = plt.subplots()
ax.plot(t,w[0],label='w_m')
ax.plot(t,w[1],label='w_c')
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

# #vary sigma_m
# mu_m,mu_c=1,1
# sigma_c,sigma=0.5,0.5
# l=0.1
# w0=np.ones(2)*1e-1
# t_switch=20
# t_end=80
# n_steps=1001
# dt=t_end/(n_steps-1)

# t=np.linspace(0,t_end,n_steps)
# t2=t_end-t_switch

# fig,ax = plt.subplots()
# sigma_m_list=[0.5,1,2,5,10,20]
# sns.set_palette('plasma',len(sigma_m_list))
# for sm in sigma_m_list:
#     w=np.zeros((2,n_steps))
#     w[0,:int(t_switch/dt)+1]=m1_sol(t[:int(t_switch/dt)+1],w0[0],mu_m,sm)
#     w[1,:int(t_switch/dt)+1]=c1_sol(t[:int(t_switch/dt)+1],w0[1],sigma_c)
#     sol=solve_ivp(f,(0,t2),w[:,int(t_switch/dt)],t_eval=np.linspace(0,t2,int(t2/dt)+1),args=(mu_m,mu_c,sm,sigma_c,l))
#     w[:,int(t_switch/dt):]=sol.y
#     acc=np.zeros(n_steps)
#     acc[:int(t_switch/dt)]=accuracy(w[:,:int(t_switch/dt)],mu_m,0,sm,sigma_c,sigma)
#     acc[int(t_switch/dt):]=accuracy(w[:,int(t_switch/dt):],mu_m,mu_c,sm,sigma_c,sigma)
#     plt.axvline(t_switch,linestyle='dashed',color='red')
#     ax.plot(t,acc,label='SNR_m='+str(1/sm))
# ax.set_xlabel('t')
# ax.set_ylabel('accuracy')
# ax.set_xticks([0,20,40,60,80])
# ax.set_xticklabels([0,'t_switch=20',40,60,80])
# ax.legend()
# plt.title('SNR_c=2')

# #vary lambda
# mu_m,mu_c=1,1
# sigma_m,sigma_c,sigma=5,0.5,0.5
# w0=np.ones(2)*1e-1
# t_switch=20
# t_end=80
# n_steps=1001
# dt=t_end/(n_steps-1)

# t=np.linspace(0,t_end,n_steps)
# t2=t_end-t_switch

# fig,ax = plt.subplots()
# lambda_list=[0.01,0.02,0.05,0.1,0.2,0.5]
# sns.set_palette('plasma',len(lambda_list))
# for l in lambda_list:
#     w=np.zeros((2,n_steps))
#     w[0,:int(t_switch/dt)+1]=m1_sol(t[:int(t_switch/dt)+1],w0[0],mu_m,sigma_m)
#     w[1,:int(t_switch/dt)+1]=c1_sol(t[:int(t_switch/dt)+1],w0[1],sigma_c)
#     sol=solve_ivp(f,(0,t2),w[:,int(t_switch/dt)],t_eval=np.linspace(0,t2,int(t2/dt)+1),args=(mu_m,mu_c,sigma_m,sigma_c,l))
#     w[:,int(t_switch/dt):]=sol.y
#     acc=np.zeros(n_steps)
#     acc[:int(t_switch/dt)]=accuracy(w[:,:int(t_switch/dt)],mu_m,0,sigma_m,sigma_c,sigma)
#     acc[int(t_switch/dt):]=accuracy(w[:,int(t_switch/dt):],mu_m,mu_c,sigma_m,sigma_c,sigma)
#     plt.axvline(t_switch,linestyle='dashed',color='red')
#     ax.plot(t,acc,label='lambda='+str(l))
# ax.set_xlabel('t')
# ax.set_ylabel('accuracy')
# ax.set_xticks([0,20,40,60,80])
# ax.set_xticklabels([0,'t_switch=20',40,60,80])
# ax.legend()

# #vary init
# mu_m,mu_c=1,1
# sigma_m,sigma_c,sigma=5,0.5,0.5
# l=0.1
# t_switch=20
# t_end=80
# n_steps=1001
# dt=t_end/(n_steps-1)

# t=np.linspace(0,t_end,n_steps)
# t2=t_end-t_switch

# fig,ax = plt.subplots()
# init_list=[1e-3,1e-2,1e-1,1]
# sns.set_palette('plasma',len(init_list))
# for init_scale in init_list:
#     w0=np.ones(2)*init_scale
#     w0[0]=1e-1
#     w=np.zeros((2,n_steps))
#     w[0,:int(t_switch/dt)+1]=m1_sol(t[:int(t_switch/dt)+1],w0[0],mu_m,sigma_m)
#     w[1,:int(t_switch/dt)+1]=c1_sol(t[:int(t_switch/dt)+1],w0[1],sigma_c)
#     sol=solve_ivp(f,(0,t2),w[:,int(t_switch/dt)],t_eval=np.linspace(0,t2,int(t2/dt)+1),args=(mu_m,mu_c,sigma_m,sigma_c,l))
#     w[:,int(t_switch/dt):]=sol.y
#     acc=np.zeros(n_steps)
#     acc[:int(t_switch/dt)]=accuracy(w[:,:int(t_switch/dt)],mu_m,0,sigma_m,sigma_c,sigma)
#     acc[int(t_switch/dt):]=accuracy(w[:,int(t_switch/dt):],mu_m,mu_c,sigma_m,sigma_c,sigma)
#     plt.axvline(t_switch,linestyle='dashed',color='red')
#     ax.plot(t,acc,label='wc_0='+str(init_scale))
# ax.set_xlabel('t')
# ax.set_ylabel('accuracy')
# ax.set_xticks([0,20,40,60,80])
# ax.set_xticklabels([0,'t_switch=20',40,60,80])
# ax.legend()
