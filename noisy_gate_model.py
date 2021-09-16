# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 12:50:30 2021

@author: leoto
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import seaborn as sns

#w=(wm,gm,wc,gc)

def f_noise(w,mu_m,mu_c,sigma_m,sigma_c,l,sigma_ksi):
    fw=np.zeros(4)
    ksi=np.random.normal(0,sigma_ksi,4)
    fw[0]=-(mu_m**2+sigma_m**2)*w[1]**2*w[0]-mu_m*mu_c*w[1]*w[2]*w[3]+mu_m*w[1]
    fw[1]=-(mu_m**2+sigma_m**2)*w[0]**2*w[1]-mu_m*mu_c*w[0]*w[2]*w[3]+mu_m*w[0]-l*np.sign(w[1])
    fw[2]=-(mu_c**2+sigma_c**2)*w[3]**2*w[2]-mu_m*mu_c*w[0]*w[1]*w[3]+mu_c*w[3]
    fw[3]=-(mu_c**2+sigma_c**2)*w[2]**2*w[3]-mu_m*mu_c*w[0]*w[1]*w[2]+mu_c*w[2]-l*np.sign(w[3])
    return fw+ksi

def accuracy(w,mu_m,mu_c,sigma_m,sigma_c,sigma):
    sigma_tot2=sigma**2+(w[0]*w[1]*sigma_m)**2+(w[2]*w[3]*sigma_c)**2
    return (1+erf((w[0]*w[1]*mu_m+w[2]*w[3]*mu_c)/np.sqrt(2*sigma_tot2)))/2

def integrate(w0,mu_m,mu_c,sigma_m,sigma_c,l,sigma_ksi,dt,t_end):
    nsteps=int(t_end/dt)
    w=np.zeros((4,nsteps+1))
    w[:,0]=w0
    for i in range(nsteps):
        w[:,i+1]=w[:,i]+dt*f_noise(w[:,i],mu_m,mu_c,sigma_m,sigma_c,l,sigma_ksi)
    return w

# mu_m,mu_c=1,1
# sigma_m,sigma_c,sigma,sigma_ksi=5,0.5,0.5,0.
# l=0.1
# w0=np.array([0.05,0.5,0.05,0.5])
# t_switch=20
# t_end=100
# dt=0.05
# n_steps=int(t_end/dt)

# t=np.linspace(0,t_end,n_steps+1)
# w=np.zeros((4,n_steps+1))
# w[:,:int(t_switch/dt)+1]=integrate(w0,mu_m,0,sigma_m,sigma_c,l,sigma_ksi,dt,t_switch)

# t2=t_end-t_switch
# w[:,int(t_switch/dt):]=integrate(w[:,int(t_switch/dt)],mu_m,mu_c,sigma_m,sigma_c,l,sigma_ksi,dt,t2)

# acc=np.zeros(n_steps)
# acc[:int(t_switch/dt)]=accuracy(w[:,:int(t_switch/dt)],mu_m,0,sigma_m,sigma_c,sigma)
# acc[int(t_switch/dt):]=accuracy(w[:,int(t_switch/dt):],mu_m,mu_c,sigma_m,sigma_c,sigma)

# fig,ax = plt.subplots()
# ax.plot(t,w[0],label='w_m')
# ax.plot(t,w[1],label='g_m')
# ax.plot(t,w[2],label='w_c')
# ax.plot(t,w[3],label='g_c')
# plt.axvline(t_switch,linestyle='dashed',color='red')
# ax2=ax.twinx()
# ax2.plot(t,acc,color='black',label='accuracy')
# ax.set_xlabel('t')
# ax.set_ylabel('w')
# ax2.set_ylabel('accuracy')
# ax.set_xticks([0,20,40,60,80,100])
# ax.set_xticklabels([0,'t_switch=20',40,60,80,100])
# h1, l1 = ax.get_legend_handles_labels()
# h2, l2 = ax2.get_legend_handles_labels()
# ax.legend(h1+h2, l1+l2)

# #vary sigma_m
# mu_m,mu_c=1,1
# sigma_c,sigma,sigma_ksi=0.5,0.5,0
# l=0.1
# w0=np.array([0.05,0.5,0.05,0.5])
# t_switch=20
# t_end=400
# dt=0.02
# n_steps=int(t_end/dt+1)

# t=np.linspace(0,t_end,n_steps)
# t2=t_end-t_switch

# fig,ax = plt.subplots()
# sigma_m_list=[0.5,1,2,5,10,20]
# sns.set_palette('plasma',len(sigma_m_list))
# for sigma_m in sigma_m_list:
#     w=np.zeros((4,n_steps))
#     w[:,:int(t_switch/dt)+1]=integrate(w0,mu_m,0,sigma_m,sigma_c,l,sigma_ksi,dt,t_switch)
#     w[:,int(t_switch/dt):]=integrate(w[:,int(t_switch/dt)],mu_m,mu_c,sigma_m,sigma_c,l,sigma_ksi,dt,t2)
#     acc=np.zeros(n_steps)
#     acc[:int(t_switch/dt)]=accuracy(w[:,:int(t_switch/dt)],mu_m,0,sigma_m,sigma_c,sigma)
#     acc[int(t_switch/dt):]=accuracy(w[:,int(t_switch/dt):],mu_m,mu_c,sigma_m,sigma_c,sigma)
#     ax.plot(t,acc,label='SNR_m='+str(1/sigma_m))
# plt.axvline(t_switch,linestyle='dashed',color='red')
# ax.set_xlabel('t')
# ax.set_ylabel('accuracy')
# ax.set_xticks([0,20,100,200,300,400])
# ax.set_xticklabels([0,'             t_switch=20',100,200,300,400])
# ax.legend()
# plt.title('SNR_c=2')

# #vary lambda
# mu_m,mu_c=1,1
# sigma_m,sigma_c,sigma,sigma_ksi=5,0.5,0.5,0
# w0=np.array([0.05,0.5,0.05,0.5])
# t_switch=20
# t_end=200
# dt=0.02
# n_steps=int(t_end/dt+1)

# t=np.linspace(0,t_end,n_steps)
# t2=t_end-t_switch

# fig,ax = plt.subplots()
# lambda_list=[0.02,0.05,0.1,0.2]
# sns.set_palette('plasma',len(lambda_list))
# for l in lambda_list:
#     w=np.zeros((4,n_steps))
#     w[:,:int(t_switch/dt)+1]=integrate(w0,mu_m,0,sigma_m,sigma_c,l,sigma_ksi,dt,t_switch)
#     w[:,int(t_switch/dt):]=integrate(w[:,int(t_switch/dt)],mu_m,mu_c,sigma_m,sigma_c,l,sigma_ksi,dt,t2)
#     acc=np.zeros(n_steps)
#     acc[:int(t_switch/dt)]=accuracy(w[:,:int(t_switch/dt)],mu_m,0,sigma_m,sigma_c,sigma)
#     acc[int(t_switch/dt):]=accuracy(w[:,int(t_switch/dt):],mu_m,mu_c,sigma_m,sigma_c,sigma)
#     ax.plot(t,acc,label='lambda='+str(l))
# plt.axvline(t_switch,linestyle='dashed',color='red')
# ax.set_xlabel('t')
# ax.set_ylabel('accuracy')
# ax.set_xticks([0,20,50,100,150,200])
# ax.set_xticklabels([0,'   t_switch=20',50,100,150,200])
# ax.legend()

# #vary init
# mu_m,mu_c=1,1
# sigma_m,sigma_c,sigma,sigma_ksi=5,0.5,0.5,0
# l=0.1
# sigma0_w,sigma0_g=0.05,0.5
# t_switch=20
# t_end=200
# dt=0.02
# n_steps=int(t_end/dt+1)

# t=np.linspace(0,t_end,n_steps)
# t2=t_end-t_switch

# fig,ax = plt.subplots()
# sns.set_palette('tab10')
# for k in range(10):
#     w=np.zeros((4,n_steps))
#     w[0,0],w[1,0],w[2,0],w[3,0]=np.random.normal(0,sigma0_w),np.random.normal(0,sigma0_g),np.random.normal(0,sigma0_w),np.random.normal(0,sigma0_g)
#     print(w[:,0])
#     w[:,:int(t_switch/dt)+1]=integrate(w[:,0],mu_m,0,sigma_m,sigma_c,l,sigma_ksi,dt,t_switch)
#     w[:,int(t_switch/dt):]=integrate(w[:,int(t_switch/dt)],mu_m,mu_c,sigma_m,sigma_c,l,sigma_ksi,dt,t2)
#     acc=np.zeros(n_steps)
#     acc[:int(t_switch/dt)]=accuracy(w[:,:int(t_switch/dt)],mu_m,0,sigma_m,sigma_c,sigma)
#     acc[int(t_switch/dt):]=accuracy(w[:,int(t_switch/dt):],mu_m,mu_c,sigma_m,sigma_c,sigma)
#     ax.plot(t,acc)
# plt.axvline(t_switch,linestyle='dashed',color='red')
# ax.set_xlabel('t')
# ax.set_ylabel('accuracy')
# ax.set_xticks([0,20,50,100,150,200])
# ax.set_xticklabels([0,'   t_switch=20',50,100,150,200])   

#with noise
mu_m,mu_c=1,1
sigma_m,sigma_c,sigma,sigma_ksi=5,0.5,0.5,0.05
l=0.1
w0=np.array([0.05,0.5,0.05,0.5])
t_switch=20
t_end=200
dt=0.02
n_steps=int(t_end/dt+1)

t=np.linspace(0,t_end,n_steps)
t2=t_end-t_switch

fig,ax = plt.subplots()
sns.set_palette('tab10')
for k in range(10):
    w=np.zeros((4,n_steps))
    w[:,:int(t_switch/dt)+1]=integrate(w0,mu_m,0,sigma_m,sigma_c,l,sigma_ksi,dt,t_switch)
    w[:,int(t_switch/dt):]=integrate(w[:,int(t_switch/dt)],mu_m,mu_c,sigma_m,sigma_c,l,sigma_ksi,dt,t2)
    acc=np.zeros(n_steps)
    acc[:int(t_switch/dt)]=accuracy(w[:,:int(t_switch/dt)],mu_m,0,sigma_m,sigma_c,sigma)
    acc[int(t_switch/dt):]=accuracy(w[:,int(t_switch/dt):],mu_m,mu_c,sigma_m,sigma_c,sigma)
    ax.plot(t,acc)
plt.axvline(t_switch,linestyle='dashed',color='red')
ax.set_xlabel('t')
ax.set_ylabel('accuracy')
ax.set_xticks([0,20,50,100,150,200])
ax.set_xticklabels([0,'   t_switch=20',50,100,150,200])   