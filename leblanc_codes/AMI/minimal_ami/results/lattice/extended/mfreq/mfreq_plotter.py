# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 14:24:58 2023

@author: Rayan
"""


import numpy as np
import matplotlib.pyplot as plt
axis_font = {'fontname':'Arial', 'size':'21'}


mfreq0_b5_v2 = np.loadtxt('mfreq0_b5_v2.txt')
mfreq_b5_v2 = np.loadtxt('mfreq_b5_v2.txt')
mfreq_b5_v2n = np.loadtxt('mfreq_b5_vn2.txt')
mfreq_b5_v2n_udn = np.loadtxt('mfreq_b5_vn2_udn.txt')

freq0_b5_v2 = np.loadtxt('freq0_b5_v2.txt')
freq_b5_v2 = np.loadtxt('freq_b5_v2.txt')
freq_b5_v2n = np.loadtxt('freq_b5_vn2.txt')
freq_b5_v2n_udn = np.loadtxt('freq_b5_vn2_udn.txt')

mfreq0_b5_v3 = np.loadtxt('mfreq0_b5_v3.txt')
mfreq_b5_v3 = np.loadtxt('mfreq_b5_v3.txt')
mfreq_b5_v3n = np.loadtxt('mfreq_b5_vn3.txt')
mfreq_b5_v3n_udn = np.loadtxt('mfreq_b5_vn3_udn.txt')


freq0_b5_v3 = np.loadtxt('freq0_b5_v3.txt')
freq_b5_v3 = np.loadtxt('freq_b5_v3.txt')
freq_b5_v3n = np.loadtxt('freq_b5_vn3.txt')
freq_b5_v3n_udn = np.loadtxt('freq_b5_vn3_udn.txt')

plt.figure(figsize=[13,8])
plt.subplot(1,2,1)
plt.errorbar(mfreq0_b5_v2[:,0],mfreq0_b5_v2[:,1],yerr=mfreq0_b5_v2[:,2]*1.5,fmt='X-',label=r'$U=2,V=0$')
plt.errorbar(mfreq_b5_v2[:,0],mfreq_b5_v2[:,1],yerr=mfreq_b5_v2[:,2]*1.5,fmt='X-',label=r'$U=2,V_{ud}=V_{uu} =0.25$')
plt.errorbar(mfreq_b5_v2n[:,0],mfreq_b5_v2n[:,1],yerr=mfreq_b5_v2n[:,2]*1.5,fmt='X-',label=r'$U=2,V_{ud}=V_{uu} =-0.25$')
plt.errorbar(mfreq_b5_v2n_udn[:,0],mfreq_b5_v2n_udn[:,1],yerr=mfreq_b5_v2n_udn[:,2]*1.5,fmt='X-',label=r'$U=2,V_{ud}=-0.25,V_{uu} =0.25$')
plt.xlabel('Matsubara_freq',**axis_font)
plt.ylabel(r'$ Im \Sigma (q = [\pi,0])$',**axis_font)
plt.legend()
plt.subplot(1,2,2)
plt.errorbar(mfreq0_b5_v2[:,0],freq0_b5_v2[:,1],yerr=freq0_b5_v2[:,2]*1.5,fmt='X-',label=r'$U=2,V=0$')
plt.errorbar(mfreq_b5_v2[:,0],freq_b5_v2[:,1],yerr=freq_b5_v2[:,2]*1.5,fmt='X-',label=r'$U=2,V_{ud}=V_{uu} =0.25$')
plt.errorbar(mfreq_b5_v2n[:,0],freq_b5_v2n[:,1],yerr=freq_b5_v2n[:,2]*1.5,fmt='X-',label=r'$U=2,V_{ud}=V_{uu} =-0.25$')
plt.errorbar(mfreq_b5_v2n_udn[:,0],freq_b5_v2n_udn[:,1],yerr=freq_b5_v2n_udn[:,2]*1.5,fmt='X-',label=r'$U=2,V_{ud}=-0.25,V_{uu} =0.25$')
plt.xlabel('Matsubara_freq',**axis_font)
plt.ylabel(r'$ Im \Sigma (q = [\pi/2,\pi/2])$',**axis_font)
plt.legend()
plt.savefig('b5_mfreq_U2.pdf', dpi=2400, bbox_inches='tight')
plt.show()

plt.figure(figsize=[12,8])

plt.subplot(1,2,1)
plt.xlabel('Matsubara_freq',**axis_font)
plt.ylabel(r'$ Im \Sigma (q = [\pi,0])$',**axis_font)
plt.errorbar(mfreq0_b5_v2[:,0],mfreq0_b5_v3[:,1],yerr=mfreq0_b5_v3[:,2]*1.5,fmt='X-',label=r'$U=3.5,V=0$')
plt.errorbar(mfreq_b5_v2[:,0],mfreq_b5_v3[:,1],yerr=mfreq_b5_v3[:,2]*1.5,fmt='X-',label=r'$U=3.5,V_{ud}=V_{uu} =0.438$')
plt.errorbar(mfreq_b5_v2n[:,0],mfreq_b5_v3n[:,1],yerr=mfreq_b5_v3n[:,2]*1.5,fmt='X-',label=r'$U=3.5,V_{ud}=V_{uu} =-0.438$')
plt.errorbar(mfreq_b5_v2n_udn[:,0],mfreq_b5_v3n_udn[:,1],yerr=mfreq_b5_v3n_udn[:,2]*1.5,fmt='X-',label=r'$U=3.5,V_{ud}=-0.438,V_{uu} =0.438$')
plt.legend()
plt.subplot(1,2,2)
plt.ylabel(r'$ Im \Sigma (q = [\pi/2,\pi/2])$',**axis_font)
plt.errorbar(mfreq0_b5_v3[:,0],freq0_b5_v3[:,1],yerr=freq0_b5_v3[:,2]*1.5,fmt='X-',label=r'$U=3.5,V=0$')
plt.errorbar(mfreq_b5_v2[:,0],freq_b5_v3[:,1],yerr=freq_b5_v3[:,2]*1.5,fmt='X-',label=r'$U=3.5,V_{ud}=V_{uu} =0.438$')
plt.errorbar(mfreq_b5_v2n[:,0],freq_b5_v3n[:,1],yerr=freq_b5_v3n[:,2]*1.5,fmt='X-',label=r'$U=3.5,V_{ud}=V_{uu} =-0.438$')
plt.errorbar(mfreq_b5_v2n_udn[:,0],freq_b5_v3n_udn[:,1],yerr=freq_b5_v3n_udn[:,2]*1.5,fmt='X-',label=r'$U=3.5,V_{ud}=-0.438,V_{uu} =0.438$')
plt.legend()
plt.xlabel('Matsubara_freq',**axis_font)
plt.savefig('b5_mfreq_U34.pdf', dpi=2400, bbox_inches='tight')
plt.show()



