# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 14:24:58 2023

@author: Rayan
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
axis_font = {'fontname':'Arial', 'size':'16'}
axis_font1 = {'fontname':'Arial', 'size':'13'}

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
plt.legend(fontsize=13)
plt.subplot(1,2,2)
plt.errorbar(mfreq0_b5_v2[:,0],freq0_b5_v2[:,1],yerr=freq0_b5_v2[:,2]*1.5,fmt='X-',label=r'$U=2,V=0$')
plt.errorbar(mfreq_b5_v2[:,0],freq_b5_v2[:,1],yerr=freq_b5_v2[:,2]*1.5,fmt='X-',label=r'$U=2,V_{ud}=V_{uu} =0.25$')
plt.errorbar(mfreq_b5_v2n[:,0],freq_b5_v2n[:,1],yerr=freq_b5_v2n[:,2]*1.5,fmt='X-',label=r'$U=2,V_{ud}=V_{uu} =-0.25$')
plt.errorbar(mfreq_b5_v2n_udn[:,0],freq_b5_v2n_udn[:,1],yerr=freq_b5_v2n_udn[:,2]*1.5,fmt='X-',label=r'$U=2,V_{ud}=-0.25,V_{uu} =0.25$')
plt.xlabel('Matsubara_freq',**axis_font)
plt.ylabel(r'$ Im \Sigma (q = [\pi/2,\pi/2])$',**axis_font)
plt.legend(fontsize=13)
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
plt.legend(fontsize=12)
plt.subplot(1,2,2)
plt.ylabel(r'$ Im \Sigma (q = [\pi/2,\pi/2])$',**axis_font)
plt.errorbar(mfreq0_b5_v3[:,0],freq0_b5_v3[:,1],yerr=freq0_b5_v3[:,2]*1.5,fmt='X-',label=r'$U=3.5,V=0$')
plt.errorbar(mfreq_b5_v2[:,0],freq_b5_v3[:,1],yerr=freq_b5_v3[:,2]*1.5,fmt='X-',label=r'$U=3.5,V_{ud}=V_{uu} =0.438$')
plt.errorbar(mfreq_b5_v2n[:,0],freq_b5_v3n[:,1],yerr=freq_b5_v3n[:,2]*1.5,fmt='X-',label=r'$U=3.5,V_{ud}=V_{uu} =-0.438$')
plt.errorbar(mfreq_b5_v2n_udn[:,0],freq_b5_v3n_udn[:,1],yerr=freq_b5_v3n_udn[:,2]*1.5,fmt='X-',label=r'$U=3.5,V_{ud}=-0.438,V_{uu} =0.438$')
plt.legend(fontsize=12)
plt.xlabel('Matsubara_freq',**axis_font)
plt.savefig('b5_mfreq_U34.pdf', dpi=2400, bbox_inches='tight')
plt.show()


mfreq0 = np.loadtxt('mfreq0_b5_u2v2.txt')
mfreq1 = np.loadtxt('mfreq1_b5_u2v2.txt')

Vuu = np.reshape(mfreq0[:,0], (7,7))*2
Vup = np.reshape(mfreq0[:,1], (7,7))*2


imag0 = np.reshape(mfreq0[:,2], (7,7))
imag1 = np.reshape(mfreq1[:,2], (7,7))
err = np.reshape(mfreq0[:,3], (7,7))
delta = (imag1-imag0)


levels = np.linspace(-delta.max()/1.5, delta.max()/1.5, 50)  # Adjust the number of levels as needed


colors = [(0.0, 'blue'),(0.3, 'blue'), (0.5, 'white'), (1, 'red')]
cmap = LinearSegmentedColormap.from_list('custom_bwr', colors)
contour1 = plt.pcolormesh(Vuu, Vup, delta, cmap=cmap, alpha=1,vmin=-delta.max()/1.1,vmax=delta.max()/1.1)
plt.ylabel(r'$V_{\uparrow \downarrow}$',**axis_font)
plt.xlabel(r'$V_{\uparrow \uparrow}$',**axis_font)
plt.text(-0.2,0.4,r'U/t =2 $ q=[\pi,0]$',**axis_font1)

cut1 = np.diag(delta)
cutv1=np.diag(Vuu)
cut2 = np.diag(np.fliplr(delta))
cutv2=np.diag(Vup)


cut1_err = np.diag(err)
cut2_err = np.diag(np.fliplr(err))
cbar = plt.colorbar(contour1)
cbar.set_label(r'$Im  \Sigma(i\omega_{1})-\Sigma(i\omega_{0}) $', rotation=270, labelpad=15,**axis_font1)  # LaTeX label


# Show the plot
plt.show()

print(np.diag(Vuu),np.diag(Vup))

# plt.errorbar(np.diag(Vuu),np.diag(imag0),yerr=label='mfreq0')
plt.errorbar(cutv1, cut1,yerr=cut1_err,label=r'$V_{\uparrow \downarrow} = V_{\uparrow \uparrow}$')
plt.errorbar(cutv2, cut2,yerr=cut2_err,label=r'$V_{\uparrow \downarrow} = -V_{\uparrow \uparrow}$')
plt.xlabel(r'$V_{\uparrow \uparrow}$',**axis_font)
plt.legend(fontsize=12)
plt.ylabel(r'$Im   \Sigma(i\omega_{1})-\Sigma(i\omega_{0}) $',**axis_font)
plt.legend()
plt.show()





