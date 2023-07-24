# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 18:17:25 2023

@author: Rayan
"""

import numpy as np
import matplotlib.pyplot as plt
axis_font = {'fontname':'Arial', 'size':'17'}




notp_11_v0 = np.loadtxt('o4_11_no_tp03_U3V0.txt')
notp_11_v15 = np.loadtxt('o4_11_no_tp03_U3V15.txt')

notp_33_v0 = np.loadtxt('o4_33_no_tp03_U3V0.txt')
notp_33_v15 = np.loadtxt('o4_33_no_tp03_U3V15.txt')

q =range(len(notp_11_v0[:,0]))


plt.figure(figsize=[12,7])
plt.suptitle('Tperp = 0.5, Tperp_p=0.0',**axis_font)
plt.subplot(1,2,1)

plt.errorbar(q, notp_11_v0[:,0],yerr=notp_11_v0[:,1],label=r'U=3,V=0, $\Sigma_{00}$')
plt.errorbar(q, notp_11_v15[:,0],yerr=notp_11_v15[:,1],label='U=3,V=1.5,$\Sigma_{00}$')
plt.errorbar(q, notp_33_v0[:,0],yerr=notp_33_v0[:,1],label='U=3,V=0,$\Sigma_{11}$')
plt.errorbar(q, notp_33_v15[:,0],yerr=notp_33_v15[:,1],label='U=3,V=1.5,$\Sigma_{11}$')
plt.ylabel(r'Re $\Sigma$',**axis_font)
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend(fontsize=11)
plt.subplot(1,2,2)

plt.errorbar(q, notp_11_v0[:,2],yerr=notp_11_v0[:,3],label=r'U=3,V=0, $\Sigma_{00}$')
plt.errorbar(q, notp_11_v15[:,2],yerr=notp_11_v15[:,3],label='U=3,V=1.5,$\Sigma_{00}$')
plt.errorbar(q, notp_33_v0[:,2],yerr=notp_33_v0[:,3],label='U=3,V=0,$\Sigma_{11}$')
plt.errorbar(q, notp_33_v15[:,2],yerr=notp_33_v15[:,3],label='U=3,V=1.5,$\Sigma_{11}$')
plt.ylabel(r'Im $\Sigma$',**axis_font)
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])

plt.savefig('bilayer_tp05_tpp0.pdf', dpi=2400, bbox_inches='tight')
plt.show()


notp_11_v0 = np.loadtxt('o4_11_tp03_U3V0.txt')
notp_11_v15 = np.loadtxt('o4_11_tp03_U3V15.txt')

notp_33_v0 = np.loadtxt('o4_33_tp03_U3V0.txt')
notp_33_v15 = np.loadtxt('o4_33_tp03_U3V15.txt')

q =range(len(notp_11_v0[:,0]))

plt.figure(figsize=[12,7])
plt.suptitle('Tperp = 0.5, Tperp_p=0.3',**axis_font)
plt.subplot(1,2,1)

plt.errorbar(q, notp_11_v0[:,0],yerr=notp_11_v0[:,1],label=r'U=3,V=0, $\Sigma_{00}$')
plt.errorbar(q, notp_11_v15[:,0],yerr=notp_11_v15[:,1],label='U=3,V=1.5,$\Sigma_{00}$')
plt.errorbar(q, notp_33_v0[:,0],yerr=notp_33_v0[:,1],label='U=3,V=0,$\Sigma_{11}$')
plt.errorbar(q, notp_33_v15[:,0],yerr=notp_33_v15[:,1],label='U=3,V=1.5,$\Sigma_{1}$')
plt.ylabel(r'Re $\Sigma$',**axis_font)
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend(fontsize=11)
plt.subplot(1,2,2)

plt.errorbar(q, notp_11_v0[:,2],yerr=notp_11_v0[:,3],label=r'U=3,V=0, $\Sigma_{00}$')
plt.errorbar(q, notp_11_v15[:,2],yerr=notp_11_v15[:,3],label='U=3,V=1.5,$\Sigma_{00}$')
plt.errorbar(q, notp_33_v0[:,2],yerr=notp_33_v0[:,3],label='U=3,V=0,$\Sigma_{11}$')
plt.errorbar(q, notp_33_v15[:,2],yerr=notp_33_v15[:,3],label='U=3,V=1.5,$\Sigma_{11}$')
plt.ylabel(r'Im $\Sigma$',**axis_font)
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])

plt.savefig('bilayer_tp05_tpp03.pdf', dpi=2400, bbox_inches='tight')
plt.show()


