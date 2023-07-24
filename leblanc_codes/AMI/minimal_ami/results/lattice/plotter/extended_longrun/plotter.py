# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 00:18:59 2023

@author: Rayan
"""

import numpy as np
import matplotlib.pyplot as plt
axis_font = {'fontname':'Arial', 'size':'21'}


db_v0_u2_im = np.loadtxt('db_v0_u2_im.txt')
db_v4_u2_im = np.loadtxt('db_v4_u2_im.txt')
db_v8_u2_im = np.loadtxt('db_v8_u2_im.txt')

db_v0_u2_re = np.loadtxt('db_v0_u2_re1.txt')
db_v4_u2_re = np.loadtxt('db_v4_u2_re.txt')
db_v8_u2_re = np.loadtxt('db_v8_u2_re.txt')


mband_v0_u2 = np.loadtxt('mband_v0_u2.txt')
mband_v8_u2 = np.loadtxt('mband_v8_u2.txt')
mband_v4_u2 = np.loadtxt('mband_v4_u2.txt')
 

james_v0_u2 = np.loadtxt('james_v0_u2.txt')


x = np.arange(0,25,2)
x1 = np.arange(0,25-24/36,24/36)




print(len(x1))

plt.figure(figsize=[22,15])
plt.subplot(2,3,4)
plt.ylabel(r'Im $\Sigma$',**axis_font)
plt.plot(db_v0_u2_im[:,0],db_v0_u2_im[:,1],'o-',label='Dual-Boson')
plt.errorbar(x,mband_v0_u2[:,2],yerr=mband_v0_u2[:,3],fmt='x-',label='mband-AMI-o4')
plt.errorbar(x1,james_v0_u2[:,4],yerr=james_v0_u2[:,5],fmt='.-',label='James-AMI-o4')
plt.legend(fontsize='14')
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.text(13.5,-0.10,'U=2.0\nV=0.0',fontsize=15)

plt.subplot(2,3,5)
plt.plot(db_v8_u2_im[:,0],db_v8_u2_im[:,1],'o-',label='Dual-Boson')
plt.errorbar(x,mband_v8_u2[:,2],yerr=mband_v8_u2[:,3],fmt='x-',label='mband-AMI-o4')
plt.legend(fontsize='14')
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.text(13.5,-0.10,'U=2.0\nV=0.25',fontsize=15)

plt.subplot(2,3,6)
plt.plot(db_v4_u2_im[:,0],db_v4_u2_im[:,1],'o-',label='Dual-Boson')
plt.errorbar(x,mband_v4_u2[:,2],yerr=mband_v4_u2[:,3],fmt='x-',label='mband-AMI-o4')
plt.text(13.5,-0.14,'U=2.0\nV=0.50',fontsize=15)
plt.legend(fontsize='14')
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])



plt.subplot(2,3,1)
plt.ylabel(r'Re $\Sigma$',**axis_font)
plt.plot(db_v0_u2_re[:,0],db_v0_u2_re[:,1],'o',label='Dual-Boson')
plt.errorbar(x,mband_v0_u2[:,0],yerr=mband_v0_u2[:,1],fmt='x-',label='mband-AMI-o4')
plt.errorbar(x1,james_v0_u2[:,2]+1,yerr=james_v0_u2[:,3],fmt='.-',label='James-AMI-o4')
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend(fontsize='14')
plt.text(13,1.00,'U=2.0\nV=0.00',fontsize=15)

plt.subplot(2,3,2)
plt.plot(db_v8_u2_re[:,0],db_v8_u2_re[:,1],'o-',label='Dual-Boson')
plt.errorbar(x,mband_v8_u2[:,0],yerr=mband_v8_u2[:,1],fmt='x-',label='mband-AMI-o4')
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend(fontsize='14')
plt.text(13,0.95,'U=2.0\nV=0.25',fontsize=15)

plt.subplot(2,3,3)
plt.plot(db_v4_u2_re[:,0],db_v4_u2_re[:,1],'o-',label='Dual-Boson')
plt.errorbar(x,mband_v4_u2[:,0],yerr=mband_v4_u2[:,1],fmt='x-',label='mband-AMI-o4')
plt.text(13,0.95,'U=2.0\nV=0.50',fontsize=15)
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend(fontsize='14')
plt.savefig('DBvsMband.pdf', dpi=2400, bbox_inches='tight')
plt.show()



plt.figure(figsize=[8,5])
plt.subplot(1,2,1)
plt.ylabel(r'Re $\Sigma$',**axis_font)
plt.errorbar(x,mband_v0_u2[:,0],yerr=mband_v0_u2[:,1],fmt='x-',label='U=2,V=0')
plt.errorbar(x,mband_v8_u2[:,0],yerr=mband_v8_u2[:,1],fmt='x-',label='U=2,V=0.25')
plt.errorbar(x,mband_v4_u2[:,0],yerr=mband_v4_u2[:,1],fmt='x-',label='U=2,V=0.50')
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend(fontsize='10')
plt.subplot(1,2,2)
plt.ylabel(r'Im $\Sigma$',**axis_font)
plt.errorbar(x,mband_v0_u2[:,2],yerr=mband_v0_u2[:,3],fmt='x-',label='U=2,V=0')
plt.errorbar(x,mband_v8_u2[:,2],yerr=mband_v8_u2[:,3],fmt='x-',label='U=2,V=0.25')
plt.errorbar(x,mband_v4_u2[:,2],yerr=mband_v4_u2[:,3],fmt='x-',label='U=2,V=0.50')
plt.xticks([0, 8, 16, 24], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend(fontsize='10')
plt.tight_layout()
plt.savefig('MbandU2.pdf', dpi=2400, bbox_inches='tight')
plt.show()