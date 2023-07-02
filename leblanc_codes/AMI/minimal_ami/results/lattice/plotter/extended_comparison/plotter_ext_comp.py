# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 17:45:07 2023

@author: Rayan
"""

import numpy as np
import matplotlib.pyplot as plt
axis_font = {'fontname':'Arial', 'size':'21', 'weight':'bold'}

o2_U2_v0 = np.loadtxt('result_o2_U2_v0.txt')
o2_U4_v0 = np.loadtxt('result_o2_U4_v0.txt')

o2_U2_v05 = np.loadtxt('result_o2_U2_v05.txt')
o2_U4_v10 = np.loadtxt('result_o2_U4_v1.txt')


o2_U2_v025 = np.loadtxt('result_o2_U2_v025.txt')
o2_U4_v05 = np.loadtxt('result_o2_U4_v05.txt')

x = range(len(o2_U2_v0[:,0]))


plt.figure(figsize=[10,8])
plt.suptitle('truncated at order 2, beta=4')
plt.subplot(2,2,1)
plt.errorbar(x,o2_U2_v0[:,0],yerr=o2_U2_v0[:,1],fmt='x-',label='U=2,V=0')
plt.errorbar(x,o2_U2_v025[:,0],yerr=o2_U2_v025[:,1],fmt='.-',label='U=2,V=0.25')
plt.errorbar(x,o2_U2_v05[:,0],yerr=o2_U2_v05[:,1],fmt='o-',label='U=2,V=0.50')
plt.ylabel(r'Re $\Sigma$',**axis_font)
plt.xticks([0, 12, 24, 36], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend()
plt.subplot(2,2,2)
plt.errorbar(x,o2_U4_v0[:,0],yerr=o2_U4_v0[:,1],fmt='x-',label='U=4,V=0')
plt.errorbar(x,o2_U4_v05[:,0],yerr=o2_U4_v05[:,1],fmt='.-',label='U=4,V=0.50')
plt.errorbar(x,o2_U4_v10[:,0],yerr=o2_U4_v10[:,1],fmt='o-',label='U=4,V=1.0')
plt.xticks([0, 12, 24, 36], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend()
plt.subplot(2,2,3)
plt.errorbar(x,o2_U2_v0[:,2],yerr=o2_U2_v0[:,3],fmt='x-',label='U=2,V=0')
plt.errorbar(x,o2_U2_v025[:,2],yerr=o2_U2_v025[:,3],fmt='.-',label='U=2,V=0.25')
plt.errorbar(x,o2_U2_v05[:,2],yerr=o2_U2_v05[:,3],fmt='o-',label='U=2,V=0.50')
plt.xticks([0, 12, 24, 36], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.ylabel(r'Imag $\Sigma$',**axis_font)
plt.subplot(2,2,4)
plt.errorbar(x,o2_U4_v0[:,2],yerr=o2_U4_v0[:,3],fmt='x-',label='U=4,V=0')
plt.errorbar(x,o2_U4_v05[:,2],yerr=o2_U4_v05[:,3],fmt='.-',label='U=4,V=0.50')
plt.errorbar(x,o2_U4_v10[:,2],yerr=o2_U4_v10[:,3],fmt='o-',label='U=4,V=1.0')
plt.xticks([0, 12, 24, 36], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.savefig('extended_self_b4_o2.pdf', dpi=2400, bbox_inches='tight')
plt.show()


o4 = np.loadtxt('sigma_U4_o4_v0.dat')





o3_U2_v0 = np.loadtxt('result_o3_U2_v0.txt')
o3_U4_v0 = np.loadtxt('result_o3_U4_v0.txt')

o3_U2_v05 = np.loadtxt('result_o3_U2_v05.txt')
o3_U4_v10 = np.loadtxt('result_o3_U4_v1.txt')


o3_U2_v025 = np.loadtxt('result_o3_U2_v025.txt')
o3_U4_v05 = np.loadtxt('result_o3_U4_v05.txt')



plt.figure(figsize=[10,8])
plt.suptitle('truncated at order 3, beta=4')
plt.subplot(2,2,1)
plt.errorbar(x,o3_U2_v0[:,0],yerr=o3_U2_v0[:,1],fmt='x-',label='U=2,V=0')
plt.errorbar(x,o3_U2_v025[:,0],yerr=o3_U2_v025[:,1],fmt='.-',label='U=2,V=0.25')
plt.errorbar(x,o3_U2_v05[:,0],yerr=o3_U2_v05[:,1],fmt='o-',label='U=2,V=0.50')
plt.ylabel(r'Re $\Sigma$',**axis_font)
plt.xticks([0, 12, 24, 36], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend()
plt.subplot(2,2,2)
plt.errorbar(x,o3_U4_v0[:,0],yerr=o3_U4_v0[:,1]+o4[:,0],fmt='x-',label='U=4,V=0')
plt.errorbar(x,o3_U4_v05[:,0],yerr=o3_U4_v05[:,1],fmt='.-',label='U=4,V=0.50')
plt.errorbar(x,o3_U4_v10[:,0],yerr=o3_U4_v10[:,1],fmt='o-',label='U=4,V=1.0')
plt.xticks([0, 12, 24, 36], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend()
plt.subplot(2,2,3)
plt.errorbar(x,o3_U2_v0[:,2],yerr=o3_U2_v0[:,3],fmt='x-',label='U=2,V=0')
plt.errorbar(x,o3_U2_v025[:,2],yerr=o3_U2_v025[:,3],fmt='.-',label='U=2,V=0.25')
plt.errorbar(x,o3_U2_v05[:,2],yerr=o3_U2_v05[:,3],fmt='o-',label='U=2,V=0.50')
plt.xticks([0, 12, 24, 36], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.ylabel(r'Imag $\Sigma$',**axis_font)
plt.subplot(2,2,4)
plt.errorbar(x,o3_U4_v0[:,2],yerr=o3_U4_v0[:,3],fmt='x-',label='U=4,V=0')
plt.errorbar(x,o3_U4_v05[:,2],yerr=o3_U4_v05[:,3],fmt='.-',label='U=4,V=0.50')
plt.errorbar(x,o3_U4_v10[:,2],yerr=o3_U4_v10[:,3],fmt='o-',label='U=4,V=1.0')
plt.xticks([0, 12, 24, 36], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.savefig('extended_self_b4_o3.pdf', dpi=2400, bbox_inches='tight')
plt.show()


o4 = np.loadtxt('sigma_U4_o4_v0.dat')

plt.subplot(1,2,1)
plt.plot(o4[:,0]+2)
plt.subplot(1,2,2)
plt.plot(o4[:,2])
plt.show()
