# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 14:50:59 2023

@author: Rayan
"""

import numpy as np
import matplotlib.pyplot as plt
axis_font = {'fontname':'Arial', 'size':'21', 'weight':'bold'}

o12_ami = np.loadtxt('james_o12.txt')
o12_mb = np.loadtxt('mband_o12.txt')
plt.title('beta=4.0,mu=0.01',**axis_font)
plt.errorbar(np.arange(13), o12_ami[:,2],yerr=o12_ami[:,3],label='AMI-o12')
plt.errorbar(np.arange(13), o12_mb[:,2],yerr=o12_mb[:,3],label='mband-o12')
plt.ylabel(r'$P_{d}$',**axis_font)
plt.xticks([0, 4, 8, 12], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.legend(fontsize=11)
plt.show()