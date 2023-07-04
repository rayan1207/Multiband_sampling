# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 17:27:29 2023

@author: Rayan
"""

import numpy as np
import matplotlib.pyplot as plt
axis_font = {'fontname':'Arial', 'size':'15'}
data = [1,2,4,8,15,30,60]
node1 = [6500,13000,26000,52000,92500,188500,373000]

node4 = [6000,11500,24000,46500,88000,165000,345000]

node8 = [5000,9000,17500,34000,63500,127500,257000]

data = np.asarray(data)
node1 = np.asarray(node1)
node4= np.asarray(node4)
node8 =np.asarray(node8)
plt.subplot(1,2,1)
plt.plot(data,node1,'x-', label='core = 1')
plt.plot(data,node4,'^-',label='core = 4')
plt.plot(data,node8,'h-',label='core = 8')
plt.xlabel('Cpu-minutes',**axis_font)
plt.ylabel('sample/core',**axis_font)
plt.legend(fontsize=9)
plt.subplot(1,2,2)
plt.plot(data,node1/data,'x-', label='core = 1')
plt.plot(data,node4/data,'^-',label='core = 4')
plt.plot(data,node8/data,'h-',label='core = 8')
plt.xlabel('Cpu-minutes',**axis_font)
plt.ylabel('sample / (core*time)',**axis_font)
plt.tight_layout()
plt.legend(fontsize=9)
plt.show()
