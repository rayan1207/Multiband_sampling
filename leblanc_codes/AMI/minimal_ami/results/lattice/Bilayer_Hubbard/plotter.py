# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 02:00:39 2023

@author: Rayan
"""

import numpy as np
import matplotlib.pyplot as plt

data =np.loadtxt('filtered_output.txt')

print(data[:,9])

plt.plot(data[:,9],'x')
plt.show()