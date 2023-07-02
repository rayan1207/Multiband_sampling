# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 03:43:24 2023

@author: Rayan
"""

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import os
axis_font = {'fontname':'Arial', 'size':'15'}
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)


script_folder = os.path.dirname(os.path.abspath(__file__))
folder_path = script_folder  # Use the script folder as the folder path

data_list1=[]
for filename in os.listdir(folder_path):
    if filename.endswith('.txt'):
        file_path = os.path.join(folder_path, filename)
        file_data = np.loadtxt(file_path)
        data_list1.append(file_data)



num =100
def sum_rows(matrix):
    row_sums = np.sum(matrix, axis=1)
    return row_sums

def iterate(data_list):
    real = np.zeros(num)
    imag  = np.zeros(num)
    for i in data_list:
       real=real+ sum_rows(np.reshape(i[:,2],(num,int(len(i[:,1])/num))))
       imag =imag+ sum_rows(np.reshape(i[:,3],(num,int(len(i[:,1])/num))))
    return [real,imag]


plot= iterate(data_list1)

freq=[]
for i in range(num):
    freq.append((2*i+1)*np.pi/50)
    

#np.savetxt('sigma11_o4.dat',plot)


plt.subplot(1,2,1)
plt.plot(freq,plot[0],'x-')
plt.subplot(1,2,2)
plt.plot(freq,plot[1],'.-')
plt.show()


print(plot[0])
