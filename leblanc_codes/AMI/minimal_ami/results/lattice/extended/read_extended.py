# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 15:34:58 2023

@author: Rayan
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
file_path = "real_freq/Rf_an_v8_o4_neg.txt"
# data = np.loadtxt(file_path)es\
U =4
V= -4
J = 0

def symmetrize(filtered_df):
    filtered_df1 = filtered_df[(filtered_df[5] == 0.000)]
    filtered_df2 = filtered_df[(filtered_df[4] == 3.141590) & (filtered_df[5] > 0.00)]
    filtered_df3 = filtered_df[(filtered_df[5] == filtered_df[4]) & (filtered_df[5] < 3.14)].iloc[::-1]
    
    symmetry_cut_df = pd.concat([filtered_df1, filtered_df2, filtered_df3])

    return symmetry_cut_df

def spectral(real,imag,omega,qx,qy, i =0):
    e = -2*(np.cos(qx)+np.cos(qy))
    imag =imag-i
    den = (omega - e - real)**2 + (imag)**2
    num = imag
    return -num/(den*np.pi)
Uindex= np.arange(0,4,1)
Uvalues = np.append(np.ones(2)*V, np.ones(2)*U)

print(Uvalues)


data_U = {
    'Type_of_U': Uindex,
    'U': Uvalues
}



df_U = pd.DataFrame(data_U)
df = pd.read_csv(file_path, delim_whitespace=True, skiprows=1, header=None)
# df = df.iloc[:, :-2]

for i in range(4):
    df = df.merge(df_U, left_on=14 + i, right_on='Type_of_U', how='left')
    df.fillna(1, inplace=True)
    print(df)
    df[10] = df[10] * df['U']
    df[11] = df[11] * df['U']
    df[12] = df[12] * df['U']
    df[13] = df[13] * df['U']
    df = df.drop(columns=['U', 'Type_of_U'])

print(df)

# Group by the first 10 columns and sum the values in columns 9 and 10
df =df[(df[0] < 4)]
grouped_df = df.groupby(df.columns[1:10].tolist())[[10, 11, 12, 13]].sum().reset_index()
# grouped_df= df.groupby(df.columns[:10].tolist())[[10, 11, 12, 13]].sum().reset_index()

# filtered_df = grouped_df[(grouped_df[8] == 1) & (grouped_df[9] == 1)]
filtered_df = grouped_df[ (grouped_df[8] == 1) & (grouped_df[9] == 1) ]
# filtered_df = grouped_df[(grouped_df[8] == 1) & (grouped_df[9] == 1) & (grouped_df[0] == 3)]

filtered_df = symmetrize(filtered_df)
print(filtered_df)




data= filtered_df.to_numpy()

print(data)

# result = np.column_stack((np.asarray(data[:,6]),np.asarray(data[:,11]), np.asarray(data[:,12])/1.5))
# # result = np.column_stack((np.asarray(data[:,10]+U/2), np.asarray(data[:,11]),np.asarray(data[:,12]), np.asarray(data[:,13])))
# np.savetxt('mfreq/freq_b5_vn3_udn.txt',result)

 
plt.subplot(1,2,1)
plt.errorbar(range(len(data[:,9])),data[:,9],yerr=data[:,10]/1.5,fmt='x-')
plt.subplot(1,2,2)
plt.errorbar(range(len(data[:,9])),data[:,11],yerr=data[:,12]/1.5,fmt='x-')
plt.show()
# plt.subplot(1,2,1)
# plt.errorbar(range(len(data[:,9])),data[:,10]+U/2,yerr=data[:,11],fmt='x-')
# plt.subplot(1,2,2)
# plt.errorbar(range(len(data[:,9])),data[:,12],yerr=data[:,13],fmt='x-')
# plt.show()


# def plotspec(data):
#     y = spectral(data[:,9],data[:,11],data[:,5],data[:,4],data[:,3],0.125)
#     print(y)
#     plt.subplot(1,2,1)
#     plt.errorbar(data[:,5],y,yerr=data[:,10]*0,fmt='x-')
#     plt.subplot(1,2,2)
#     plt.errorbar(data[:,5],data[:,11],yerr=data[:,12],fmt='x-')
#     plt.show()
#     plt.close()

# plotspec(data)
# grouped_df.to_csv('data.txt', sep=' ', header=False, index=False, float_format='%.7f')
# grouped_df_ord.to_csv('data_ord.txt', sep=' ', header=False, index=False, float_format='%.7f')


# print(grouped_df)