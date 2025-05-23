import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file_path = "output1.txt"
# data = np.loadtxt(file_path)

U =-2
data_U = {
    'Type_of_U': [0, 1],
    'U': [U, U]
}

df_U = pd.DataFrame(data_U)
df = pd.read_csv(file_path, delim_whitespace=True, skiprows=1, header=None)

for i in range(4):
    df = df.merge(df_U, left_on=15 + i, right_on='Type_of_U', how='left')
    df.fillna(1, inplace=True)
    # print(df)
    df[11] = df[11] * df['U']
    df[12] = df[12] * df['U']
    df[13] = df[13] * df['U']
    df[14] = df[14] * df['U']
    df = df.drop(columns=['U', 'Type_of_U'])

#print(df)

# Group by the first 10 columns and sum the values in columns 9 and 10

# df = df[(df[0] >1)]
grouped_df = df.groupby(df.columns[1:11].tolist())[[ 11, 12, 13,14]].sum().reset_index()
# grouped_df_ord = df.groupby(df.columns[:10].tolist())[[10, 11, 12, 13]].sum().reset_index()

filtered_df = grouped_df[(grouped_df[9] ==12) & (grouped_df[10] == 12)  ] 

# def symmetrize(filtered_df):
#     filtered_df1 = filtered_df[(filtered_df[5] == 0.000)]
#     filtered_df2 = filtered_df[(filtered_df[4] == 3.141590) & (filtered_df[5] > 0.00)]
#     filtered_df3 = filtered_df[(filtered_df[5] == filtered_df[4]) & (filtered_df[5] < 3.14)].iloc[::-1]
    
#     symmetry_cut_df = pd.concat([filtered_df1, filtered_df2, filtered_df3])

#     return symmetry_cut_df

# def spectral(real,imag,omega,qx,qy, i =0):
#     e = -2*(np.cos(qx)+np.cos(qy))
#     den = (omega - e - real)**2 + (+imag)**2
#     num = imag
#     return -num/(den*np.pi)



# spectral(data[:,9] )

# filtered_df = symmetrize(filtered_df)
print(filtered_df)
data= filtered_df.to_numpy()
# filtered_df.to_csv('filtered_output.txt', sep=' ', header=False, index=False, float_format='%.7f')
# result = np.column_stack((np.asarray(data[:,3]),np.asarray(data[:,4]), np.asarray(data[:,9]),np.asarray(data[:,10])))
# # # result = np.column_stack((np.asarray(data[:,10]+U/2), np.asarray(data[:,11]),np.asarray(data[:,12]), np.asarray(data[:,13])))
# np.savetxt('pairing/mband_o12_pipi_mucut.txt',result)

# data= filtered_df.to_numpy()

# plt.errorbar( range(len(data1[:,0])),data1[:,0],yerr=data1[:,1],fmt='.-')
plt.errorbar( range(len(data[:,0])),data[:,10],yerr=data[:,11],fmt='.-')
plt.show()

# def plotspec(data):
#     y = spectral(data[:,9],data[:,11],data[:,5],data[:,4],data[:,3],0.125)
    
#     plt.subplot(1,2,1)
#     plt.errorbar(data[:,5],y,yerr=data[:,10]*0,fmt='x-')
#     plt.subplot(1,2,2)
#     plt.errorbar(data[:,5],data[:,11],yerr=data[:,12],fmt='x-')
#     plt.show()




# plotspec(data)