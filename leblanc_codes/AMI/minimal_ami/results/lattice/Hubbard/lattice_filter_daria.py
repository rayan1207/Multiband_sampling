import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file_path = "result_fullo1234_U1.txt"
# data = np.loadtxt(file_path)

U =4
data_U = {
    'Type_of_U': [0, 1],
    'U': [U, U]
}

df_U = pd.DataFrame(data_U)
df = pd.read_csv(file_path, delim_whitespace=True, skiprows=1, header=None)

ext_vars=np.loadtxt('../../../loader/ext_vars.dat')


for i in range(4):
    df = df.merge(df_U, left_on=14 + i, right_on='Type_of_U', how='left')
    df.fillna(1, inplace=True)
    # print(df)
    df[10] = df[10] * df['U']
    df[11] = df[11] * df['U']
    df[12] = df[12] * df['U']
    df[13] = df[13] * df['U']
    df = df.drop(columns=['U', 'Type_of_U'])



# Group by the first 10 columns and sum the values in columns 9 and 10

grouped_df = df.groupby(df.columns[1:10].tolist())[[10, 11, 12, 13]].sum().reset_index()
# grouped_df = df.groupby(df.columns[:10].tolist())[[10, 11, 12, 13]].sum().reset_index()

filtered_df = grouped_df[ (grouped_df[8] == 1) & (grouped_df[9] == 1)  ]
# filtered_df = grouped_df[(grouped_df[8] == 1) & (grouped_df[9] == 1)  & (grouped_df[0] == 4)]
# filtered_df = filtered_df.sort_values(by=[ext_vars[:, 5], ext_vars[:, 6]], axis=0)


def symmetrize(filtered_df):
    filtered_df1 = filtered_df[(filtered_df[5] == 0.000)]
    filtered_df2 = filtered_df[(filtered_df[4] == 3.141590) & (filtered_df[5] > 0.00)]
    filtered_df3 = filtered_df[(filtered_df[5] == filtered_df[4]) & (filtered_df[5] < 3.14)].iloc[::-1]
    
    symmetry_cut_df = pd.concat([filtered_df1, filtered_df2, filtered_df3])

    return symmetry_cut_df

filtered_df = symmetrize(filtered_df)
filtered_df.to_csv('filtered_output.txt', sep=' ', header=False, index=False, float_format='%.7f')

data= filtered_df.to_numpy()
plt.subplot(1,2,1)
plt.plot(data[:,9]+U/2,'x-')
plt.subplot(1,2,2)
plt.plot(data[:,11],'x-')
plt.show()


result = np.column_stack((np.asarray(data[:,9]+U/2), np.asarray(data[:,10]),np.asarray(data[:,11]), np.asarray(data[:,12])))
# np.savetxt('result_o3_U4_v0.txt',result)



# grouped_df.to_csv('data.txt', sep=' ', header=False, index=False, float_format='%.7f')
# grouped_df_ord.to_csv('data_ord.txt', sep=' ', header=False, index=False, float_format='%.7f')


# print(grouped_df)
