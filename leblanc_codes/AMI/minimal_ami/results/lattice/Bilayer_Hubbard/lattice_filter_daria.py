import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
file_path = "output.txt"
# data = np.loadtxt(file_path)
U = 3
V = 0
J = 0

def symmetrize(filtered_df):
    filtered_df1 = filtered_df[(filtered_df[5] == 0.000)]
    filtered_df2 = filtered_df[(filtered_df[4] == 3.141590) & (filtered_df[5] > 0.00)]
    filtered_df3 = filtered_df[(filtered_df[5] == filtered_df[4]) & (filtered_df[5] < 3.14)].iloc[::-1]
    
    symmetry_cut_df = pd.concat([filtered_df1, filtered_df2, filtered_df3])

    return symmetry_cut_df



Uindex= np.arange(0,12,1)
Uvalues = np.append(np.ones(4)*U, np.ones(8)*V)


print(Uvalues)
data_U = {
    'Type_of_U': Uindex,
    'U': Uvalues
}



df_U = pd.DataFrame(data_U)
df = pd.read_csv(file_path, delim_whitespace=True, skiprows=1, header=None)

for i in range(4):
    df = df.merge(df_U, left_on=14 + i, right_on='Type_of_U', how='left')
    df.fillna(1, inplace=True)
    print(df)
    df[10] = df[10] * df['U']
    df[11] = df[11] * df['U']
    df[12] = df[12] * df['U']
    df[13] = df[13] * df['U']
    df = df.drop(columns=['U', 'Type_of_U'])

#print(df)

# Group by the first 10 columns and sum the values in columns 9 and 10

grouped_df = df.groupby(df.columns[1:10].tolist())[[10, 11, 12, 13]].sum().reset_index()
grouped_df_ord = df.groupby(df.columns[:10].tolist())[[10, 11, 12, 13]].sum().reset_index()

filtered_df = grouped_df[(grouped_df[8] == 1) & (grouped_df[9] == 1)]

filtered_df = symmetrize(filtered_df)
print(filtered_df)

filtered_df.to_csv('filtered_output.txt', sep=' ', header=False, index=False, float_format='%.7f')



data= filtered_df.to_numpy()
plt.plot(data[:,9],'x-')
plt.show()


# grouped_df.to_csv('data.txt', sep=' ', header=False, index=False, float_format='%.7f')
# grouped_df_ord.to_csv('data_ord.txt', sep=' ', header=False, index=False, float_format='%.7f')


# print(grouped_df)