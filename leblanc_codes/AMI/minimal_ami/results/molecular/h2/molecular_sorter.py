# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 02:44:11 2023

@author: Rayan
"""

import pandas as pd
import glob
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
dfs = []

current_dir = os.getcwd()

for filename in os.listdir(current_dir):
    if filename.endswith('.txt'):
        file_path = os.path.join(current_dir, filename)
        df = pd.read_csv(file_path, delimiter=' ', header=None)
        dfs.append(df)

merged_df = pd.concat(dfs, ignore_index=False)
# merged_df = merged_df[(merged_df[4] == 2) & (merged_df[5] == 2)]
grouped_df = merged_df.groupby(df.columns[:2].tolist())[[2,3]].sum().reset_index()

data= grouped_df.to_numpy()

plt.subplot(1,2,1)
plt.plot(data[:,1],data[:,2])
plt.subplot(1,2,2)
plt.plot(data[:,1],data[:,3])
plt.show()

