# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 20:09:21 2023

@author: Rayan
"""

import pandas as pd

# Read the data into a DataFrame
df = pd.read_csv('output.txt', delim_whitespace=True, skiprows=1, header=None)

# Exclude the first row

# Add a new column with original row numbers
df['Row Number'] = df.index

# Sort the DataFrame by the absolute value of column 12 and then column 14
sorted_df = df.iloc[df[[11]].abs().max(axis=1).argsort()]
print(sorted_df)

# Write the sorted DataFrame to a file
sorted_df.to_csv('read_max.txt',sep=' ', header=False, index=False, float_format='%.7f')
