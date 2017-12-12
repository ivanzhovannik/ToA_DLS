# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 17:08:51 2017

@author: User
"""
import oct2py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter as sgf

oc = oct2py.Oct2Py()
data = {}

I = {}
ACF_ = {}
ACF = {}
ACF_cut = {}
ACF_fil = {}
ACF_fil_der = {}
ACF_base_der = {}
ACF_base_I = {}
ACF_fil_base ={}
num_of_samples = 1
sample = {1: '5.txt', 2: '2.txt', 3: '3.txt', 4: '4.txt',
          5: '5.txt', 6: '6.txt', 7: '7.txt', 8: '8.txt',
          9: '9.txt', 10: '10.txt'}

# use the following to merge the data from several files to one pd.DataFrame
'''
data_ = {} # to merge the data
num_of_samples_merge = 2
for i in range(1, num_of_samples_merge + 1):
    data_[i] = pd.read_csv(''.join(['09.12.17/sampling_time', sample[i]]),
                           names=['ToA'])
    if i == 2:
        new_data = np.r_[data_[i-1]['ToA'], data_[i]['ToA']]
    elif i >= 2:
        new_data = np.r_[new_data, data_[i]['ToA']]
data[1] = pd.DataFrame(data = new_data, columns=['ToA'])
'''
# skip the lines above to skip merging

for i in range(1, num_of_samples + 1):
    data[i] = pd.read_csv(''.join(['09.12.17/sampling_time', sample[i]]), names=['ToA'])
    I[i] = 1e9 / data[i]['ToA']
    ACF_base_I[i] = np.mean(I[i]) ** 2
    dToA = 1 # ns
    dN = 100 # ns
    ACF_[i] = oc.TTTCU(data[i]['ToA'].as_matrix(), dToA, dN)
    ACF[i] = np.c_[ACF_[i][:ACF_[i].shape[0] // 2, 0], ACF_[i][ACF_[i].shape[0] // 2:, 0]]
    
    # cut the values with lag times less than 500 ns
    ACF_cut[i] = ACF[i][np.where(ACF[i][:, 0] > 1000)]
    
    # filter the curves with SavGol filter
    ACF_fil[i] = np.c_[ACF_cut[i][:, 0], sgf(ACF_cut[i][:, 1], 35, 5)]
    
    # find the first derivatives of the ACF curves
    ACF_fil_der[i] = np.gradient(ACF_fil[i][:, 1])
    ACF_base_der[i] = ACF_fil[i][np.argmin(
            (np.mean(
                    ACF_fil[i][np.where(np.diff(np.sign(ACF_fil_der[i]))), 0]
                          ) - ACF_fil[i][:, 0]) ** 2
                                      ), 1]
    
    # cut the baselines 
    ACF_fil_base[i] = ACF_fil[i]
    ACF_fil_base[i][:, 1] = ACF_fil_base[i][:, 1] / ACF_base_der[i] - 1



plt.figure(figsize=(15, 10))
# plt.xlim([1e3, np.max(ACF[1][:, 0])])
# plt.ylim([1000, 3000])
for i in range(1, num_of_samples + 1):
    # plt.semilogx(ACF_cut[i][:, 0], ACF_cut[i][:, 1])
    # plt.semilogx(ACF_fil[i][:, 0], ACF_fil[i][:, 1])
    plt.semilogx(ACF_fil_base[i][:, 0], ACF_fil_base[i][:, 1])
plt.savefig('fig3.png')
print(ACF_base_I, ACF_base_der, sep='\n')

for i in range(1, num_of_samples + 1):
    data_to_write = pd.DataFrame(data = ACF_fil_base[i],
                                 columns=['lag_time', 'ACF'])
    data_to_write.to_excel(''.join(['result ', str(i),'.xlsx']),
                 sheet_name=str(i),
                 index=False)
