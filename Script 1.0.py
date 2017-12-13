# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 17:08:51 2017

@author: User
"""

# import modules
import oct2py # to run code in MATLAB
import numpy as np # MATLAB-like data and operations
import pandas as pd # fancy excel-like tables
import matplotlib.pyplot as plt # plot figures as in MATLAB
from scipy.signal import savgol_filter as sgf

# initialize oct2py to execute code in MATLAB
oc = oct2py.Oct2Py()

# initialize all dictionaries
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
num_of_samples = 4

# initialize samples and its keys
'''
sample = {1: '5.txt', 2: '2.txt', 3: '3.txt', 4: '4.txt',
          5: '5.txt', 6: '6.txt', 7: '7.txt', 8: '8.txt',
          9: '9.txt', 10: '10.txt'}
'''
sample = {1: '2.txt', 2: '3.txt', 3: '6.txt', 4: '7.txt'}

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

for i in range(1, num_of_samples + 1): # for each of the samples
    
    # read the data to matlab
    data[i] = pd.read_csv(''.join(['au 13.12.17 1/au 0', sample[i]]), names=['ToA'])
    
    
    I[i] = 1e9 / data[i]['ToA'] # ToA in [ns], transform it to intencities
    ACF_base_I[i] = np.mean(I[i]) ** 2 # ACF base line option #1
    
    # run code in MATLAB
    dToA = 1 # [ns] for ToA from device
    dN = 100 # [ns] for intencities    
    ACF_[i] = oc.TTTCU(data[i]['ToA'].as_matrix(), dToA, dN) # execute
    ACF[i] = np.c_[ACF_[i][:ACF_[i].shape[0] // 2, 0], ACF_[i][ACF_[i].shape[0] // 2:, 0]]
    
    # cut the values with lag times less than 500 ns
    ACF_cut[i] = ACF[i][np.where(ACF[i][:, 0] > 1000)]
    
    # filter the curves with SavGol filter
    ACF_fil[i] = np.c_[ACF_cut[i][:, 0], ACF_cut[i][:, 1]] # sgf(ACF_cut[i][:, 1], 35, 5)]
    
    # find the first derivatives of the ACF curves
    ACF_fil_der[i] = np.gradient(ACF_fil[i][:, 1])
    
    # compute the base line: SOME ERRORS OCCURES
    ACF_base_der[i] = ACF_fil[i][np.argmin(
            (np.mean(
                    ACF_fil[i][np.where(np.diff(np.sign(ACF_fil_der[i]))), 0] # average all the time points where derivative changes its sign
                          ) - ACF_fil[i][:, 0]) ** 2 # find the index of the closest point to one was found
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
