# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 16:55:16 2017

@author: User
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def count_size(lag_time, ACF, device='arn'):
    lambd = 632.8 * 1e-9
    theta_mal = 173 * np.pi / 180
    theta_arn = 90 * np.pi / 180
    n_water = 1.33
    q_mal = 4 * np.pi * n_water * np.sin(theta_mal / 2) / lambd
    q_arn = 4 * np.pi * n_water * np.sin(theta_arn / 2) / lambd

    k_B = 1.38 * 1e-23
    visc = 8.89 * 1e-4
    T_C = 27
    T = 273 + T_C
    
    help_arn = (q_arn ** 2) * k_B * T / (6 * np.pi * visc)
    help_mal = (q_mal ** 2) * k_B * T / (6 * np.pi * visc)
    
    weight = np.shape(lag_time)[0] * (ACF ** 2) / np.sum(ACF ** 2)
    Cum = np.polyfit(lag_time, np.log(np.abs(ACF + 1e-240)), 2, w=weight ** 2)
    
    
    if device == 'arn':
        size = -2 * help_arn / Cum[1]
        PdI = 2 * Cum[0] / (Cum[1] ** 2)
    elif device == 'mal':
        size = -2 * help_mal / Cum[1]
        PdI = 2 * Cum[0] / (Cum[1] ** 2)
    return {'size': size, 'PdI': PdI, 'P': Cum}


# malvern data
ACF_mal = pd.read_excel('mal 01.11.17.xlsx', sheetname=['ACF'])
ACF_mal['ACF']['X Lag Time'] *= 1e3

# data dicts
ACF_arn = {}
ACF_zeta = {}
ACF_arn_I = {}
ACF_zeta_I = {}

num_of_samples = 5

plt.figure(figsize=(20, 15))
for i in range(1, num_of_samples + 1):
    
    # ToA ACFs
    ACF_arn[i] = pd.read_excel(''.join(['ARN ACF ToA 01.11.17/result ',
                                       str(i),
                                       '.xlsx']))
    ACF_zeta[i] = pd.read_excel(''.join(['ZetaDLS ACF ToA 01.11.17/result ',
                                       str(i),
                                       '.xlsx']))
    
    # intensity ACFs
    ACF_arn_I[i] = pd.read_csv(''.join(['ARN/I ', str(i), '.txt']),
                             delimiter=' ', names=['lag_time', 'ACF'])
    ACF_zeta_I[i] = pd.read_csv(''.join(['ZetaDLS/I ', str(i), '.txt']),
                             delimiter=' ', names=['lag_time', 'ACF'])
    ACF_arn_I[i]['lag_time'] *= 1e3
    ACF_zeta_I[i]['lag_time'] *= 1e3
    
    plt.subplot(2, 3, i)
    plt.semilogx(ACF_zeta[i]['lag_time'], ACF_zeta[i]['ACF'] / np.max(ACF_zeta[i]['ACF']), 'b')
    plt.semilogx(ACF_arn[i]['lag_time'], ACF_arn[i]['ACF'] / np.max(ACF_arn[i]['ACF']), 'y')
    
    plt.semilogx(ACF_zeta_I[i]['lag_time'], ACF_zeta_I[i]['ACF'] / np.max(ACF_zeta_I[i]['ACF']), ':b')
    plt.semilogx(ACF_arn_I[i]['lag_time'], ACF_arn_I[i]['ACF'] / np.max(ACF_arn_I[i]['ACF']), ':y')
    
    plt.semilogx(ACF_mal['ACF']['X Lag Time'],
                 ACF_mal['ACF'][''.join(['Record ',
                                  str(144 + i),
                                  ': Au ', str(i)])] / 
                np.max(ACF_mal['ACF'][''.join(['Record ',
                                  str(144 + i),
                                  ': Au ', str(i)])]))
    
    plt.xlabel('lag time, ns', {'fontsize': 14})
    plt.title(' '.join(['SAMPLE ', str(i)]), {'fontsize': 14})
    print('\n'.join(['Zeta_ToA:',
                                   str(count_size(ACF_zeta[i]['lag_time'],
                                                  ACF_zeta[i]['ACF'],
                                                  device='arn')),
                                   'ARN_ToA:',
                                   str(count_size(ACF_arn[i]['lag_time'],
                                                  ACF_arn[i]['ACF'],
                                                  device='arn')),
                                   'Zeta_I:',
                                   str(count_size(ACF_zeta_I[i]['lag_time'],
                                                  ACF_zeta_I[i]['ACF'],
                                                  device='arn')),
                                   'ARN_I:',
                                   str(count_size(ACF_arn_I[i]['lag_time'],
                                                  ACF_arn_I[i]['ACF'],
                                                  device='arn')),
                                   'Mal:',
                                   str(count_size(ACF_mal['ACF']['X Lag Time'],
                                                  ACF_mal['ACF'][''.join(['Record ',
                                                         str(144 + i),
                                                         ': Au ', str(i)])],
                                                  device='mal'))]))
    print('\n\n')
plt.savefig('ACFs.png')



