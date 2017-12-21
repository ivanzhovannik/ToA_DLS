'''
that is a file for the TTTCU algorithm
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tkinter import Tk, filedialog # dialog windows


def GenerateLagTimeVect(Tmax, B=2, lag_time_min=3e3, ToA_st=1):
    '''
    The function generates a lag time vector for a correlation function 
    Tmax = upper correlation time
    B =  integer base number
    '''
    
    tau_ind_lim = 500 # limit the index
    tauRaw = np.zeros((tau_ind_lim,))
    N = lag_time_min / ToA_st
    tauRaw[0] = N
    j = 0 # the while loop index
    
    # maximal lag time must be a way smaller than Tmax
    while tauRaw[j] < Tmax / 100:
        j += 1
        tauRaw[j] = tauRaw[j - 1] + N * 2 ** ((j - 1) // B)
    return tauRaw[:j + 1]


def TTTCU(dToA, ToA_st=1, lag_time_min=3e3):
    '''
    This function generates correlation function from a vector of
    single-photon arrival times (ToA) in nanoseconds [ns]
    
    ToA - vector of ToA in differential form
    ToA_st - sampling time for ToA (it is usually equal to 1 ns) [ns]
    lag_time_min - minimal lag time for the correlation function [ns]
    '''
    
    # check the types of the arguments
    assert(np.ndarray == type(dToA))
    
    ToA_shape = dToA.shape
    ToA = np.zeros(ToA_shape)
    calc = 0
    for i, elem in enumerate(dToA):
        calc += elem
        ToA[i] = calc
    tauVect = GenerateLagTimeVect(ToA[-1],
                                  lag_time_min=lag_time_min,
                                  ToA_st=ToA_st)
    
    # initialize a vector for the correlation function
    ACF = np.zeros(tauVect.shape)
    
    # calculate the ACF for each of the lag times
    for i, tau in enumerate(tauVect):
        ToA_shifted = ToA + tau
        ACF[i] = ToA_shape[0] - 0.5 * len(np.unique(np.r_[ToA, ToA_shifted]))
    return np.c_[tauVect, ACF]


def ToA2ACF(ToA_st=1, lag_time_min=3e3, sav_gol_param=None):
    '''
    This function takes several ToA files customly chosen from a folder and
    transforms them to ACF via TTTCU algorithm, cuts the ACFs' baselines, and
    averages them (with possible filtration)
    
    Inputs:
        ToA_st - sampling time for ToA (it is usually equal to 1 ns) [ns]
        lag_time_min - minimal lag time for the correlation function [ns]
        sav_gol_param - Savitzky-Golay filter parameters dictionary (see below)
    
    Outputs: 
        ACF - output averaged correlation function
        data - ToAs dictionary
    
    Importantly, ToAs should be stored in a txt file WITHOUT any header.
    
    Savitzki-Golay filter parameters should be stored in a dictionary as
    follows:
        {'poly': polynomial_order, 'window': window_size}
        
    The exp ** -1 decay for the ACF of 1 nm particles starts for lag times more
    than 4 miliseconds. Hence, it is important to find a base lines only for
    these lag times (more than 4 ms). We use 1 ms to have a safety factor.
    '''
    
    # open a dialog window and select some text files with ToAs
    root = Tk()
    root.fileNames = filedialog.askopenfilenames(initialdir="/",
                                                 title="Select file",
                                                 filetypes=(("ToA text files",
                                                             "*.txt"),
        ("all files","*.*"))) # get the file names
    fileNames = root.fileNames # store the file names
    root.destroy() # close the browser
    
    # convert fileNames to an array of file names
    if type(fileNames) == str:
        fileNames = list(fileNames)
    elif type(fileNames) != list and type(fileNames) != tuple:
        raise TypeError('"fileNames" must be a list of file names or a string')
    
    # data storages
    data = {}
    ACF = {}
    ACF_cut = {}
    ACF_fil_der = {}
    ACF_base_der = {}
    ACF_base = {}
    ACF_res_y = 0
    
    # apply the algorithm to each of the files
    for i, fname in enumerate(fileNames):
        
        # read the data from the current file
        data[i] = pd.read_csv(fname, names=['ToA'])
        
        # run the algorithm: the single ACF calculation
        ACF[i] = TTTCU(data[i]['ToA'].as_matrix(),
                               ToA_st=ToA_st,
                               lag_time_min=lag_time_min)
        plt.figure()
        plt.plot(ACF[i][:, 0], ACF[i][:, 1])
        
        # cut the values with lag times less than 3000 ns
        ACF_cut[i] = ACF[i][np.where(ACF[i][:, 0] > 3000)]
        
        # find the first derivatives of the ACF curves
        ACF_fil_der[i] = np.gradient(ACF_cut[i][:, 1])
        
        # compute the base line
        ACF_base_der[i] = ACF_cut[i][np.argmin(
            (np.mean(
                    ACF_cut[i][np.where(np.diff(np.sign(ACF_fil_der[i]))) and # average all the time points where derivative changes its sign
                           np.where(ACF_cut[i][:, 0] > 1e6), 0] # baseline cannot be found for lag times less the 1 ms
                    ) - ACF_cut[i][:, 0]) ** 2 # find the index of the closest point to one was found
                                                ), 1]
        
        # cut the baselines 
        ACF_base[i] = ACF_cut[i]
        ACF_base[i][:, 1] = ACF_base[i][:, 1] / ACF_base_der[i] - 1

        # for the future averaging
        if i != 0:
            indmax = np.min([ACF_res_y.shape[0],ACF_base[i].shape[0]])
            ACF_res_y = ACF_base[i][:indmax, 1] + ACF_res_y[:indmax]
        else:
            ACF_res_y = ACF_base[i][:, 1]
    
    # average ACF by dividing to the number of files
    ACF_res = np.c_[ACF_base[0][:, 0], ACF_res_y / (i + 1)]
        
    return ACF_res, data
