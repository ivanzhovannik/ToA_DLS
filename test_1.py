'''
that is a test file for the TTTCU algorithm
'''
import numpy as np
import matplotlib.pyplot as plt


def GenerateLagTimeVect(Tmax, B = 2, lag_time_min = 3e3, ToA_st = 1):
    '''
    The function generates a lag time vector for a correlation function 
    Tmax = upper correlation time
    B =  integer base number
    '''
    
    tau_ind_lim = 500 # limit the index
    tauRaw = np.zeros((500,))
    N = lag_time_min / ToA_st
    tauRaw[0] = N
    j = 0 # the while loop index
    
    # maximal lag time must be a way smaller than Tmax
    while tauRaw[j] < Tmax / 100:
        j += 1
        tauRaw[j] = tauRaw[j - 1] + N * 2 ** ((j - 1) // B)
    return tauRaw[:j + 1]


def TTTCU(dToA, ToA_st = 1, lag_time_min = 3e3):
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
    tauVect = 
    return 1
