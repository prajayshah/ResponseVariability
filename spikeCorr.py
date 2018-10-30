import numpy as np
import math
from scipy.stats import pearsonr


def spikeCorr(sTimes, fs, delta, tList=None):

    '''
    :param sTimes:      Numpy array of spike times coded as 0s and 1s
    :param fs:          Sampling rate
    :param delta:       Window function in ms as per Galan et al 2008
    :param tList:
    :return: sCorr:     Spike correlations between all possible pairs in Matrix
    :return: sCorrVec:  Spike train correlations as a vector
    '''

    if tList is None:
        tList = []

    if type(sTimes) == list:
        nTrials = len(sTimes)
        nPoints = len(sTimes[0])
    else:
        nTrials, nPoints = sTimes.shape

    dlength = math.floor((delta/1000.)*fs)
    dfunc = [1.] * int(dlength)
    sconv = [[0.] * nPoints] * nTrials
    sCorrVec = []

    for i in range(0, nTrials):
        sconv[i] = np.convolve(sTimes[i], dfunc, 'same')

    if len(tList) == 0:
        c = 0
        sCorr = np.array([[0.] * nTrials] * nTrials)
        sCorrVec = [0.] * int(nTrials*(nTrials-1)/2)
        for i in range(0, nTrials):
            for j in range((i+1), nTrials):
                r, prob = pearsonr(sconv[i], sconv[j])
                sCorr[i, j] = r
                sCorrVec[c] = sCorr[i][j]
                c += 1
    else:
        # more than one cell passed to the function
        tListLength = len(tList)
        for i in range(0, tListLength):
            for j in range((i+1), tListLength):
                if i == 0:
                    start_1 = 0
                else:
                    start_1 = tList[i] + 1
                end_1 = start_1+tList[i]-1
                start_2 = sum(tList[0:(j-1)])+1
                end_2 = start_2+tList[j]-1

                sconv1 = sconv[start_1:end_1, :]
                sconv2 = sconv[start_2:end_2, :]
                sCorr = corr2cells(sconv1, sconv2)

    return sCorr, sCorrVec


def corr2cells(sconv1, sconv2):
    nTrials1 = len(sconv1[0])
    nTrials2 = len(sconv2[0])

    sCorr = [[0]*nTrials1] * nTrials2

    for i in range(0, nTrials1):
        for j in range(0, nTrials2):
            sCorr[i, j] = pearsonr(sconv1[i, :].T, sconv2[j, :].T)

    return sCorr