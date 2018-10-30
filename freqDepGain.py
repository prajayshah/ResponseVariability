import numpy as np
from scipy.signal import correlate
from goertzel import *
from detect_peaks import *

#fs = 1/siSec

def freqDepGain(V, I, fs, ap):
    nPoints = len(V)
    r = np.array([0]*nPoints)
    spInd = detect_peaks(V, mph=ap['rv']['threshold'])
    r[spInd] = fs

    # compute the correlations as per eqs 10 & 11 from Higgs and Spain, 2008, J Neurosci
    csr = correlate(r, I, mode='full')
    tau = np.argmax(csr)    ## TROUBLESHOOT!!!! not same as matlab for some reason...
    css = np.correlate(I, I, mode='full')
    css = css[:len(css) / 2]    # use only first half of correlation? - TROUBLESHOOT!!

    freq = np.arange(1., ap['rv']['frange'], 0.1*np.log(10))

    N = len(csr)    # length of the CSR function

    # calculate the dft for each specific frequency
    CSS = np.zeros(shape = (len(freq), 3))
    CSR = np.zeros(shape = (len(freq), 3))

    for i in range(0, len(freq)):
        w = wind(freq[i], tau/fs)

        windowed_css = css*w # select a certain window by multiplying a gaussian window by the data window
        windowed_csr = csr*w # select a certain window by multiplying a gaussian window by the data window
        freqInd = round(freq[i]/fs*N)+1 # frequency index is the input to the goertzel function
        CSS[i] = np.array(list(goertzel(windowed_css, fs, (int(freqInd), int(freqInd)))[1][0])) * 2./N # get the css dft for this specific frequency
        CSR[i] = np.array(list(goertzel(windowed_csr, fs, (int(freqInd)-fs/len(csr), int(freqInd)+fs/len(csr)))[1][2])) * 2./N # get the csr dft for this specific frequency


def wind(f, tau):
    w = np.exp((-np.power(f,2)*np.power(tau,2))/2)

    return w

