## Computes the response variabilty as the correlation between spike trains adapted from Galan et al., 2008
#       INPUT
#       dap         vector of structures pointing to the different files obtained from 'excel_read'


#       OUTPUT
#       sCorr       matrix of correlation values
#       sCorrVec    Vector of the relevant values from sCorr
#       spBins
#       R           results object that collects all analysis and will be saved

import os
import pickle
import math
import numpy as np
# import matlab.engine
from neo import io
from sl_sync_params import *
from spikeCorr import *
from doCalc import *

# eng = matlab.engine.start_matlab()

def respVariability(dap, excelName, doplot = True, dosave = False):

    dir = dap[0]['Dir'] + '/'

    # dap = [18125012, 18125013]    # vector which contains the files of cells that need to be analysed

    ap = sl_sync_params()     # dictionary containing all of default analysis parameters - consider making ap it's own class object with attributes

    nFiles = len(dap)

    files = []
    for i in dap:
        files.append(dap[i]['fname'][0][:-1])

    # get voltage and current trace
    print('Loading data...')
    h = {}
    si = {}     # sampling intervals for each cell in us
    d = {}
    V = {}
    I = {}
    for f in files:
        fpath = dir + str(f) + '.abf' # combine file name with '.abf'
        if str(f) + '.abf' in os.listdir(dir):
            a = io.AxonIO(filename=fpath)
            h[f] = a.read_header()
            bl = a.read_block(lazy=False, cascade=True)
                # - .segments represent sweeps (one segment = one sweep)
                # - .analogsignals for each segment: numpy array of voltage recordings and current input,
                #   length of recording block, and sampling rate
            vchan = list(h[f]['nADCSamplingSeq']).index(4) # set channel 4 as voltage channel
            ichan = list(h[f]['nADCSamplingSeq']).index(14) # set channel 14 as command channel
            V[f] = []  # numpy array of voltage recordings for all sweeps/segments - rows = sweeps, columns = data
            for i in range(0, len(bl.segments)):
                a = bl.segments[i].analogsignals[vchan].__array__().tolist()
                V[f].append([item for x in a for item in x])
            V[f] = np.array(V[f])
            I[f] = []  # numpy array of voltage recordings for all sweeps/segments - rows = sweeps, columns = data
            for i in range(0, len(bl.segments)):
                a = bl.segments[i].analogsignals[ichan].__array__().tolist()
                I[f].append([item for x in a for item in x])
            I[f] = np.array(I[f])

            # save data block for each cell
            d[f] = bl

            # sampling rate for each sweep
            si[f] = h[f]['fADCSampleInterval']*h[f]['nADCNumChannels']
        else:
            print('File does not exist: %s' % fpath)

    # Make sure all cells are sampled at the same rate - otherwise it gets a bit messy
    nsi = []
    for i in si.itervalues():
        nsi.append(i)
    if len(set(nsi)) > 1:
        raise ValueError('Data acquired at different sampling rates')
    si = nsi[0]

    # compute some constants
    siSec = si*1e-6
    fs = 1/siSec
    wPoints = math.floor((ap['rv']['window']/1000.)/siSec)

    STA = {}; G = {}; phi = {}; freq = {}; spBins = {}
    for f in files:
        print('Cell #%s: Performing analysis...' % f)
        sta = []; g = []; PHI = []; FREQ = []; SPbINS = []
        nTrials = len(d[f].segments)
        for i in range(0, nTrials):
            V_1 = V[f][i]
            I_1 = I[f][i]
            STA_, G_, phi_, freq_, spBins_ = doCalc(V_1, I_1, si, ap)
            sta.append(STA_); g.append(G_); PHI.append(phi_); FREQ.append(freq_); SPbINS.append(spBins_)
    STA[f] = sta; G[f] = g; phi[f] = PHI; freq[f] = FREQ; spBins[f] = np.array(SPbINS)

    # repackage all the results after the for loop - I dont think this is neccessary yet - work on plotting and saving functions, and see what exactly is needed for this part here
    # STA_a = []
    # c = 0
    # for i in range(0, nTrials):
    #     spBins_a[i, :] = spBins[i]
    #     if len(G[i]) > 0:
    #         c += 1
    #         G_a[c,:] = G[i]
    #         phi_a[c,:] = phi[i]
    #         STA_a = np.transpose(np.array(np.transpose(STA_a), np.transpose(STA[i])))
    #


    # spBins = np.array(spBins)
    # G = G_a
    # phi = phi_a
    # STA = STA_a
    # freq = freq[0]

    # compute the average spike rate across all trials
    for i in spBins:
        avgSpikeRate = sum(sum(spBins[i]))/nTrials*ap['rv']['length']
        print('Average Spike Rate: %d Hz' % avgSpikeRate)

    # If the average spike rate is less than an arbitrary amount send a warning and dont compute the correlations

    if avgSpikeRate < ap['rv']['minSpikeRate']:
        print('Average spike rate too low to perform analyses')

    # Do spike correlations
    sCorr = {}; sCorrVec = {}
    for f in files:
        sCorr_, sCorrVec_ = spikeCorr(spBins[f], fs/wPoints, ap['rv']['delta'])        # note that the pearson correlation values are slightly discordant between Matlab and python
        sCorr[f] = sCorr_; sCorrVec[f] = sCorrVec_

    # save results
    pickle.dump([sCorr, sCorrVec, V, I, spBins, dap, fs, wPoints, STA,
                 avgSpikeRate, G, phi, freq],
                open(str(os.path.join(dap[0]['Dir'], '%s_results.pkl' % excelName)), 'wb'))


    return sCorr, sCorrVec, spBins #, R
##
