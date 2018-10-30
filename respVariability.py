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
from neo import io
from sl_sync_params import *
from spikeCorr import *
from doCalc import *


def respVariability(dap, excelName, doplot = True, dosave = False):


    # set directory for these cells to be the directory of the first cell (since the rest of the cells should also be in the same directory)
    dir = dap[list(dap.keys())[0]]['Dir'] + '/'

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
            # h[f] = a.read_header()
            bl = a.read_block(lazy=False, signal_group_mode='split-all', units_group_mode='split-all')
                # - .segments represent sweeps (one segment = one sweep)
                # - .analogsignals for each segment: numpy array of voltage recordings and current input,
                #   length of recording block, and sampling rate
            iclamp = 0 # channel 4 as voltage channel
            current_in = 1 # channel 14 as command channel
            V[f] = []  # numpy array of voltage recordings for all sweeps/segments - rows = sweeps, columns = data
            for i in range(0, len(bl.segments)):
                a = bl.segments[i].analogsignals[iclamp].__array__().tolist()
                V[f].append([item for x in a for item in x])
            V[f] = np.array(V[f])


            if 'Offset' in dap[f]['cond.names'][0]:
                V[f] = V[f] + dap[f]['cond.times']

                # compute and report average RMP
                num_pts_RMP = int(ap['rv']['prefix'] * 10000)
                nTrials = len(V[f])
                x = 0
                for i in range(0, nTrials):
                    trial_average = np.average(V[f][i][:num_pts_RMP])
                    x += trial_average
                cell_average = x / nTrials
                print('Average RMP %s: %d' % (f, cell_average))

            else:
                # compute the average RMP - if no offset specified
                num_pts_RMP = ap['rv']['prefix'] * 10000
                num_pts_RMP = int(num_pts_RMP)
                nTrials = len(V[f])
                x = 0
                for i in range(0, nTrials):
                    trial_average = np.average(V[f][i][:num_pts_RMP])
                    x += trial_average
                cell_average = x / nTrials
                # adjust based on RMP specified in excel file
                # retrieve RMP from excel file (stored in dap) and calculate offset required
                offset = float(dap[f]['sf']['Tags'][dap[f]['sf']['Tags'].index('RMP')+4:dap[f]['sf']['Tags'].index('RMP')+9]) - cell_average
                V[f] = V[f] + offset

                # calculate and report new cell RMP average
                x = 0
                for i in range(0, nTrials):
                    trial_average = np.average(V[f][i][:num_pts_RMP])
                    x += trial_average
                cell_average = x / nTrials
                print('Average RMP %s: %d' % (f, cell_average))


            I[f] = []  # numpy array of stimulus for all sweeps/segments - rows = sweeps, columns = data
            for i in range(0, len(bl.segments)):
                a = bl.segments[i].analogsignals[current_in].__array__().tolist()
                I[f].append([item for x in a for item in x])
            I[f] = np.array(I[f])

            # save data block for each cell
            d[f] = bl

            # sampling interval for each sweep in us
            si[f] = 1e6/10000 # sampling rate was 10000 Hz
            #si[f] = h[f]['fADCSampleInterval']*h[f]['nADCNumChannels']
        else:
            print('File does not exist: %s' % fpath)

    # Make sure all cells are sampled at the same rate - otherwise it gets a bit messy
    nsi = []
    for i in si.values():
        nsi.append(i)
    if len(set(nsi)) > 1:
        raise ValueError('Data acquired at different sampling rates')
    si = nsi[0]

    # compute some constants
    siSec = si*1e-6 # sampling interval in seconds
    fs = 1/siSec # sampling rate (Hz)
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




    # compute the average spike rate across all trials
    avgSpikeRate = {}
    for f in files:

        spikecount = 0
        for i in range(0, len(spBins[f])):
            spikecount += sum(spBins[f][i])

        avgSpikeRate[f] = spikecount / (len(V[f]) * ap['rv']['length'])
        print('Average Spike Rate: %d Hz' % avgSpikeRate[f])

    sCorr = {}; sCorrVec = {}
    for f in files:
        # If the average spike rate is less than an arbitrary amount send a warning and dont compute spike correlation
        if avgSpikeRate[f] < ap['rv']['minSpikeRate']:
            print('Average spike rate too low to perform analyses')
        else:
            # Do spike correlations
            sCorr_, sCorrVec_ = spikeCorr(spBins[f], fs/wPoints, ap['rv']['delta'])        # note that the pearson correlation values are slightly discordant between Matlab and python
            sCorr[f] = sCorr_; sCorrVec[f] = sCorrVec_

    # save results
    pickle.dump([sCorr, sCorrVec, V, I, spBins, dap, fs, wPoints, STA,
                 avgSpikeRate, G, phi, freq],
                open(str(os.path.join(dap[[*dap][0]]['Dir'], '%s_results.pkl' % excelName)), 'wb'))


    return sCorr, sCorrVec, spBins #, R
##
