from pynwb import NWBHDF5IO

import math
import numpy as np
import pickle
import pandas as pd
import datetime

from sl_sync_params_v2 import *
from spikeCorr import *
from doCalc import *
from rvCompare import *

def readNWBpatchClamp(fpath):

    # read nwb file for the chosen file
    io = NWBHDF5IO(fpath, 'r')
    nwbfile = io.read()

    # current input
    ccss = nwbfile.get_stimulus('ccss')
    current_stimulus = ccss.data[()]

    # current output
    ccs = nwbfile.get_acquisition('ccs')
    current_clamp = ccs.data[()]

    io.close()

    return nwbfile, ccss, ccs, current_stimulus, current_clamp


# read NWB file of cell and return necessary dap object
def NWBread(fpath, export_dir=str):

    nwbfile, ccss, ccs, current_stimulus, current_clamp = readNWBpatchClamp(fpath=fpath)

    # Populate a dictionary with file parameters for each of the fields
    ap = sl_sync_params_v2()

    ap['Dir'] = fpath
    ap['fname'] = nwbfile.identifier
    ap['comment'] = nwbfile.experiment_description
    ap['ExportDir'] = export_dir

    # ---- Condition/stimulus information
    ap['cond.gain'] = ccss.gain
    ap['cond.dc'] = ccss.description
    ap['cond.protocol'] = nwbfile.protocol

    # ---- Cell information
    ap['cell.date'] = nwbfile.session_start_time  # date that cell was recorded
    ap['cell.type'] = nwbfile.experiment_description.split()[-1]  # relevant description of cell type
    ap['cell.species'] = nwbfile.experiment_description.split()[0]    # species of the cell
    ap['cell.exp_condition'] = nwbfile.experiment_description.split()[1]    # experimental condition of the cell


    # ---- Special fields
    ap['sf.RMP_offset'] = int(nwbfile.notes.split()[-1])  # offset of the RMP


    return ap, nwbfile, current_clamp, current_stimulus

# perform response variability calculations
def respVariability(fpath, export_dir=''):


    #
    ap, nwbfile, current_clamp, current_stimulus = NWBread(fpath, export_dir)

    cell_id = ap['fname']

    # get voltage and current trace
    print('Loading data...')
    h = {}
    si = {}     # sampling intervals for each cell in us
    d = {}
    V = {}
    I = {}


    V[cell_id] = current_clamp

    # RMP offset correction
    num_pts_RMP = int(ap['rv']['prefix'] * 10000)
    nTrials = len(V[cell_id])
    if type(ap['sf.RMP_offset']) == int:  # perform RMP offset correction if specified
        V[cell_id] = V[cell_id] + ap['sf.RMP_offset']

    # compute and report resting membrane potential (RMP)
    x = 0
    for i in range(0, nTrials):
        trial_average = np.average(V[cell_id][i][:num_pts_RMP])
        x += trial_average
    cell_average = x / nTrials
    print('Average RMP %s: %d' % (cell_id, cell_average))


    I[cell_id] = current_stimulus


    # sampling interval for each sweep in us
    si[cell_id] = 1e6/10000 # sampling rate was 10000 Hz


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

    print('Cell #%s: Performing analysis...' % cell_id)
    sta = []; g = []; PHI = []; FREQ = []; SPbINS = []

    for i in range(0, nTrials):
        V_1 = V[cell_id][i]
        I_1 = I[cell_id][i]
        STA_, G_, phi_, freq_, spBins_ = doCalc(V_1, I_1, si, ap)
        sta.append(STA_); g.append(G_); PHI.append(phi_); FREQ.append(freq_); SPbINS.append(spBins_)

    STA[cell_id] = sta; G[cell_id] = g; phi[cell_id] = PHI; freq[cell_id] = FREQ; spBins[cell_id] = np.array(SPbINS)


    # compute the average spike rate across all trials
    avgSpikeRate = {}
    spikecount = 0
    for i in range(0, len(spBins[cell_id])):
        spikecount += sum(spBins[cell_id][i])

    avgSpikeRate[cell_id] = spikecount / (len(V[cell_id]) * ap['rv']['length'])
    print('Average Spike Rate: %d Hz' % avgSpikeRate[cell_id])

    # compute spike correlations
    sCorr = {}; sCorrVec = {}
    # If the average spike rate is less than an arbitrary amount send a warning and dont compute spike correlation
    if avgSpikeRate[cell_id] < ap['rv']['minSpikeRate']:
        print('Average spike rate too low to perform analyses')
    else:
        # Do spike correlations
        sCorr_, sCorrVec_ = spikeCorr(spBins[cell_id], fs/wPoints, ap['rv']['delta'])        # note that the pearson correlation values are slightly discordant between Matlab and python
        sCorr[cell_id] = sCorr_; sCorrVec[cell_id] = sCorrVec_

    # save results
    pickle.dump([sCorr, sCorrVec, V, I, spBins, ap, fs, wPoints, STA,
                 avgSpikeRate, G, phi, freq],
                open((export_dir + cell_id + '_results.pkl'), 'wb'))


    # ------- UPDATE .csv file of Response Variability Results
    df = pd.read_csv('/Volumes/HD1/White_noise/ResponseVariabilityCells.csv')

    df_append = df.append({'cell_id': cell_id,
                           'exp_condition': ap['cell.exp_condition'],
                           'cell_type': ap['cell.type'],
                           'gain': ap['cond.gain'],
                           'dc': ap['cond.dc'],
                           'firing rate': avgSpikeRate[cell_id],
                           'analysis_date': datetime.datetime.now().strftime("%I:%M%p %B %d, %Y")
                           }, ignore_index = True)


    df_append.to_csv('/Volumes/HD1/White_noise/ResponseVariabilityCells.csv')

##

def runCell(fpath, export_dir):
    respVariability(fpath, export_dir=export_dir)


