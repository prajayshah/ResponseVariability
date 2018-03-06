# combines the results from multiple cells

import os
import pickle

def rvCompare():

    # specify and find files of interest
    dFile = []
    dFile.append('C:\Users\praja\OneDrive\UT MDPhD\White noise\Human tissue\Jan 25, 2018\Cell 3\Gain 20')  # insert full paths to cells

    fpaths = []
    for path in dFile:
        for i in os.listdir(path):  ## try to use glob.glob with this
            if i.endswith('.pkl'):
                fpaths.append(os.path.join(path, i))

    # load up the results files
    sCorr, sCorrVec, V, I, spBins, dap, fs, wPoints, STA, \
    avgSpikeRate, G, phi, freq = pickle.load(open(fpaths[0], 'rb'))

    for i in fpaths[1:]:

        # load up the results files
        sCorr_, sCorrVec_, V_, I_, spBins_, dap_, fs_, wPoints_, STA_, \
        avgSpikeRate_, G_, phi_, freq_ = pickle.load(open(i, 'rb'))
        for j in spBins_.keys():
            spBins[j] = spBins_[j]


    allspBins = []
    for i in spBins:
        allspBins.append(spBins[i])












