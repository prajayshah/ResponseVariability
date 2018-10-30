import os
from runCell_v2 import *
from rvPlotResults import *

def rvRunDir(dname, dosave = False):

    '''
    Runs all the files in a directory, and if multiple taus were used will plot the correlation as a function of tau

    :param dname:       Relative directory name
    :param dosave:      Saves the figures to the specified directory
    :return: sCorr      Correlation matrices
    :return: sCorrVec   List form of sCorr
    :return: tau        taus from different files
    :return: sigma      Sigmas of the different files
    '''

    #DATA_DIR = 'C:\Users\praja\OneDrive\UT MDPhD'

    #dpath = os.path.join(DATA_DIR, dname)
    dpath = dname # delete once stable, and replace dname in function definition
    d = []
    for i in os.listdir(dpath):         ## try to use glob.glob with this
        if i.endswith('.xlsx'):
            d.append(i)

    nFiles = len(d)

    if nFiles == 0:
        raise KeyError('No excel files found!')

    # analyse cell
    sCorr = {}; sCorrVec = {}; spBins = {}; R = {}
    for i in range(0, nFiles):
        print('Opening... %s' % d[i])

        sCorr_, sCorrVec_, spBins_ = runCell_v2(dpath, d[i][:-5], False, False)
        sCorr[i] = sCorr_; sCorrVec[i] = sCorrVec_; spBins[i] = spBins_


        print('finished analyzing %s .... ' % d[i])
        print('')
        print('')


    pFiles = []
    for i in os.listdir(dpath):  ## try to use glob.glob with this
        if i.endswith('.pkl'):
            pFiles.append(i)

    # plot results
    # for file in pFiles:
    #     rvPlotResults(dpath, file[:-4], dosave=True)


    print('')
    print('Done!')

    # if dosave == True:

