# combines and compares the results from multiple cells and produces plots
# TODO test rvCompare on multiple cells

import os
import pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
from spikeCorr import *
from sl_sync_params import *
from rand_cmap import *

def rvCompare(files):

    '''

    :param files:       List of directory paths which contain file results of cells to compare
    :return:            Plots

    '''

    ap = sl_sync_params()

    # specify and find files of interest - only for developing code
    # files = []
    # files.append('C:\Users\praja\OneDrive\UT MDPhD\White noise\Human tissue\Jan 25, 2018\Cell 3\Gain 20')  # insert full paths to cells

    ##
    fpaths = []
    for path in files:
        for i in os.listdir(path):  ## try to use glob.glob with this
            if i.endswith('.pkl'):
                fpaths.append(os.path.join(path, i))

    # load up the results files for the first cell (to initialize the necessary results which are in dictionary formats)
    sCorr, sCorrVec, V, I, spBins, dap, fs, wPoints, STA, \
    avgSpikeRate, G, phi, freq = pickle.load(open(fpaths[0], 'rb'))

    cell_name = spBins.keys()[0]
    nTrials = spBins[cell_name].shape[0]
    cell2trials = [cell_name] * nTrials

    # use a loop to gather results from remaining cells, and add to spBins dictionary as an additional result under the appropriate cell name
    for i in fpaths[1:]:

        # load up the results files
        sCorr_, sCorrVec_, V_, I_, spBins_, dap_, fs_, wPoints_, STA_, \
        avgSpikeRate_, G_, phi_, freq_ = pickle.load(open(i, 'rb'))
        for j in spBins_.keys():
            spBins[j] = spBins_[j]
            cell2trials.extend([j] * spBins_[j].shape[0])
            nTrials += spBins_[j].shape[0]
            nPoints = spBins_[j].shape[1]
            V[j] = V_[j]
            I[j] = I_[j]


    cells = spBins.keys(); colors = dict(zip(cells, ['orange', 'brown']))       # specify which colors to use for plotting individual cells
    # TODO add random selection of colors from a pre-specified pallate

    # random color map generation - I think this is a good one to use from here on out
    color = rand_cmap(len(cells), type='bright', first_color_black=False, verbose=False)
    color_cells = {}
    for i in range(len(cells)):
        color_cells[cells[i]] = color(i)


    # for cell in spBins.keys():
    #     colors[cell] = [colormap(i) for i in np.linspace(0, 0.9, len(spBins.keys()))][cells.index(cell)]
    row_colors = []
    for i in cell2trials:
        row_colors.append(color_cells[i])


    ## plot Raster of all cells
    #plt.style.use('ggplot')
    f, ax = plt.subplots(len(cells), 1, sharex=True)
    T = np.array(range(0, nPoints - 1)) / (fs / wPoints)
    # gs1 = gridspec.GridSpec(4, 4)
    # gs1.update(wspace=0.025, hspace=0.05)  # set the spacing between axes.
    #

    f.subplots_adjust(hspace=0, wspace=0)

    for cell in cells:
        r = np.where(spBins[cell] == 1)[0] #+ trials
        c = np.where(spBins[cell] == 1)[1]
        ax[cells.index(cell)].plot(T[c], r + 1, 'k.', ms=5, color=color_cells[cell])
        #ax[cells.index(cell)].plot(T[c], r + 1, 'k.', ms=5, color=color(cells.index(cell)))

        # previous attempts at controlling color of cell raster - can delete as long as current one is stable across MANY more cells
        # ax[cells.index(cell)].plot(T[c], r + 1, 'k.', ms=5)
        # ax[cells.index(cell)].plot(T[c], r + 1, 'k.', ms=5, color=colors[cell])

    for i in ax:
        i.set_ylabel('Trial #')
        trials = len(spBins[cells[ax.tolist().index(i)]])
        i.set_ylim([0, trials + 1])
    ax[-1].set_xlabel('Time (secs)')
    f.suptitle('Raster of %d cells' % len(cells))
    plt.show()

    ## organize trials into a single array for subsequent PCA analysis
    allspBins = []
    for i in cells:
        for j in range(len(spBins[i])):
            allspBins.append(spBins[i][j])

    # PCA calculation on all trials together
    pca = PCA(n_components=len(allspBins)).fit(allspBins)
    result = pd.DataFrame(pca.transform(allspBins), columns=['PCA%i' % i for i in range(len(allspBins))])

    # Plot PCA
    fig = plt.figure()
    ax = Axes3D(fig)
    trials = 0  # ticker that counts trials from a cell to help make sure the next cell goes on the trials after the preceeding cell
    for cell in cells:
        ntrial = spBins[cell].shape[0]
        ax.scatter(result.iloc[trials:trials+ntrial]['PCA0'], result.iloc[trials:trials+ntrial]['PCA1'], result.iloc[trials:trials+ntrial]['PCA2'], cmap="Set2_r", s=60,
                   color=color_cells[cell])
        trials += ntrial
    # make simple, bare axis lines through space:
    xAxisLine = ((min(result['PCA0']), max(result['PCA0'])), (0, 0), (0, 0))
    ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r')
    yAxisLine = ((0, 0), (min(result['PCA1']), max(result['PCA1'])), (0, 0))
    ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r')
    zAxisLine = ((0, 0), (0, 0), (min(result['PCA2']), max(result['PCA2'])))
    ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r')
    # label the axes
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    ax.set_title("PCA (%d cells)" % len(cells))

    plt.show()


    ## calculate spike correlations across all trials in comparison
    sCorr, sCorrVec = spikeCorr(allspBins, fs/wPoints, ap['rv']['delta'])

    ## plot correlations of all trials
    sns.clustermap(sCorr, cmap='Greens', row_cluster=False, col_cluster= False, vmin=0, vmax=0.5, row_colors=row_colors)















