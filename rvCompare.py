# combines the results from multiple cells

import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
from spikeCorr import *
from sl_sync_params import *

def rvCompare():

    ap = sl_sync_params()

    # specify and find files of interest
    dFile = []
    dFile.append('C:\Users\praja\OneDrive\UT MDPhD\White noise\Human tissue\Jan 25, 2018\Cell 3\Gain 20')  # insert full paths to cells

    fpaths = []
    for path in dFile:
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


    cells = spBins.keys(); colors = dict(zip(cells, ['orange', 'yellow', 'brown']))
    # for cell in spBins.keys():
    #     colors[cell] = [colormap(i) for i in np.linspace(0, 0.9, len(spBins.keys()))][cells.index(cell)]
    row_colors = []
    for i in cell2trials:
        row_colors.append(colors[i])



    # plot Raster
    plt.figure()
    plt.style.use('ggplot')
    T = np.array(range(0, nPoints - 1)) / (fs / wPoints)
    trials = 0   # ticker that counts trials from a cell to help make sure the next cell goes on the trials after the preceeding cell
    for cell in spBins.keys():
        r = np.where(spBins[cell] == 1)[0] + trials
        c = np.where(spBins[cell] == 1)[1]
        plt.plot(T[c], r + 1, 'k.', ms=5, color = colors[cell])
        trials += spBins[cell].shape[0]     # update ticker with number of trials from the current cell
    plt.ylim(0, nTrials + 1)
    plt.ylabel('Trials')
    plt.xlabel('Time (secs)')
    plt.title('Raster of %d cells' % len(spBins.keys()))
    plt.show()

    # organize trials into a single array
    allspBins = []
    for i in cells:
        for j in range(len(spBins[i])):
            allspBins.append(spBins[i][j])

    # PCA calculation on all trials together
    pca = PCA(n_components=10).fit(allspBins)
    result = pd.DataFrame(pca.transform(allspBins), columns=['PCA%i' % i for i in range(10)])

    # Plot PCA
    fig = plt.figure()
    ax = Axes3D(fig)
    trials = 0  # ticker that counts trials from a cell to help make sure the next cell goes on the trials after the preceeding cell
    for cell in cells:
        ntrial = spBins[cell].shape[0]
        ax.scatter(result.iloc[trials:trials+ntrial]['PCA0'], result.iloc[trials:trials+ntrial]['PCA1'], result.iloc[trials:trials+ntrial]['PCA2'], cmap="Set2_r", s=60,
                   color=colors[cell])
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

    # spike correlations across all trials in comparison
    sCorr, sCorrVec = spikeCorr(allspBins, fs/wPoints, ap['rv']['delta'])

    # plot correlations of all trials
    plt.figure()
    sns.clustermap(sCorr, cmap='Greens', row_cluster=True, col_cluster= False, vmin=0, vmax=0.5, row_colors=row_colors)
    # add title and cell labels on axes














