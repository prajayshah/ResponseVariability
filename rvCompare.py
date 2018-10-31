# combines and compares the results from multiple cells and produces plots

import os
import pickle
import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from mpl_toolkits.mplot3d import Axes3D
from spikeCorr import *
from sl_sync_params_v2 import *
from generate_random_color import *

def rvCompare(files = [], cells_of_interest = []):

    '''

    :param: files       List of directory paths which contain file results of cells to compare
    :return:            Plots

    '''

    ap_2 = sl_sync_params_v2()

    # specify and find files of interest - only for developing code
    files = []; cells_of_interest = ['18o22014', '18220020']
    files.append('/Volumes/HD1/White_noise/Human_tissue/Epilepsy cases/results/')  # insert full paths to cells
    files.append('/Volumes/HD1/White_noise/Human_tissue/Tumor cases/results/')

    ##
    fpaths = []
    for path in files:
        for i in os.listdir(path):  ## try to use glob.glob with this
            if i.endswith('.pkl'):
                fpaths.append(os.path.join(path, i))

    # load up the results files for the first recording (to initialize the necessary results which are in dictionary formats)
    initial = fpaths[0]
    print('Loading...%s' % initial)
    sCorr, sCorrVec, V, I, spBins, dap, fs, wPoints, STA, avgSpikeRate, G, phi, freq = pickle.load(open(initial, 'rb'))
    ap = {dap['fname']: dap}

    cell_name = dap["fname"]
    nTrials = spBins[cell_name].shape[0]
    cell2trials = [cell_name] * nTrials
    nPoints = spBins[cell_name].shape[1]

    # use a loop to gather results from remaining recordings, and add to spBins dictionary as an additional result under the appropriate recording name
    for i in fpaths[1:]:

        # load up the results files
        print('Loading...%s' % i)
        sCorr_, sCorrVec_, V_, I_, spBins_, dap_, fs_, wPoints_, STA_, \
        avgSpikeRate_, G_, phi_, freq_ = pickle.load(open(i, 'rb'))
        for j in spBins_.keys():
            spBins[j] = spBins_[j]
            cell2trials.extend([j] * spBins_[j].shape[0])
            nTrials += spBins_[j].shape[0]
            V[j] = V_[j]
            I[j] = I_[j]
            ap[j] = dap_
            freq[j] = freq_[j]
            avgSpikeRate[j] = avgSpikeRate_[j]


    # cells of particular interest

    cells_to_analyze = cells_of_interest
    if len(cells_to_analyze) == 0:
        cells_to_analyze = [*spBins]

    ## organize trials into a single array for analysis
    allspBins = []
    for i in cells_to_analyze:
        for j in range(len(spBins[i])):
            allspBins.append(spBins[i][j])

    # ==================================================================================================================
    # Choosing random colors to assign to all cells - to be used uniquely for each cell throughout plotting
    # ==================================================================================================================


    # generate_random_colors and assign to individual recordings
    colors = []
    for i in range(0, len(cells_to_analyze)):
        colors.append(generate_new_color(colors, pastel_factor=0.2))

    color_cells = {}
    for i in range(0, len(cells_to_analyze)):
        color_cells[cells_to_analyze[i]] = colors[i]


    ##
    ## PLOTTING
    ##
    # ==================================================================================================================
    # Plot 1: Raster of all cells
    print('Plotting Raster')
    # ==================================================================================================================
    plt.style.use('seaborn-white')
    fig, ax = plt.subplots(len(cells_to_analyze), 1, sharex=True, squeeze=False)
    T = np.array(range(0, nPoints - 1)) / (fs / wPoints)

    fig.subplots_adjust(hspace=0, wspace=0)

    for cell in cells_to_analyze:
        r = np.where(spBins[cell] == 1)[0] #+ trials
        c = np.where(spBins[cell] == 1)[1]
        ax[cells_to_analyze.index(cell), 0].plot(T[c], r + 1, 'k.', ms=5, color=color_cells[cell])


    for i in range(len(ax)):

        # # set color scheme of yaxis labels based on cell sample
        # if cells_to_analyze[i] in L2_3:
        #     ax[i].set_ylabel('%d: %dHz (L2/3)' % (i + 1, avgSpikeRate[cells_to_analyze[i]]), rotation=0)
        # elif cells_to_analyze[i] in HuL5:
        #     ax[i].set_ylabel('%d: %dHz (L5)' % (i + 1, avgSpikeRate[cells_to_analyze[i]]), rotation=0)
        # elif cells_to_analyze[i] in MsL5:
        #     ax[i].set_ylabel('%d: %dHz (Ms)' % (i + 1, avgSpikeRate[cells_to_analyze[i]]), rotation=0)
        # if cells_to_analyze[i] in human_cells:
        #     ax[i].yaxis.label.set_color('red')


        #ax[i].yaxis.set_label_coords(-0.05, 0.40)

        trials = len(spBins[cells_to_analyze[i]])
        ax[i,0].set_ylim(0, trials + 1)
        ax[i,0].set_yticklabels([])
        ax[i,0].spines['left'].set_visible(False)
        ax[i,0].spines['right'].set_visible(False)
        #ax[i,0].spines['top'].set_visible(False)
    ax[len(cells_to_analyze)-1,0].set_xlabel('Time (secs)')
    ax[len(cells_to_analyze)-1,0].set_xlim(T.min(), T.max())
    fig.suptitle('Raster of %d cells' % len(cells_to_analyze))
    fig.show()

    # ==================================================================================================================
    # Plot 2: PCA of all cells
    print('Plotting PCA')
    # ==================================================================================================================

    # PCA calculation on all trials together
    pca = PCA(n_components=10).fit(allspBins)
    result = pd.DataFrame(pca.transform(allspBins), columns=['PCA%i' % i for i in range(10)])

    # Plot PCA
    plt.style.use('seaborn-white')
    fig = plt.figure()
    ax = Axes3D(fig)
    trials = 0  # ticker that counts trials from a cell to help make sure the next cell goes on the trials after the preceeding cell
    for cell in cells_to_analyze:
        if cell in []:
            ntrial = spBins[cell].shape[0]
            ax.scatter(result.iloc[trials:trials + ntrial]['PCA0'], result.iloc[trials:trials + ntrial]['PCA1'],
                       result.iloc[trials:trials + ntrial]['PCA2'],
                       cmap="Set2_r", s=60, marker='v',
                       color=color_cells[cell])
            # draw border around human cells (in order to allow comparison with mouse cells) - disabled
            # if cell in human_cells:
            #     ntrial = spBins[cell].shape[0]
            #     ax.scatter(result.iloc[trials:trials + ntrial]['PCA0'], result.iloc[trials:trials + ntrial]['PCA1'],
            #                result.iloc[trials:trials + ntrial]['PCA2'],
            #                cmap="Set2_r", s=60, marker='v', facecolors='none', edgecolors='black', linewidths=1.5,
            #                color=color_cells[cell])


        else:
            ntrial = spBins[cell].shape[0]
            ax.scatter(result.iloc[trials:trials+ntrial]['PCA0'], result.iloc[trials:trials+ntrial]['PCA1'], result.iloc[trials:trials+ntrial]['PCA2'],
                       cmap="Set2_r", s=60,
                       color=color_cells[cell])
            # draw border around human cells (in order to allow comparison with mouse cells) - disabled for now
            # if cell in human_cells:
            #     ntrial = spBins[cell].shape[0]
            #     ax.scatter(result.iloc[trials:trials + ntrial]['PCA0'], result.iloc[trials:trials + ntrial]['PCA1'],
            #                result.iloc[trials:trials + ntrial]['PCA2'],
            #                cmap="Set2_r", s=60, facecolors='none', edgecolors='black', linewidths=1.5,
            #                color=color_cells[cell])
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

    # make the axis graph panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    ax.set_title("PCA (%d cells)" % len(cells_to_analyze))

    fig.show()

    # ==================================================================================================================
    # Plot 3: Pearson correlation of all cells
    print('Plotting correlation')
    # ==================================================================================================================

    ## calculate spike correlations across all trials in comparison - NOTE: allspBins was set earlier
    sCorr, sCorrVec = spikeCorr(allspBins, fs/wPoints, ap_2['rv']['delta'])


    ## determine row colors to use
    row_colors = []
    cell2trials = []
    for i in cells_to_analyze:
        cell2trials.extend([i] * spBins[i].shape[0])
    for i in cell2trials:
        row_colors.append(color_cells[i])

    # make row_colors based on human/mouse
    colors2 = []
    for i in range(3):
        colors2.append(generate_new_color(colors2, pastel_factor=0.2))

    # color_cells2 = {}
    # for i, cell in enumerate(all_cells):
    #     if cell in HuL5:
    #         color_cells2[cell] = colors2[0]
    #     elif cell in L2_3:
    #         color_cells2[cell] = colors2[1]
    #     elif cell in MsL5:
    #         color_cells2[cell] = colors2[2]

    # row_colors2 = []
    # for i in cell2trials:
    #     row_colors2.append(color_cells2[i])

    ## plot correlations of all trials
    sCorr[np.isnan(sCorr)] = 0 # needed to avoid Linkage 'Z' contains negative counts error - that seems to be due to either NaN values in the dataframe or very long floats that are close to zero
    sCorr_t = sCorr + sCorr.T - np.diag(np.diag(sCorr))

    # # with heirachical clustering
    # sns.clustermap(sCorr_t, cmap='Greens', row_cluster=True, col_cluster=True, vmin=0, vmax=1, row_colors=[row_colors, row_colors2])
    # plt.show()
    # # w/o heirarchical clustering
    sns.clustermap(sCorr_t, cmap='Greens', row_cluster=False, col_cluster=False, vmin=0, vmax=1, row_colors=[row_colors])
    plt.show()



    # ==================================================================================================================
    # Plot 4: tSNE of all cells
    print('Plotting tSNE')
    # ==================================================================================================================

    # PCA calculation on all trials together
    tsne = TSNE(n_components=3, perplexity=30).fit_transform(allspBins)
    result = pd.DataFrame(tsne, columns=['tSNE%i' % i for i in range(3)])

    # Plot tSNE
    plt.style.use('seaborn-white')
    fig = plt.figure()
    ax = Axes3D(fig)
    trials = 0  # ticker that counts trials from a cell to help make sure the next cell goes on the trials after the preceeding cell
    for cell in cells_to_analyze:
        if cell in []:
            ntrial = spBins[cell].shape[0]
            ax.scatter(result.iloc[trials:trials + ntrial]['tSNE0'], result.iloc[trials:trials + ntrial]['tSNE1'],
                       result.iloc[trials:trials + ntrial]['tSNE2'],
                       cmap="Set2_r", s=60, marker='v',
                       color=color_cells[cell])
            # # draw border around human cells (in order to allow comparison with mouse cells)
            # if cell in human_cells:
            #     ntrial = spBins[cell].shape[0]
            #     ax.scatter(result.iloc[trials:trials + ntrial]['tSNE0'], result.iloc[trials:trials + ntrial]['tSNE1'],
            #                result.iloc[trials:trials + ntrial]['tSNE2'],
            #                cmap="Set2_r", s=60, marker='v', facecolors='none', edgecolors='black', linewidths=1.5,
            #                color=color_cells[cell])


        else:
            ntrial = spBins[cell].shape[0]
            ax.scatter(result.iloc[trials:trials + ntrial]['tSNE0'], result.iloc[trials:trials + ntrial]['tSNE1'],
                       result.iloc[trials:trials + ntrial]['tSNE2'],
                       cmap="Set2_r", s=60,
                       color=color_cells[cell])
            # draw border around human cells (in order to allow comparison with mouse cells)
            # if cell in human_cells:
            #     ntrial = spBins[cell].shape[0]
            #     ax.scatter(result.iloc[trials:trials + ntrial]['tSNE0'], result.iloc[trials:trials + ntrial]['tSNE1'],
            #                result.iloc[trials:trials + ntrial]['tSNE2'],
            #                cmap="Set2_r", s=60, facecolors='none', edgecolors='black', linewidths=1.5,
            #                color=color_cells[cell])
        trials += ntrial

    # make simple, bare axis lines through space:
    xAxisLine = ((min(result['tSNE0']), max(result['tSNE0'])), (0, 0), (0, 0))
    ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r')
    yAxisLine = ((0, 0), (min(result['tSNE1']), max(result['tSNE1'])), (0, 0))
    ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r')
    zAxisLine = ((0, 0), (0, 0), (min(result['tSNE2']), max(result['tSNE2'])))
    ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r')
    # label the axes
    ax.set_xlabel("tSNE 1")
    ax.set_ylabel("tSNE 2")
    ax.set_zlabel("tSNE 3")

    # make the axis graph panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    ax.set_title("tSNE (%d cells)" % len(cells_to_analyze))

    fig.show()

    # ==================================================================================================================
    print('Done!')
    # ==================================================================================================================















