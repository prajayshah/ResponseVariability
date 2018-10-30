from sl_sync_params import *
import pickle
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D

def rvPlotResults(Dir, fname, dosave = None):

    '''
    Plots the results from response variability analysis

    :param Dir:
    :param fname:
    :param dosave:

    '''

    ap = sl_sync_params()

    # define path to save plots

    # ----

    # example file:
    # fname = 'cell2_3_50_results'
    # Dir = 'C:\Users\praja\OneDrive\UT MDPhD\White noise\Human tissue\Jan 29, 2018\Cell 2\Gain 20'

    # load up the results files
    sCorr, sCorrVec, V, I, spBins, dap, fs, wPoints, STA, \
    avgSpikeRate, G, phi, freq = pickle.load(open(str(os.path.join(Dir, fname + '.pkl')), 'rb'))

    files = []
    for i in dap:
        files.append(dap[i]['fname'][0][:-1])

    # Plot results if the average spike rate is greater than the set amount

    for f in files:
        if avgSpikeRate[f] > ap['rv']['minSpikeRate']:
            nTrials, nPoints = V[f].shape
            sampling_rate = fs  # in Hz
            t = np.arange(0, nPoints) * (1.0 / sampling_rate)
            ind = np.where((t >= ap['rv']['prefix']) & (t <= (ap['rv']['prefix'] + ap['rv']['length'])))


            # plot traces of all trials - add titles, and make trace lines thinner
            plt.style.use('ggplot')
            fig, axes = plt.subplots(nTrials + 1, 1, sharex=True)
            axes[0].plot(t[ind], I[f][0][ind], color='gray')
            axes[0].set_ylabel('pA')
            axes[0].set_xlabel('seconds')
            for i in range(0, nTrials):
                axes[i+1].plot(t[ind], V[f][i][ind], color='black')
                axes[i+1].set_ylabel('mV')
            plt.xlim(min(t[ind]), max(t[ind]))
            plt.suptitle("Cell# %s" % f)
            fig.show()


            # plot correlation per cell
            plt.figure()
            sns.heatmap(sCorr[f], cmap='Greens', vmin=0, vmax=0.5)


            # plot STA
            plt.figure()
            sta_segments = pd.DataFrame(STA[f][0]).T
            for i in range(1, len(STA[f])):
                sta_segments = pd.concat([pd.DataFrame(STA[f][i]).T, sta_segments], 1)

            sta_segments.reset_index(drop=True)
            tSTA = np.linspace(0, ap['rv']['STAwindow'], sta_segments.shape[0])-ap['rv']['STAwindow']

            plt.style.use('seaborn-darkgrid')
            for column in sta_segments:
                plt.plot(pd.DataFrame(tSTA), sta_segments.iloc[:, column], marker='', color='grey', linewidth=1, alpha=0.1)

            staMean = np.mean(sta_segments.T)
            plt.plot(pd.DataFrame(tSTA), staMean, marker='', color='black', linewidth=3, alpha=1)
            plt.ylim(min(staMean), max(staMean))
            plt.title("Spike Triggered Average, n = % spikes" % sta_segments.shape[1], fontsize=12, fontweight=1)
            plt.xlabel("Time before peak of spike (ms)")
            plt.ylabel("Stimulus (pA)")


            # plot Raster
            plt.figure()
            nTrials, nPoints = spBins[f].shape
            T = np.array(range(0, nPoints-1))/(fs/wPoints)
            r, c = np.where(spBins[f] == 1)
            plt.style.use('ggplot')
            plt.plot(T[c], r+1, 'k.', ms = 5)   # note: r+1 is just to raise the plot off of the 0 line
            plt.ylim(0, nTrials+1)
            plt.ylabel('Trial #')
            plt.xlabel('Time (secs)')
            plt.title('Raster (avg freq = %s Hz)' % avgSpikeRate)


            # plot PCA
            # calculate PCA with as many components as trials in this cell
            pca = PCA(n_components=len(spBins[f])).fit(spBins[f])
            result = pd.DataFrame(pca.transform(spBins[f]), columns=['PCA%i' % i for i in range(len(spBins[f]))])
            # Plot initialisation
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter(result['PCA0'], result['PCA1'], result['PCA2'], cmap="Set2_r", s=60)
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
            ax.set_title("PCA %s" % f)

            fig.show()
        else:
            print("plots not made for cell %s" % f)

    if dosave == True:
        pass
    else:
        pass





