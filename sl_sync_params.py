# The ap dictionary contains analysis parameters for analysis
# use by setting analysis parameters here and then calling this function

import numpy as np

def sl_sync_params():
    ap = {}

    #######

    # ---- Directory information
    ap['Dir'] = ''
    ap['ExportDir'] = ''
    ap['fname'] = ''
    ap['comment'] = ''

    # ---- Condition information
    ap['cond.fname'] = []
    ap['cond.names'] = []
    ap['cond.times'] = []
    ap['chlabels'] = []

    # ---- Channel information
    ap['ch'] = 1        # Channel to analyse


    # ---- Special fields
    ap['sf'] = {'names': ['Comment', 'ExportDir', 'Tags', 'GroupInc'],
                'vals': []}


    # Response Variability
    ap['rv'] = {'window': 1,            # Time window for binning spikes
                'prefix': 0.0467,       # length of time in seconds before the noise starts
                'length': 2.5,          # length of noise sequence in seconds
                'currentLength': 0.3,   # length of time over which to compute DC current level
                'threshold': 20,        # Threshold for finding spikes
                'spikeWidth': 2,        # width of the spike
                'delta': 3,             # delta as per Galan for performing spike train correlations
                'minSpikeRate': 6,      # Minimum spike rate in Hz. If lower than this abort analysis.
                'STAwindow': 30,        # Time window to keep for STA
                'frange': 1000,         # Range of frequencies to compute frequency dependent gain
                'computeGain': False,    # Toggle to computation to save time
                'f': np.array([0, 1, 1, 1, 1, 0]),  # Filter to convolve traces to remove spurious peaks
                'dimRed': 'pca',        # Dimensionality reduction techniques 'pca' or 'tsne'
                'tsne.Perplexity': 32}  # tsne perplexity
    #

    #######

    return ap

