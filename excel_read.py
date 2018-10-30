import os
from sl_sync_params import *
import pandas as pd

def excel_read(dir, fname):
    '''
    Read in an excel spreadsheet and parse the data out of it.

    :param dir:         location of excel spreadsheet
    :param fname:       name of the excel spreadsheet
    :return:
    '''

    def_params = sl_sync_params()
    sf_names = def_params['sf']['names']

    full_path = os.path.join(dir, fname + '.xlsx')

    ss = pd.read_excel(full_path, sheet_name='Sheet1')

    # read in SHEET 2 - channel assignments
    ch = pd.read_excel(full_path, sheet_name='Sheet2')

    if ss.count().sum() < 3:
        # There must be more than just a Dir and (un-named) filename field
        raise ValueError('Not enough fields.')

    # get number of conditions and files - what is number of conditions? which conditions?
    nfiles, a = ss.shape
    ncond = (a-6)//3

    # get the condition names - needs work!
    cond_names = []
    for i in range(0, ncond):
        # pass
        cond_names.append(ss.columns[2])

    # find any special fields
    sf_index = []
    for i in sf_names:
        sf_index.append(ss.columns.tolist().index(i))


    # Populate a dictionary with file parameters for each of the fields

    ap = {}
    for i in range(0, nfiles):      ## be careful with this one - it might not retrieve the correct row with >1 nfiles
        ap[i] = sl_sync_params()
        ap[i]['Dir'] = ss.iloc[:, 0].values.tolist()
        ap[i]['fname'] = ss.iloc[:, 1].values.tolist()
        ap[i]['cond.names'] = []

        c = 0
        for j in range(0, ncond):
            ap[i]['cond.times'].append(ss.iloc[0,2+j])
            ap[i]['cond.names'].append(cond_names)
            ap[i]['cond.fname'].append(ss.iloc[0,1])

        # get the special fields
        # c = 0
        # ap[i]['sf']['names'] = []
        # for j in range(0, len(sf_names)):
        #     for x in sf_index:
        #         ap[i]['sf']['names'].append(sf_names[j])
        #         ap[i]['sf']['vals'].append(ss.values[0][x])

        # get the special fields
        c = 0
        #ap[i]['sf']['names'] = [sf_names]
        # for name in sf_names:
        #     ap[i]['sf'][name] = ss.values[0][sf_names.index(name)]
        for x in sf_index:
            ap[i]['sf'][sf_names[x-sf_index[0]]] = ss.values[0][x]
            #ap[i]['sf']['names'].append(sf_names)
            #ap[i]['sf']['vals'].append(ss.values[0][x])

    # Set the channel to analyze
    ap[i]['ch'] = ss['AnChan'].values[0]

    # Get the channel labels from the table (sheet 2). If there are multiple files
    # for the same slice, it is assumed that the channel assignments did not change
    # during the experiment

    chfnames = []
    for j in range(0, len(ch)):
        chfnames.append(ch.iloc[0, 0])

    ## reindex stuff that still needs to be coded in

    ap[i]['chlabels'] = ch.iloc[0:, 1:]

    ## rename ap with appropriate keyname
    fname = ap[i]['fname'][0][:-1]
    ap[fname] = ap.pop(i)

    return ss, ap






