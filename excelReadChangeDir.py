from excel_read import *

def excelReadChangeDir(Dir, filename):
    '''
    Reads in an excel data sheet and changes the ap directory to the directory where the excel spreadsheet is.
    By definition, the excel spreadsheet must be in the same directory as the data file.

    :param Dir:         directory where data file is
    :param filename:    name of data file (without extensions) - must be a string
    :return:            an updated ap dictionary
    '''

    ss, ap = excel_read(Dir, filename)

    for i in ap.keys():
        ap[i]['Dir'] = Dir

    return ap
