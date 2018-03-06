from excelReadChangeDir import *
from respVariability import *

def runCell_v2(Dir, excelName, doplot = True, dosave = False):

    dap = excelReadChangeDir(Dir, excelName)
    sCorr, sCorrVec, spBins = respVariability(dap, excelName)


    return sCorr, sCorrVec, spBins #, R

