#!/usr/bin/env python

from os import *
from pyfits import *
from glob import glob

## -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

SkyFITSFileNameList = glob("*-bkg_region-sky.fits")

for SkyFITSFileName in SkyFITSFileNameList :
    DetFITSFileName = SkyFITSFileName.replace("sky", "det")

    print SkyFITSFileName, DetFITSFileName

    if not path.isfile(DetFITSFileName) :
        print "Error:", DetFITSFileName, "does not exits"
        continue

    SkyHDUList = open(SkyFITSFileName)
    SkyRegionList = SkyHDUList[1].data
    DetHDUList = open(DetFITSFileName, mode="update")
    DetRegionList = DetHDUList[1].data

    SkyR = SkyRegionList.field("R")
    DetR = DetRegionList.field("R")
    for i in xrange(len(SkyRegionList)) :  DetR[i] = SkyR[i]

    print SkyRegionList.field("R")
    print DetRegionList.field("R")

    SkyHDUList.close()
    DetHDUList.close()
