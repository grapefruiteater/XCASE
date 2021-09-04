#!/usr/bin/env python

from sys import *
import os
import pyfits
from glob import glob
import commands
import numpy as np

if len(argv) != 2 : 
   print "INPUT ERROR in python";
   print "python DeleteRegionPS.py Region.reg"
   print " Region.reg    -  regionfile"
   print ""
   print " output "
   print "   mos1S001-bkg_region-det.fits(update)  "
   print "   mos1S001-bkg_region-sky.fits(update)  "
   print "   mos2S002-bkg_region-det.fits(update)  "
   print "   mos2S002-bkg_region-sky.fits(update)  "
   print "   pnS003-bkg_region-det.fits(update)  "
   print "   pnS003-bkg_region-sky.fits(update)  "
   print ""
   print " e.g : python DeleteRegionPS.py Region.reg"
   exit();

RegionFile = argv[1]

## -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

RegionList = open(RegionFile,"r")
PSList = RegionList.readlines()
mark = PSList.index("fk5\n")

mos1DetX=np.arange(len(PSList)-mark-1,dtype=np.float)
mos1DetY=np.arange(len(PSList)-mark-1,dtype=np.float)
mos2DetX=np.arange(len(PSList)-mark-1,dtype=np.float)
mos2DetY=np.arange(len(PSList)-mark-1,dtype=np.float)
pnDetX=np.arange(len(PSList)-mark-1,dtype=np.float)
pnDetY=np.arange(len(PSList)-mark-1,dtype=np.float)

for Prefix in ("mos1S001", "mos2S002", "pnS003") :
    print "#     Selcted Point Source\n"
    print Prefix
    PSList=PSList
    n=mark
    for i in xrange(n+1,len(PSList)):
        skyX = PSList[i].split(",")[0]
        skyX = skyX.split("(")[1]
        degx = float(skyX.split(":")[0])
        minx = float(skyX.split(":")[1])
        secx = float(skyX.split(":")[2])
        skyY = PSList[i].split(",")[1]
        degy = float(skyY.split(":")[0])
        miny = float(skyY.split(":")[1])
        secy = float(skyY.split(":")[2]) 
        skyX = degx*15.+minx*15./60+secx*15./60/60
        skyY = degy-miny/60-secy/60/60
        print skyX,skyY
        output_esky2det = commands.getoutput("bin/esky2det.csh %f %f %s-obj-image-sky.fits | head -1" % (skyX, skyY,Prefix))
        DetX = float(output_esky2det.split()[0])
        DetY = float(output_esky2det.split()[1])
        print DetX,DetY
        if Prefix=="mos1S001":
            mos1DetX[i-n-1]=DetX
            mos1DetY[i-n-1]=DetY
        if Prefix=="mos2S002":
            mos2DetX[i-n-1]=DetX
            mos2DetY[i-n-1]=DetY
        if Prefix=="pnS003":
            pnDetX[i-n-1]=DetX
            pnDetY[i-n-1]=DetY

for Prefix in ("mos1S001", "mos2S002", "pnS003") :

        SkyFITSFileName = "%s-bkg_region-sky.fits"%Prefix

        DetFITSFileName = SkyFITSFileName.replace("sky", "det")

        print SkyFITSFileName, DetFITSFileName

        if not os.path.isfile(DetFITSFileName) :
            print "Error:", DetFITSFileName, "does not exits"
            continue

        SkyHDUList = pyfits.open(SkyFITSFileName, mode="update")
        SkyRegionList = SkyHDUList[1].data
        DetHDUList = pyfits.open(DetFITSFileName, mode="update")
        DetRegionList = DetHDUList[1].data
        DetX=DetRegionList.field("DETX")
        DetY=DetRegionList.field("DETY")        
        SkyR = SkyRegionList.field("R")
        DetR = DetRegionList.field("R")
        print "#     R Threshold\n"
        for i in xrange(len(SkyR)):
            if 0<SkyR[i][0]<=15.*20:
                print SkyR[i][0]
                SkyR[i][0]=15.*20
                DetR[i][0]=15.*20
        print "#     Mask Delete\n"
        if Prefix=="mos1S001":
            MinList=[]
            Mos1DetX=np.arange(len(DetX),dtype=np.float)
            Mos1DetY=np.arange(len(DetX),dtype=np.float)
            for i in xrange(len(DetX)):
                Mos1DetX[i]=DetX[i][0]
                Mos1DetY[i]=DetY[i][0]
            for i in xrange(len(mos1DetX)):
                Sub=abs((Mos1DetX-mos1DetX[i])**2+(Mos1DetY-mos1DetY[i])**2)
                Min=np.argmin(Sub)
                print Min,mos1DetX[i],Mos1DetX[Min]
                MinList.append(Min)
            MinList=sorted(MinList, reverse=True)
            for i in xrange(len(MinList)):
                print MinList[i]  
                DetR[MinList[i]] = 0.
                SkyR[MinList[i]] = 0.
           
        if Prefix=="mos2S002":
            MinList=[]
            Mos2DetX=np.arange(len(DetX),dtype=np.float)
            Mos2DetY=np.arange(len(DetX),dtype=np.float)
            for i in xrange(len(DetX)):
                Mos2DetX[i]=DetX[i][0]
                Mos2DetY[i]=DetY[i][0]
            for i in xrange(len(mos2DetX)):
                Sub=abs((Mos2DetX-mos2DetX[i])**2+(Mos2DetY-mos2DetY[i])**2)
                Min=np.argmin(Sub)
                print Min,mos2DetX[i],Mos2DetX[Min]
                MinList.append(Min)
            MinList=sorted(MinList, reverse=True)
            for i in xrange(len(MinList)):
                print MinList[i]
                DetR[MinList[i]] = 0.
                SkyR[MinList[i]] = 0.
                
        if Prefix=="pnS003":
            MinList=[]
            PNDetX=np.arange(len(DetX),dtype=np.float)
            PNDetY=np.arange(len(DetX),dtype=np.float)
            for i in xrange(len(DetX)):
                PNDetX[i]=DetX[i][0]
                PNDetY[i]=DetY[i][0]
            for i in xrange(len(pnDetX)):
                Sub=abs((PNDetX-pnDetX[i])**2+(PNDetY-pnDetY[i])**2)
                Min=np.argmin(Sub)
                print Min,pnDetX[i],PNDetX[Min]
                MinList.append(Min)
            MinList=sorted(MinList, reverse=True)
            for i in xrange(len(MinList)):
                DetR[MinList[i]] = 0.
                SkyR[MinList[i]] = 0.

        SkyHDUList.close()
        DetHDUList.close()




