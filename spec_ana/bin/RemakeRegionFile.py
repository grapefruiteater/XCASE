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

RegionFile=argv[1]
mos1=commands.getoutput("head -1 log/prefix.log | gawk '{print $1}'")
mos2=commands.getoutput("head -2 log/prefix.log | tail -1 | gawk '{print $1}'")
pn=commands.getoutput("tail -1 log/prefix.log | gawk '{print $1}'")


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

print "#     Selcted Point Source\n"
for Prefix in (mos1, mos2, pn) :
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
        if skyY[0]=="+":
            skyY = degy+miny/60+secy/60/60.
        else :
            skyY = degy-miny/60-secy/60/60.
        output_esky2det = commands.getoutput("bin/esky2det.csh %f %f %s-obj-image-det.fits | gawk '{if(NF==2) print $0}'" % (skyX, skyY,Prefix))
        DetX = float(output_esky2det.split()[0])
        DetY = float(output_esky2det.split()[1])
        #print Prefix,DetX,DetY
        if Prefix==mos1:
            mos1DetX[i-n-1]=DetX
            mos1DetY[i-n-1]=DetY
        if Prefix==mos2:
            mos2DetX[i-n-1]=DetX
            mos2DetY[i-n-1]=DetY
        if Prefix==pn:
            pnDetX[i-n-1]=DetX
            pnDetY[i-n-1]=DetY
       
#astrometry between mos1,mos2,pn is slightly different.
# we here large margin for radius 10 arcsec
Rth=200.;
        
for Prefix in (mos1, mos2, pn) :
#for Prefix in (pn,) :
        SkyFITSFileName = "%s-bkg_region-sky.fits"%Prefix

        DetFITSFileName = SkyFITSFileName.replace("sky", "det")

        print SkyFITSFileName, DetFITSFileName
        if not os.path.isfile(DetFITSFileName) :
            print "Error:", DetFITSFileName, "does not exits"
            break;
        SkyHDUList = pyfits.open(SkyFITSFileName, mode="update")
        SkyRegionList = SkyHDUList[1].data
        DetHDUList = pyfits.open(DetFITSFileName, mode="update")
        DetRegionList = DetHDUList[1].data
        DetX=DetRegionList.field("DETX")   
        DetY=DetRegionList.field("DETY")   
        SkyR = SkyRegionList.field("R")
        DetR = DetRegionList.field("R")
        DetComp = DetRegionList.field("COMPONENT")

        if Prefix==mos1:
           for i in np.arange(DetComp.size) :
               if np.any(np.sqrt((DetX[i][0]-mos1DetX)**2+(DetY[i][0]-mos1DetY)**2)<Rth) == False :
                  DetR[i][0] = 0.
                  SkyR[i][0] = 0.

        if Prefix==mos2:
           for i in np.arange(DetComp.size) :
               if np.any(np.sqrt((DetX[i][0]-mos2DetX)**2+(DetY[i][0]-mos2DetY)**2)<Rth) == False :
                  DetR[i][0] = 0.
                  SkyR[i][0] = 0.

        if Prefix==pn:
           for i in np.arange(DetComp.size) :
               if np.any(np.sqrt((DetX[i][0]-pnDetX)**2+(DetY[i][0]-pnDetY)**2)<Rth) == False :
                  DetR[i][0] = 0.
                  SkyR[i][0] = 0.
            

        SkyHDUList.close()
        DetHDUList.close()

