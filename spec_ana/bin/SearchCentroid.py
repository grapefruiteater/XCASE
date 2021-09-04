#!/usr/bin/env python
from sys import *
import os
import numpy as np
import commands
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM

kpctocm=3.085677581e21
Mpctocm=3.085677581e24
Mp=1.6726231e-24
Msun=1.989e33
Munit=1e14

#--------------
# Z=0.2
mu=0.5964
ratio=1./1.174417  # ratio=n_H/n_e
ratio_HE=0.0872083 # ratio_HE=n_HE/n_H
rho_factor=1.92574264507411  #rho_g= rho_factor mu M_p n_e
n_factor=1+ratio*(1+ratio_HE) # n=n_factor*ne
#---------------------------------
H0=70.
Om0=0.3;
OL=0.7;
cos=FlatLambdaCDM(H0, Om0, Neff=0)
#-------------------------------------------------
mos1=commands.getoutput("head -1 log/prefix.log | gawk '{print $1}'")
mos2=commands.getoutput("head -2 log/prefix.log | tail -1 | gawk '{print $1}'")
pn=commands.getoutput("tail -1 log/prefix.log  | gawk '{print $1}'")
#-------------------------------------------------
if len(argv) != 6 :
    print ""
    print "Usage:", argv[0], "rate_image centX centY z pointsource_reg.file"
    print " rate_image            - all epic rate image "
    print " centX                 - initial value X [pixel] "
    print " centY                 - initial value Y [pixel] "
    print " z                     - redshift "
    print " pointsource_reg.file  - region file of point source"
    print "";
    print "Output"
    print "  region files to extract spectrum  "
    print "  log/Center.log : centers"
    print ""
    print " e.g : python SearchCentroid.py rate-400-2300-all.fits 450 450 0.189 points.reg"
    exit()

#-------------------------------------------------
#setup parameter : maximum radius for flux weight center
z=float(argv[4])
Nite=100 #
#option parameter : r0=dRlimit, r1=r0+dRlimit, r2=r1+dRlimit, r7
# currently not used here
dRlimit=12. #30 arcsec

FITSFileName = argv[1]
if not os.path.isfile(FITSFileName) :
    print "Error:", FITSFileName, "does not exits"
    exit()
    
import pyfits
import pywcs
import numpy

RegionFileName = open("SetUpFile/Region.dat","r")
RadiusList = RegionFileName.readlines()

Rlimitkpc=500.
amtokpc=cos.kpc_proper_per_arcmin(z).value
dA=cos.angular_diameter_distance(z).value
degtocm=3600.*amtokpc*kpctocm
cmtodeg=1./degtocm
arcmin_kpc=Rlimitkpc/amtokpc*60.
Rlimit=arcmin_kpc
print "\n Rlimit = %s [arcsec]\n"%Rlimit

PointsourceFileName = argv[5]
if os.path.isfile("points_deg.reg") :
    os.system("rm points_deg.reg")
os.system("ds9 -xpa localhost &")
os.system("sleep 10")
os.system("cat %s-obj-image-sky.fits | xpaset ds9 fits"%mos1)
os.system("xpaset -p ds9 regions load %s"%PointsourceFileName)
os.system("xpaset -p ds9 scale log")
os.system("xpaset -p ds9 cmap b")
os.system("xpaset -p ds9 regions system image ")
os.system("xpaset -p ds9 regions save ./points_deg.reg")

PointsourceFileNamedeg = "./points_deg.reg"
PL = open(PointsourceFileNamedeg,"r")
PointsourceList = PL.readlines()
PointsourceList = [t.replace(',',' ').replace('circle(','').replace(') # color=white width=2\n','').replace(')','') for i,t in enumerate(PointsourceList) if i > 2 and i < len(PointsourceList)-1]
TF=[i.replace('\n','').isdigit() for i in RadiusList]
TFcount=TF.count(True)
output_name="log/Center.log"

#----------------------------------------------------------------------

HDUList = pyfits.open(FITSFileName)
RateImage = HDUList[0].data

NbinsX, NbinsY = len(RateImage[0]), len(RateImage)

CentroidX=float(argv[2])
CentroidY=float(argv[3])

print "\n#Calculating second-order moment at the position of CAMIRA cluster"
X2Image=[];
Y2Image=[];
ZImage=[];
PRORateImage=RateImage;
for ix in xrange(1, NbinsX+1) :
    for iy in xrange(1, NbinsY+1) :
        if (ix-CentroidX)**2 + (iy-CentroidY)**2 < Rlimit**2 :
            for l in PointsourceList :
                px,py,pr=numpy.array(l.split()).astype(np.float64);
                pxt=-px+2*CentroidX;
                pyt=-py+2*CentroidY;
                if (ix-pxt)**2 + (iy-pyt)**2 < pr**2 :
                    PRORateImage[iy-1,ix-1]=0.0;
for ix in xrange(1, NbinsX+1) :
    for iy in xrange(1, NbinsY+1) :
        if (ix-CentroidX)**2 + (iy-CentroidY)**2 < Rlimit**2 :
            X2Image.append((ix-CentroidX)**2);
            Y2Image.append((iy-CentroidY)**2);
            ZImage.append(PRORateImage[iy-1, ix-1]);
Second_XMoment=np.average(X2Image,weights=ZImage);
Second_YMoment=np.average(Y2Image,weights=ZImage);

print "#Second-order moment (X,Y):", Second_XMoment,Second_YMoment,"\n"
print "#Calculating third-order moment at CAMIRA Cluster"
X3Image=[];
Y3Image=[];
ZImage=[];
for ix in xrange(1, NbinsX+1) :
    for iy in xrange(1, NbinsY+1) :
        if (ix-CentroidX)**2 + (iy-CentroidY)**2 < Rlimit**2 :
            X3Image.append((ix-CentroidX)**3);
            Y3Image.append((iy-CentroidY)**3);
            ZImage.append(PRORateImage[iy-1, ix-1]);
Third_XMoment=np.average(X3Image,weights=ZImage);
Third_YMoment=np.average(Y3Image,weights=ZImage);
print "Third-order moment (X,Y): ",Third_XMoment, Third_YMoment,"\n"

skewnessX=Third_XMoment/(Second_XMoment**1.5);
skewnessY=Third_YMoment/(Second_YMoment**1.5);

print "Skewness (Third-order moment/((Second-order moment)**1.5) (X,Y)",skewnessX,skewnessY,"\n"

#first check within Rlimit arcsec from the center
print "#Searching center as a weight of flux within %f arcsec from (CentroidX, CentroidY)" % float(Rlimit*2.5)
for i in xrange(Nite) :
    XImage=[];
    YImage=[];
    ZImage=[];
    PRORateImage=RateImage
    RecordX=CentroidX
    RecordY=CentroidY
    for ix in xrange(1, NbinsX+1) :
        for iy in xrange(1, NbinsY+1) :
            if (ix-CentroidX)**2 + (iy-CentroidY)**2 < Rlimit**2 :
                for l in PointsourceList :
                    px,py,pr=numpy.array(l.split()).astype(np.float64);
                    pxt=-px+2*CentroidX;
                    pyt=-py+2*CentroidY;
                    if (ix-pxt)**2 + (iy-pyt)**2 < pr**2 :
                        PRORateImage[iy-1,ix-1]=0.0;
    for ix in xrange(1, NbinsX+1) :
        for iy in xrange(1, NbinsY+1) :
            if (ix-CentroidX)**2 + (iy-CentroidY)**2 < Rlimit**2 :
               XImage.append(ix);
               YImage.append(iy);
               ZImage.append(PRORateImage[iy-1, ix-1]);

    CentroidX=np.average(XImage,weights=ZImage);
    CentroidY=np.average(YImage,weights=ZImage);
    print "Iteration #",i," : ",CentroidX, CentroidY
    if (RecordX - CentroidX)**2 + (RecordY - CentroidY)**2 < 1.0*1.0 :
        break

"""
##opition
# iterative search within 360 arcsec, 330 arcsec ..., 210 arcsec
Rlimit=144.
for i in xrange(Nite) :
    if i%2==1:
        Rlimit-=dRlimit
    XImage=[];
    YImage=[];
    ZImage=[];
    for ix in xrange(1, NbinsX+1) :
        for iy in xrange(1, NbinsY+1) :
            if (ix-CentroidX)**2 + (iy-CentroidY)**2 < Rlimit**2 :
                XImage.append(ix);
                YImage.append(iy);
                ZImage.append(RateImage[iy-1, ix-1]);
    CentroidX=np.average(XImage,weights=ZImage);
    CentroidY=np.average(YImage,weights=ZImage);
    print "Iteration #",i," : ",CentroidX, CentroidY
"""

wcs = pywcs.WCS(HDUList[0].header)
pixcrd = numpy.array([[CentroidX, CentroidY]], numpy.float_)
sky = wcs.wcs_pix2sky(pixcrd, 1)
print "RA,Dec =", sky[0]

print "See : log/Center.log"

output=file(output_name,"w");
blnk=[""];
radec="RA,Dec = %s %s" % (sky[0][0], sky[0][1])
pixxy="PixX,PixY = %s %s" % (CentroidX, CentroidY)
np.savetxt(output,[radec],fmt='%s');
np.savetxt(output,[pixxy],fmt='%s');

#-----------------------------------------------------------------
for Prefix in (mos1, mos2, pn) :
  
    import commands
    
    output_esky2det = commands.getoutput("bin/esky2det2.csh %f %f %s-obj-image-sky.fits  2>/dev/null| gawk '{if(NF==2) print $0}'" % (sky[0][0], sky[0][1],Prefix))
    
    DetX = float(output_esky2det.split()[0])
    DetY = float(output_esky2det.split()[1])

    det="%s : DETX,DETY = %s %s" % (Prefix,DetX, DetY)
    np.savetxt(output,[det],fmt='%s');
    
    for i in xrange(TFcount-1) :
       if RadiusList[i].replace('\n','').isdigit(): 
           Rin, Rout = RadiusList[i], RadiusList[i+1]
           RegFileName = "%s_%s-%ssecR.reg" % (Prefix, Rin.replace('\n',''), Rout.replace('\n',''))
           RegFile = file(RegFileName, "w")
           print >>RegFile, "&&((DETX,DETY) IN circle(%f,%f,%f))&&!((DETX,DETY) IN circle(%f,%f,%f))"  % (DetX, DetY, 20.*float(Rout.replace('\n','')), DetX, DetY, 20.*float(Rin.replace('\n','')),)
Sndmoment="2nd moment (X,Y) = %s %s" % (Second_XMoment,Second_YMoment)
Trdmoment="3rd moment= %s %s" % (Third_XMoment,Third_YMoment)
Skewness="Skewness(3rd-order moment/((2nd-order moment)**1.5)= %s %s" % (skewnessX,skewnessY)
np.savetxt(output,[Sndmoment],fmt='%s');
np.savetxt(output,[Trdmoment],fmt='%s');
np.savetxt(output,[Skewness],fmt='%s');

output.close();

print "Have made region files for spec analysis"
