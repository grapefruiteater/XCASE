#!/usr/bin/env python

"""
install colossus
http://www.benediktdiemer.com/code/
"""

camiracat="../CAMIRA/camira_s16a_wide_v2.dat"
memcat="../CAMIRA/camira_mem_s16a_wide_v2.dat"

import sys
import os
import numpy as np
import pyfits
import matplotlib.pyplot as plt;
import matplotlib.transforms as mtransforms
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.mlab import griddata
from matplotlib import cm
import matplotlib.gridspec as gridspec
import scipy.optimize
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
import scipy.integrate as integrate
from colossus.cosmology import cosmology
from colossus.halo import concentration

kpctocm=3.085677581e21
Mpctocm=3.085677581e24
Mp=1.6726231e-24
Msun=1.989e33
Munit=1e14
G=6.6726e-8 
kB=1.60217646e-9
degtorad=np.pi/180.;
radtodeg=180./np.pi;
degtoam=60.;

if len(sys.argv) != 7 : 
   print "pytyon cal_Ms.py ra_camira dec_camira ra_measure dec_measure r[arcmin] Delta";
   print "   ra_camira   - camira ra [deg] : exactly same value"
   print "   dec_camira  - camira dec [deg] : exactly same value"
   print "   ra_measure  - measurement center [deg] : e.g. xmm, wl"
   print "   dec_measure - measurement center [deg] : e.g. xmm, wl"
   print "   r           - within radius [arcmin]"
   print "   Delta       - Overdensity [2500c, 1000c,500c,200c,vir,180m,200m]"
   print "e.g. :  python cal_Ms.py  29.300512 -5.918173 29.293548  -5.869181 6.87286433 500"
   sys.exit();


ra_camira=float(sys.argv[1])
dec_camira=float(sys.argv[2])
ra_xmm=float(sys.argv[3])
dec_xmm=float(sys.argv[4])
r500=float(sys.argv[5])
Delta=float(sys.argv[6])

"""
ra_camira=29.300512
dec_camira=-5.918173
ra_xmm=29.293548
dec_xmm=-5.869181
r500=6.87286433 # arcmin
Delta=500
"""

#====================================================================================================================
def kappa (r,rs) :
      x=r/rs;

      if x<1 :
          f1_NFW=( 1.-2.*np.arctanh(np.sqrt((1.-x)/(1.+x)))/np.sqrt( 1.-x*x ) )/( x*x-1. );
      elif x==1 : 
 	  f1_NFW=1./3.;
      elif x > 1 :
          f1_NFW=(1.-2.*np.arctan(np.sqrt((x-1.)/(x+1.)))/np.sqrt(x*x-1))/(x*x-1.);
       
      return f1_NFW;

def filter (r,n,R0_am) :
    f1=scipy.special.gammaincc(n/2.,(r/R0_am)**2)-((r/R0_am)**n)*np.exp(-(r/R0_am)**2)
    f2=scipy.special.gammaincc(n/2.,0.)
    return f1/f2;
 

def cal_conversionF (rout,rs,n, R0_am) :
    kappaall=integrate.quad(lambda x:  kappa(x,rs)*filter(x,n,R0_am)*x, 0., rout)[0]
    norm=integrate.quad(lambda x:  kappa(x,rs)*filter(x,n,R0_am)*x, 0., 10*R0_am)[0]
    return kappaall/norm;

def cal_conversionF_off_theta (x,rs,n, R0_am,doff) :
    return integrate.quad(lambda theta: kappa(np.sqrt(x*x+doff*doff-2*x*doff*np.cos(theta)),rs)*filter(np.sqrt(x*x+doff*doff-2*x*doff*np.cos(theta)),n,R0_am),0,2*np.pi)[0]/2/np.pi

def cal_conversionF_off (rout ,rs, n, R0_am,doff) :
    kappaall=integrate.quad(lambda x: cal_conversionF_off_theta(x,rs,n,R0_am,doff)*x, 0., rout)[0]
    norm=integrate.quad(lambda x:  kappa(x,rs)*filter(x,n,R0_am)*x, 0., 10*R0_am)[0]
    return kappaall/norm;


def distance(ra1,dec1,ra2,dec2) :
    return np.arccos(np.sin(dec1*degtorad)*np.sin(dec2*degtorad)+np.cos(dec1*degtorad)*np.cos(dec2*degtorad)*np.cos((ra1-ra2)*degtorad))*radtodeg*degtoam;

#====================================================================================================================
doff=distance(ra_camira,dec_camira,ra_xmm,dec_xmm); # off-centering distance [arcmin]

rac_ori,decc_ori,z_ori,Ncor_ori,logMs_ori,z_bcg=np.loadtxt(camiracat,dtype='d',usecols=(0,1,2,3,4,5),unpack=True);
rac,decc,Ncor,z,ra,dec,logMs,wgt=np.loadtxt(memcat,dtype='d',usecols=(0,1,2,3,4,5,6,7),unpack=True);
#--------------------------------------------------------------
# cut
sel=[rac_ori==ra_camira,decc_ori==dec_camira]
flags = np.all(sel, axis = 0)

rac_ori=rac_ori[flags][0]
decc_ori=decc_ori[flags][0]
Ncor_ori=Ncor_ori[flags][0]
z_ori=z_ori[flags][0]
logMs_ori=logMs_ori[flags][0]

#---
sel=[rac==ra_camira,decc==dec_camira]
flags = np.all(sel, axis = 0)

rac=rac[flags]
decc=decc[flags]
Ncor=Ncor[flags]
z=z[flags]
ra=ra[flags]
dec=dec[flags]
logMs=logMs[flags]
wgt=wgt[flags]
#------------------------------------------------------------
# cosmology #note slightly different from Oguri+
H0=70.
Om0=0.3; #Omega_m0
OL=0.7;  #Omega_L
cos=FlatLambdaCDM(H0, Om0, Neff=0)
amtokpc=cos.kpc_proper_per_arcmin(z_ori).value
amtoMpc=amtokpc/1e3;
rho_cr=cos.critical_density(z_ori).value
Omz=cos.Om(z_ori);
#------------------------------------------------------------
# Oguri+
H0_Oguri=70.
Om0_Oguri=0.28
OL_Oguri=0.72;
cos_Oguri=FlatLambdaCDM(H0_Oguri, Om0_Oguri, Neff=0)
#------------------------------------------------------------
# filiter
n=4.;
R0=0.8*(H0/100.)**-1 # Mpch^-1
R0_am=R0/amtoMpc
#-----------------------------------------------------------------
M500=(4.*np.pi/3.)*Delta*rho_cr*(r500*amtokpc*kpctocm)**3 # g
cosmology.setCosmology('planck15')
Deltalabel="%.0lfc" % Delta
c500=concentration.concentration(M500*(100/H0)**-1/Msun, Deltalabel, z_ori, model='diemer15', statistic='median', conversion_profile='nfw')
rs=r500/c500   #arcmin
#-----------------------------------------------------------------
# 1 : conversion Factor method
F=cal_conversionF(r500,rs,n, R0_am)
Ms_radius=(10**logMs_ori)*F
logMs_radius=np.log10(Ms_radius);

#-----------------------------------------------------------------
# 2 : conversion Factor method with xmm center
F_off=cal_conversionF_off(r500,rs,n, R0_am,doff)
Ms_radius_off=(10**logMs_ori)*F_off
logMs_radius_off=np.log10(Ms_radius_off);

#-----------------------------------------------------------------
# 3 : sum member galaxies within radius : camira center
d_camira=distance(ra,dec,ra_camira,dec_camira)
#np.arccos(np.sin(dec*degtorad)*np.sin(dec_camira*degtorad)+np.cos(dec*degtorad)*np.cos(dec_camira*degtorad)*np.cos((ra-ra_camira)*degtorad))*radtodeg*degtoam;
selgal_camira=[d_camira<r500]
flags_camira = np.all(selgal_camira, axis = 0)

logMs_camira=logMs[flags_camira]
wgt_camira=wgt[flags_camira]

Ms_camiraCen=np.sum(10**(logMs_camira)*wgt_camira)
logMs_camiraCen=np.log10(Ms_camiraCen)

#----------------------------------------------------------------
# 4 : sum member galaxies within radius : given center
d_xmm=np.arccos(np.sin(dec*degtorad)*np.sin(dec_xmm*degtorad)+np.cos(dec*degtorad)*np.cos(dec_xmm*degtorad)*np.cos((ra-ra_xmm)*degtorad))*radtodeg*degtoam;
selgal_xmm=[d_xmm<r500]
flags_xmm = np.all(selgal_xmm, axis = 0)

logMs_xmm=logMs[flags_xmm]
wgt_xmm=wgt[flags_xmm]

Ms_xmmCen=np.sum(10**(logMs_xmm)*wgt_xmm)
logMs_xmmCen=np.log10(Ms_xmmCen)
#------------------------------------------------
# print

print "# logMs rDelta[arcmin] z Delta MDelta[10^14Msun] cDelta conversionFactor"
print "#",logMs_ori,r500, z_ori,Delta, M500/Msun/Munit, c500, F,F_off
print "# method logMs(<r) Ms(<r)[Msun]"
print "conversionFactor(camira) ",logMs_radius, Ms_radius
print "conversionFactor(given) ",logMs_radius_off, Ms_radius_off
print "camira_center(2D)    ", logMs_camiraCen,Ms_camiraCen
print "given_center(2D)     ", logMs_xmmCen,Ms_xmmCen

