#!/usr/bin/env python

"""
install colossus
http://www.benediktdiemer.com/code/
"""
doff_max=5. # maximum radius arcmin search camira cluster 
camiracat="../CAMIRA/camira_s16a_wide_v2.dat"

import sys
import os
import numpy as np
import math
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
from scipy.interpolate import UnivariateSpline

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

#--------------------------------------------------------------------------------
if len(sys.argv) != 7 : 
   print "pytyon cal_Ms_fromTable.py root ra dec z masstable dir";
   print "   root        - xray or wl"
   print "   ra          - center ra [deg]"
   print "   dec         - center dec [deg]"
   print "   z           - redshift"
   print "   masstable   - mass table for overdensity"
   print "   dir         - directory to read Mgas_profile.dat/fgas_profile.dat"
   print " "
   print " output  :  dir/fb_root.dat (e.g MHEresult/fb_xray.dat)"
   print "            for Delta >=500 "
   print "e.g. :  python cal_Ms_fromTable.py xray 29.293548 -5.869181 0.1284 MHEresult/MHE_Delta.dat MHEresult"
   sys.exit();

#--------------------------------------------------------------------------------
root=sys.argv[1]
ra_xmm=float(sys.argv[2])
dec_xmm=float(sys.argv[3])
z_xmm=float(sys.argv[4])
MHE_data=sys.argv[5]
dir=sys.argv[6]

if root == "xray" or root == "wl" :
   print "used",root,"data"
else :
   print "ERROR : root should be xray or wl"
   sys.exit()
#-------------------------------------------
Mgas_data=dir+"/Mgas_profile.dat"
fgas_data=dir+"/fgas_profile.dat"


Delta1_label=np.loadtxt(MHE_data,dtype='string',usecols=(0,),unpack=True);
Delta1,r1,r1_m,r1_p,r1_am,r1_am_m,r1_am_p,MHE, MHE_m,MHE_p=np.loadtxt(MHE_data,dtype='d',usecols=(1,2,3,4,5,6,7,8,9,10),unpack=True);
rprof_Mgas,Mgas_prof,Mgas_prof_m,Mgas_prof_p=np.loadtxt(Mgas_data,dtype='d',usecols=(0,1,2,3),unpack=True);
rprof_fgas,fgas_prof,fgas_prof_m,fgas_prof_p=np.loadtxt(fgas_data,dtype='d',usecols=(0,1,2,3),unpack=True);

# cal Mgas and fgas at given Mass and radius
func_Mgas=UnivariateSpline(rprof_Mgas, Mgas_prof, k=3,s=0);
func_Mgas_m=UnivariateSpline(rprof_Mgas, Mgas_prof_m, k=3,s=0);
func_Mgas_p=UnivariateSpline(rprof_Mgas, Mgas_prof_p, k=3,s=0);
func_fgas=UnivariateSpline(rprof_fgas, fgas_prof, k=3,s=0);
func_fgas_m=UnivariateSpline(rprof_fgas, fgas_prof_m, k=3,s=0);
func_fgas_p=UnivariateSpline(rprof_fgas, fgas_prof_p, k=3,s=0);

Mgas=func_Mgas(r1);
Mgas_m=func_Mgas_m(r1_m);
Mgas_p=func_Mgas_p(r1_p);

if root == "xray" : # note that there is a correlation between MHE and Mgas 
   fgas=func_fgas(r1);
   fgas_m=func_fgas_m(r1_m);
   fgas_p=func_fgas_p(r1_p);
elif root == "wl" : # note that there is no correlation between Mwl and Mgas 
   fgas=func_Mgas(r1)/MHE;
   fgas_m=func_Mgas_m(r1_m)/MHE_p;
   fgas_p=func_Mgas_p(r1_p)/MHE_m;


#====================================================================================================================

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
rac_ori,decc_ori,z_ori,Ncor_ori,logMs_ori,z_bcg=np.loadtxt(camiracat,dtype='d',usecols=(0,1,2,3,4,5),unpack=True);
#--------------------------------------------------------------
# cut
d=distance(rac_ori,decc_ori,ra_xmm,dec_xmm)
sel=[d<doff_max,np.abs(z_ori-z_xmm)<0.1*(1+z_xmm)]
flags = np.all(sel, axis = 0)

if rac_ori[flags].size == 0 :
   print "ERROR : no CAMIRA data"
   sys.exit;
elif rac_ori[flags].size > 1 :
   print "ERROR : multi CAMIRA objects detected"
   rac_ori=rac_ori[flags]
   decc_ori=decc_ori[flags]
   for i in np.arange(rac_ori.size) :
       print "ra, dec : ",rac_ori[i],dec_ori[i];
   sys.exit

rac_camira=rac_ori[flags][0]
decc_camira=decc_ori[flags][0]
Ncor_ori=Ncor_ori[flags][0]
z_ori=z_ori[flags][0]
logMs_ori=logMs_ori[flags][0]

doff=distance(rac_camira,decc_camira,ra_xmm,dec_xmm); # off-centering distance [arcmin]
#------------------------------------------------------------
# cosmology #note slightly different from Oguri+
H0=70.
Om0=0.3; #Omega_m0
OL=0.7;  #Omega_L
cos=FlatLambdaCDM(H0, Om0, Neff=0)
amtokpc=cos.kpc_proper_per_arcmin(z_xmm).value
amtoMpc=amtokpc/1e3;
rho_cr=cos.critical_density(z_xmm).value
Omz=cos.Om(z_xmm);
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
#----------------------------------------------------------------------------------------------------------------
#for i in np.arange(r1_am.size) :
# cosmology for concentration
# we use single rs even when MHE changes
cosmology.setCosmology('planck15')
c500=concentration.concentration(MHE[Delta1==500]*(100/H0)**-1*Munit, "500c", z_xmm, model='diemer15', statistic='median', conversion_profile='nfw')
rs=r1_am[Delta1==500]/c500

c500_p=concentration.concentration(MHE_p[Delta1==500]*(100/H0)**-1*Munit, "500c", z_xmm, model='diemer15', statistic='median', conversion_profile='nfw')
rs_p=r1_am_p[Delta1==500]/c500_p

c500_m=concentration.concentration(MHE_m[Delta1==500]*(100/H0)**-1*Munit, "500c", z_xmm, model='diemer15', statistic='median', conversion_profile='nfw')
rs_m=r1_am_m[Delta1==500]/c500_m

#----------------------------------------------------------------------------
# conversion Factor with offentering 


#only for Delta >= 500
Nsize=r1_am[Delta1>=500].size
flags=[Delta1>=500.]

Delta1_label=Delta1_label[flags]
Delta1=Delta1[flags]
r1=r1[flags]
r1_m=r1_m[flags]
r1_p=r1_p[flags]
r1_am=r1_am[flags]
r1_am_m=r1_am_m[flags]
r1_am_p=r1_am_p[flags]
MHE=MHE[flags]
MHE_m=MHE_m[flags]
MHE_p=MHE_p[flags]
Mgas=Mgas[flags]
Mgas_m=Mgas_m[flags]
Mgas_p=Mgas_p[flags]
fgas=fgas[flags]
fgas_m=fgas_m[flags]
fgas_p=fgas_p[flags]



F_off=np.zeros(r1_am.size);
F_off_m=np.zeros(r1_am.size);
F_off_p=np.zeros(r1_am.size);


for i in np.arange(r1_am.size) :
    F_off[i]=cal_conversionF_off(r1_am[i],rs,n, R0_am,doff)
    F_off_m[i]=cal_conversionF_off(r1_am_m[i],rs_m,n, R0_am,doff)
    F_off_p[i]=cal_conversionF_off(r1_am_p[i],rs_p,n, R0_am,doff)
    #print F_off[i],F_off_m[i],F_off_p[i]

Ms=(10**logMs_ori)*F_off/Munit
logMs=np.log10(Ms*Munit);
Ms_p=(10**logMs_ori)*F_off_p/Munit
logMs_p=np.log10(Ms_p*Munit);
Ms_m=(10**logMs_ori)*F_off_m/Munit
logMs_m=np.log10(Ms_m*Munit);

fs=Ms/MHE
fs_m=Ms_p/MHE_p #note MHE_p < MHE_m --> fgas_p > fgas_m because of Ms is almost const and we used the same radius for Ms measurement
fs_p=Ms_m/MHE_m 

fb=fgas+fs
fb_m=fgas_m+fs_m
fb_p=fgas_p+fs_p

#---------------------------------------------------------------------------
output_fgas_Delta=dir+"/fb_%s.dat" % root
file_fgas_Delta=file(output_fgas_Delta,"w");
label=np.zeros(100);
# r_m r_p ram[arcmin] ram_m ram_p Mgas[10^%dMsunh_%d^{-5/2}] Mgas_m Mgas_p
#r[Mpch_%d^{-1}] fgas[h_%d^{-3/2}] fgas_m fgas_p" % (H0,H0) 
line2=np.array([Delta1, r1, r1_m, r1_p, r1_am, r1_am_m, r1_am_p,MHE,MHE_m,MHE_p,Mgas,Mgas_m,Mgas_p,Ms,Ms_m,Ms_p,fgas,fgas_m,fgas_p,fs,fs_m,fs_p,fb,fb_m,fb_p]);
np.savetxt(file_fgas_Delta, ["# %s mass" % root ], fmt='%s');
np.savetxt(file_fgas_Delta, ["# camira catalog"], fmt='%s');
np.savetxt(file_fgas_Delta, ["# ra dec Ncor z logMs"], fmt='%s');
line3="# %0.8e %0.8e %0.8e %0.8e %0.8e" % (rac_camira,decc_camira,Ncor_ori,z_ori,logMs_ori)
np.savetxt(file_fgas_Delta, [line3], fmt='%s',delimiter=' ', newline='\n');
line4="# offset distance [arcmin] : %.8e" % (doff)
np.savetxt(file_fgas_Delta, [line4], fmt='%s',delimiter=' ', newline='\n');
np.savetxt(file_fgas_Delta, ["# 1  Delta"], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 2  r      [Mpch_%d^{-1}]" % H0], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 3  r_m    [Mpch_%d^{-1}]" % H0], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 4  r_p    [Mpch_%d^{-1}]" % H0], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 5  r      [arcmin]"], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 6  r_m    [arcmin]"], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 7  r_p    [arcmin]"], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 8  Mtot   [10^%dMsunh_%d^{-1}]" % (math.log10(Munit),H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 9  Mtot_m [10^%dMsunh_%d^{-1}]" % (math.log10(Munit),H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 10 Mtot_p [10^%dMsunh_%d^{-1}]" % (math.log10(Munit),H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 11 Mgas   [10^%dMsunh_%d^{-5/2}]" % (math.log10(Munit),H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 12 Mgas_m [10^%dMsunh_%d^{-5/2}]" % (math.log10(Munit),H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 13 Mgas_p [10^%dMsunh_%d^{-5/2}]" % (math.log10(Munit),H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 14 Ms     [10^%dMsunh_%d^{-2}]" % (math.log10(Munit),H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 15 Ms_m   [10^%dMsunh_%d^{-2}]" % (math.log10(Munit),H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 16 Ms_p   [10^%dMsunh_%d^{-2}]" % (math.log10(Munit),H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 17 fgas   [h_%d^{-3/2}]" % (H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 18 fgas_m [h_%d^{-3/2}]" % (H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 19 fgas_p [h_%d^{-3/2}]" % (H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 20 fs     [h_%d^{-1}]" % (H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 21 fs_m   [h_%d^{-1}]" % (H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 22 fs_p   [h_%d^{-1}]" % (H0)], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 23 fb     "], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 24 fb_m   "], fmt='%s');
np.savetxt(file_fgas_Delta, ["# 25 fb_p   "], fmt='%s');


np.savetxt(file_fgas_Delta, line2.T, fmt='%0.8e',delimiter=' ', newline='\n');
file_fgas_Delta.close()


print "output file :",output_fgas_Delta
