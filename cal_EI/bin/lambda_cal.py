#!/usr/bin/env python

from sys import *
import os
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
from scipy import integrate
import numpy as np
import array
from scipy.interpolate import interp1d

kpctocm=3.085677581e21
Mpctocm=3.085677581e24
Mp=1.6726231e-24
Msun=1.989e33
Munit=1e14

mu=0.5964
ratio=1./1.174417
ratio_HE=0.0872083
rho_factor=1.92574264507411 

H0=70.
Omega_m0=0.3;
Omega_L=0.7;
cos=FlatLambdaCDM(H0, Omega_m0, Neff=0)


if len(argv) != 4: 
   print "INPUT ERROR in python";
   print "python lamda_cal.py prefix z fake-spectrum"
   print " prefix               - type of detector "
   print " z                    - red shift "
   print " fake-spectrum        - qdp file "
   print ""
   print " output "
   print " result/Lamda.dat     - central cooling function = Lamda"
   print ""
   print " e.g : python lamda_cal.py mos1S001 0.1832 fake-data_mos1S001.qdp"
   exit();

prefix = argv[1]
z = float(argv[2])
mos1qdpfile = argv[3]
output_name="../result/Lambda_%s.dat"%prefix

amtokpc=cos.kpc_proper_per_arcmin(z).value
dA=cos.angular_diameter_distance(z).value

if os.path.isdir("../result") == False:
   os.mkdir("../result")
#------------------------------------------------------------------------------------------
#--------------
# radius definition
#-------------
Region=np.loadtxt("../spec_ana/SetUpFile/Region.dat",dtype='d',usecols=(0,),unpack=True);
Rdata=np.zeros(Region.size-2); # the outmost radii is background
R_mdata=np.zeros(Region.size-2);
R_pdata=np.zeros(Region.size-2);

for i in np.arange(Region.size-2) :
    Rdata[i]=(Region[i]+Region[i+1])/120.
    R_mdata[i]=Region[i]/60.
    R_pdata[i]=Region[i+1]/60.
#--------------------------------------------------------------------------------------------
outputfile=file(output_name,"w")

erangefile = open("../spec_ana/SetUpFile/eRange.dat","r")
erange = erangefile.readlines()
elow = float(erange[0])/1e3
ehigh = float(erange[1])/1e3

f1 = open(mos1qdpfile,"r")
d1 = f1.readlines()
#a=d1.pop(0)
#a=d1.pop(0)
#a=d1.pop(0)
c = d1.count("NO NO NO NO\n")
m = 0
Integral = np.zeros(c)
#----------------------------------------------------------------------------------------------
if c != Rdata.size :
   print "ERROR : the number of qdp output != the number of the annuli"
   print prefix
   print "Number of the annuli of ",mos1qdpfile,"    : ",c
   print "Number of the annuli of spectral fitting : ",Rdata.size
   exit();
#--------------------------------------------------------------------------------------------

if c == 0:
   m1e = d1[0:len(d1)]
   limm1e= array.array("f")
   limm1c= array.array("f")
   for i in xrange(len(m1e)):
      if elow-0.2 <=float(m1e[i].split()[0])<= ehigh+1.0:
         limm1e.append(float(m1e[i].split()[0]))
         limm1c.append(float(m1e[i].split()[2]))
   f1 = interp1d(limm1e, limm1c, kind='cubic')
   #f1 = interp1d(limm1e, limm1c, kind='quadratic')
   int1 = integrate.quad(lambda x: f1(x), elow, ehigh)[0]
   """
   Sum1 = f1(elow)+f1(ehigh)
   l = (ehigh-elow)/10000.
   for i in xrange(1,10000):
      if i%2. == 0:
         Sum1 += 2.*f1(i*l+elow)
      if i%2. != 0.:
         Sum1 += 4.*f1(i*l+elow)
   """
   #Integral[0]=float(Sum1*l/3.)
   Integral[i]=float(int1)
   print "    ",i
   print "    Spectrum Integral Value    = ",Integral[i]
   print " "
   print "    Cooling function Lambda    = ",Integral[i]*10**(-14)*(1+z)**2
   print "    Emission Measure Converter = ",Integral[i]/(4*np.pi*dA**2*(1+z)**2)*10**(-14)
   print " "


print c
for i in xrange(c):
   p = d1.index("NO NO NO NO\n")
   #print m,p
   m1e = d1[m:p]
   m = d1.index("NO NO NO NO\n")
   d1.pop(p)
   limm1e= array.array("f")
   limm1c= array.array("f")
   for k in xrange(len(m1e)):
      #print m1e[k].split()[0]
      if elow-0.2 <=float(m1e[k].split()[0])<= ehigh+1.0:
         #print float(m1e[k].split()[0]),float(m1e[k].split()[2])
         limm1e.append(float(m1e[k].split()[0]))
         limm1c.append(float(m1e[k].split()[2]))
   f1 = interp1d(limm1e, limm1c, kind='cubic')
   #f1 = interp1d(limm1e, limm1c, kind='quadratic')
   int1 = integrate.quad(lambda x: f1(x), elow, ehigh)[0]
   """
   Sum1 = f1(elow)+f1(ehigh)
   l = (ehigh-elow)/10000.
   for j in xrange(1,10000):
      if j%2. == 0:
         Sum1 += 2.*f1(j*l+elow)
      if j%2. != 0.:
         Sum1 += 4.*f1(j*l+elow)
   Integral[i]=Sum1*l/3.
   """ 
   Integral[i]=int1
   print "    ",i
   print "    Spectrum Integral Value    = ",Integral[i]
   print " "
   print "    Cooling function Lambda    = ",Integral[i]*10**(-14)*(1+z)**2
   print "    Emission Measure Converter = ",Integral[i]/(4*np.pi*dA**2*(1+z)**2)*10**(-14)
   print " "

#---------------------------------------------------------------------------
label="# r[arcmin] r_m[arcmin] r_p[arcmin] A[counts/sec*cm^5]"
data=np.array([Rdata,R_mdata,R_pdata,Integral])

np.savetxt(outputfile, [label], fmt='%s');
np.savetxt(outputfile, data.T, fmt='%0.8e',delimiter=' ', newline='\n');


outputfile.close();
