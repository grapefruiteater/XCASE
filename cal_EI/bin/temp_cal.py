#!/usr/bin/env python

from sys import*
import os
import commands

import pyfits,array 
import numpy as np
from scipy import integrate
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM

if len(argv) != 4 : 
   print ""
   print "INPUT ERROR in python";
   print "python ne_cal.py z save-all save-error"
   print " z         - redshift "
   print " save-all  - spectrum fit data "
   print " save-error- temp error data "
   print ""
   print " output "
   print "     result/temp.dat    - temperature profile"
   print "     result/norm.dat    - normalization profile"
   print "     result/abun.dat    - abundance profile"
   print ""
   print " e.g : python temp_cal.py 0.1832 save-all.xcm observed_temp_error.dat"
   exit();

mos1=commands.getoutput("head -1 ../spec_ana/log/prefix.log | gawk '{print $1}'")
#mos2=commands.getoutput("head -2 log/prefix.log | tail -1 | gawk '{print $1}'")
#pn=commands.getoutput("tail -1 log/prefix.log | gawk '{print $1}'")

Munit=1e14
Msun=1.989e33
kpctocm=3.085677581e21
Mpctocm=3.085677581e24
cv=2.99792458e10
el=4.80298e-10
kB=1.60217646e-9
h=6.62606876e-27
Mp=1.6726231e-24
Me=9.1093897e-28
Ktoerg=1.380658e-16
KtokeV=8.617385e-8
Msun=1.989e33
G=6.6726e-8
yrtosec=3.155815e7

H0=70.
Omega_m0=0.3;
Omega_L=0.7;
cos=FlatLambdaCDM(H0, Omega_m0, Neff=0)

mu=0.5964
ratio=1./1.174417 
ratio_HE=0.0872083
rho_factor=1.92574264507411 

z = float(argv[1])
datafile = argv[2]
errorfile = argv[3]

d_A=cos.angular_diameter_distance(z).value*1e3
amtokpc=cos.kpc_proper_per_arcmin(z).value

if os.path.isdir("../result") == False:
   os.mkdir("../result")

output_name="../result/temp.dat"
outputfile=file(output_name,"w");

output_name1="../result/norm.dat"
outputfile1=file(output_name1,"w");

output_name2="../result/abun.dat"
outputfile2=file(output_name2,"w");

output_name3="../result/constant.dat"
outputfile3=file(output_name3,"w");

spec_file = open(datafile,"r")
spec_data = spec_file.readlines()

#Region = []
RegionList=open("../spec_ana/SetUpFile/Region.dat")
Region=RegionList.readlines()

#pilist = [x for x in spec_data if mos1 in x]
#for i in xrange(len(pilist)):
#   b1 = pilist[i].split()
#   b2 = b1[2].split("-")
#   #print b2[2]
#  Region.append(b2[2])
#   if i==len(pilist)-1:
#      Region.append(b2[3])
#      #print b2[3]

Para = []
Tdata = array.array("f")
Rdata = array.array("f")
R_mdata = array.array("f")
R_pdata = array.array("f")
yerrorl = array.array("f")
yerrorh = array.array("f")
xerrorl = array.array("f")
xerrorh = array.array("f")
Norm = array.array("f")
Abun = array.array("f")
C1 = array.array("f")
C2 = array.array("f")
C3 = array.array("f")

firstnumber=int(commands.getoutput("gawk '{if($1==\"model\") print NR}' %s | head -1" % datafile))-1
#firstnumber = spec_data.index("model  gaussian + gaussian + gaussian + gaussian + gaussian + gaussian + gaussian + constant*constant(gaussian + gaussian + apec + (apec + apec + powerlaw)phabs + apec*phabs)\n")

for i in xrange(len(spec_data)-firstnumber):
    n = spec_data.pop(firstnumber)
    Para.append(n)
    
error_file = open(errorfile,"r")
errorinfo  = error_file.readlines()
for i in xrange(len(errorinfo)):
    e = errorinfo[i].split()
    yerrorl.append(float(e[1]))
    yerrorh.append(float(e[2]))
    xerrorl.append((float(Region[i+1])-float(Region[i]))/120)
    xerrorh.append((float(Region[i+1])-float(Region[i]))/120)


for i in xrange(len(errorinfo)):
    n = (float(Region[i])+float(Region[i+1]))/120.#/3600*pi/180*d_A
    Rdata.append(n)
    R_mdata.append(float(Region[i])/60.)
    R_pdata.append(float(Region[i+1])/60.)

for i in xrange(len(errorinfo)):
    p = Para[45+147*i].split()
    Tdata.append(float(p[0]))
    A = Para[46+147*i].split()
    if A[0]=='=':
       Abun.append(Abun[i-1])
    else:
       Abun.append(float(A[0]))
    N = Para[48+147*i].split()
    Norm.append(float(N[0]))
    
C1.append(float(Para[22].split()[0]))
C2.append(float(1.0))
C3.append(float(Para[120].split()[0]))
#print len(Rdata),len(Tdata),len(yerrorl),len(Abun),len(Norm)

#output
line=np.array([Rdata,R_mdata,R_pdata,Tdata,yerrorl,yerrorh]);
label="# r[arcmin] r_m[arcmin] r_p[arcmin] T[keV] T_m T_p" 
np.savetxt(outputfile, [label], fmt='%s');
np.savetxt(outputfile, line.T, fmt='%0.8e',delimiter=' ', newline='\n');

line=np.array([Rdata,R_mdata,R_pdata,Norm]);
label="# r[arcmin] r_m[arcmin] r_p[arcmin] Norm[1/cm^5]" 
np.savetxt(outputfile1, [label], fmt='%s');
np.savetxt(outputfile1, line.T, fmt='%0.8e',delimiter=' ', newline='\n');

line=np.array([Rdata,R_mdata,R_pdata,Abun]);
label="# r[arcmin] r_m[arcmin] r_p[arcmin] Abundance" 
np.savetxt(outputfile2, [label], fmt='%s');
np.savetxt(outputfile2, line.T, fmt='%0.8e',delimiter=' ', newline='\n');

line=np.array([C1,C2,C3]);
label="# MOS1 MSO2 PN" 
np.savetxt(outputfile3, [label], fmt='%s');
np.savetxt(outputfile3, line.T, fmt='%0.8e',delimiter=' ', newline='\n');
