#!/usr/bin/env python

from sys import*
import os
import commands
import numpy as np

#-------------------------------------------------
mos1=commands.getoutput("head -1 log/prefix.log | gawk '{print $1}'")
mos2=commands.getoutput("head -2 log/prefix.log | tail -1 | gawk '{print $1}'")
pn=commands.getoutput("tail -1 log/prefix.log  | gawk '{print $1}'")
#-------------------------------------------------
region=np.loadtxt("SetUpFile/Region.dat",dtype='d',usecols=(0,),unpack=True);
print ""
print region
print ""
#MOS
cal_cfile=open("bin/shoda_mos.csh","w")
cal_cfile.write("#!/bin/csh -f\n")
cal_cfile.write("set Prefix  = $argv[1]\n")
cal_cfile.write("xspec <<EOF\n")
for i in xrange(1,len(region)):
   annulus="%.0lf-%.0lf"%(region[i-1],region[i])
   sentence="data %s:%s ${Prefix}-obj-%s-grp.pi\n"%(i,int(i),annulus)
   cal_cfile.write("%s"%sentence)
cal_cfile.write("statistic chi\n")
cal_cfile.write("setpl ene\n")
for i in xrange(len(region)):
   cal_cfile.write("ignore %s:0.0-0.3,11.0-**\n"%int(i+1))
cal_cfile.write("ignore bad\n")
cal_cfile.write("method leven 10 0.01\n")
cal_cfile.write("abund angr\n")
cal_cfile.write("xsect bcmc\n")
cal_cfile.write("cosmo 70 0 0.7\n")
cal_cfile.write("xset delta -1\n")
cal_cfile.write("systematic 0\n")
cal_cfile.write("log ${Prefix}-specdata.log\n")
cal_cfile.write("show data\n")
cal_cfile.write("log none\n")
cal_cfile.write("y\n")
cal_cfile.write(" \n")
cal_cfile.write("<<EOF\n")
cal_cfile.close()

#PN
cal_cfile=open("bin/shoda_pn.csh","w")
cal_cfile.write("#!/bin/csh -f\n")
cal_cfile.write("set Prefix  = $argv[1]\n")
cal_cfile.write("xspec <<EOF\n")
for i in xrange(1,len(region)):
   annulus="%.0lf-%.0lf"%(region[i-1],region[i])
   sentence="data %s:%s ${Prefix}-obj-os-%s-grp.pi\n"%(i,int(i),annulus)
   cal_cfile.write("%s"%sentence)

cal_cfile.write("statistic chi\n")
cal_cfile.write("setpl ene\n")
for i in xrange(len(region)):
   cal_cfile.write("ignore %s:0.0-0.4,11.0-**\n"%int(i+1))
cal_cfile.write("ignore bad\n")
cal_cfile.write("method leven 10 0.01\n")
cal_cfile.write("abund angr\n")
cal_cfile.write("xsect bcmc\n")
cal_cfile.write("cosmo 70 0 0.7\n")
cal_cfile.write("xset delta -1\n")
cal_cfile.write("systematic 0\n")
cal_cfile.write("log ${Prefix}-specdata.log\n")
cal_cfile.write("show data\n")
cal_cfile.write("log none\n")
cal_cfile.write("y\n")
cal_cfile.write(" \n")
cal_cfile.write("<<EOF\n")
cal_cfile.close()

os.system('csh bin/shoda_mos.csh %s' % (mos1))
os.system('csh bin/shoda_mos.csh %s' % (mos2))
os.system('csh bin/shoda_pn.csh %s' % (pn))

for Prefix in (mos1,mos2,pn): 
   SpecLogFile = open("%s-specdata.log"%Prefix,"r")
   speclog = SpecLogFile.readlines()
   fnum = speclog.index("!XSPEC12> show data\n")
   RegionFile = open("SetUpFile/Region.dat","r")
   region = RegionFile.readlines()
   outdata = ["# region counts/sec exposure counts"]
   csum=0
   for i in xrange(len(region)-1):
      reg = "%s-%s"%("{0:4d}".format(int(region[i].replace("\n",""))),int(region[i+1].replace("\n","")))
      rate = float(speclog[fnum+4+13*i].split()[6])
      rateerr = float(speclog[fnum+4+13*i].split()[8])
      exp = float(speclog[fnum+8+13*i].split()[3])
      counts = float(speclog[fnum+4+13*i].split()[6])*float(speclog[fnum+8+13*i].split()[3])
      countserr = float(speclog[fnum+4+13*i].split()[8])*exp
      outdata.append(" %s %f +/- %f %f %f +/- %f"%(reg,rate,rateerr,exp,counts,countserr))
      csum=csum+counts
   outdata.append("Sum %s"%csum)
   output = open("calc_%s.dat"%Prefix,"w")
   for i in xrange(len(region)+1):
      output.write("%s\n"%outdata[i])
      print outdata[i]
   output.close()
      










