#!/usr/bin/env python

import sys
import os
import pyfits,glob,math
import numpy as np
import commands

if len(sys.argv) != 4 : 
   print "INPUT ERROR in python";
   print "python make_fake_xcm.py prefix z save-all.xcm"
   print " prefix       - type of detector "
   print " z            - redshift "
   print " save-all.xcm - specturm log file"
   print ""
   print " output "
   print "   fake_prefix.xcm  - fake script"
   print ""
   print " e.g : python make_fake_xcm.py mos1S001 0.1832 save-all.xcm "
   sys.exit();

prefix = sys.argv[1]
z = sys.argv[2]
filename = sys.argv[3]

Para = []
file1 = open(filename,"r")
parainfo = file1.readlines()


firstnumber=int(commands.getoutput("gawk '{if($1==\"model\") print NR}' %s | head -1" % filename))-1
#firstnumber2= parainfo.index("model  gaussian + gaussian + gaussian + gaussian + gaussian + gaussian + gaussian + constant*constant(gaussian + gaussian + apec + (apec + apec + powerlaw)wabs + apec*wabs)\n")
for i in xrange(len(parainfo)-firstnumber):
    n = parainfo.pop(firstnumber)
    Para.append(n)

Region=np.loadtxt("../spec_ana/SetUpFile/Region.dat",dtype='string',usecols=(0,),unpack=True);
pilist1=[]
for i in xrange(len(Region)-2):
   reg = "%s-%s"%(Region[i],Region[i+1])
   pilist1.append("%s-obj-%s-grp_noback.pi" % (prefix, reg));
   os.system("bin/grppha.csh %s %s" % (prefix,reg))

#cd = commands.getoutput("ls -tr ../data").split()
#pilist = [x for x in cd if 'noback' in x]
#pilist1 = [x for x in pilist if '%s'%prefix in x]
T = []
A = []
phabs = Para[44].split()[0]
for i in xrange(len(Region)-2):
    t = Para[45+147*i].split()
    T.append(t[0])
    a =  Para[46+147*i].split()
    if a[0] == "=":
        A.append(A[i-1])
    else:
        A.append(a[0])
    
outfile = []
fakedata=["#!/bin/csh -f\n","xspec <<EOF\n","statistic chi\n","\n","setpl ene\n"]
for i in xrange(len(pilist1)):
   fakedata.append("ignore %s:0.0-0.3,11.0-**\n"%(i+1))
fakedata.append("ignore bad\n")
fakedata.append("method leven 10 0.01\n")
fakedata.append("abund angr\n")
fakedata.append("xsect bcmc\n")
fakedata.append("cosmo 70 0 0.7\n")
fakedata.append("xset delta -1\n")
fakedata.append("systematic 0\n")


for i in xrange(len(pilist1)):  
   f = "data %s:%s ../data/"%(i+1,i+1)+pilist1[i]+"\n"
   fakedata.insert(i+2,f)

fakedata.append("model apec*phabs\n")
for i in xrange(len(pilist1)):
    fakedata.append(T[i]+"\n")
    fakedata.append(A[i]+"\n")
    fakedata.append(z+"\n")
    fakedata.append("1.0\n")
    fakedata.append(phabs+"\n")
fakedata.append("fakeit\n")
fakedata.append("y\n\n")
for i in xrange(len(pilist1)):
    fakedata.append(prefix+"fake-"+Region[i]+"-"+Region[i+1]+".pi\n")
    fakedata.append("100000,1,100000\n")
fakedata.append("exit\n")
fakedata.append("<<EOF\n")
output = open("fake_step1_%s.csh"%prefix,"w")
for i in xrange(len(fakedata)):
    output.write(fakedata[i])
output.close()
 
fakedata2=["#!/bin/csh -f\n","xspec <<EOF\n","statistic chi\n","\n","setpl ene\n"]
for i in xrange(len(pilist1)):
   fakedata2.append("ignore %s:0.0-0.1,11.0-**\n"%(i+1))
fakedata2.append("ignore bad\n")
fakedata2.append("method leven 10 0.01\n")
fakedata2.append("abund angr\n")
fakedata2.append("xsect bcmc\n")
fakedata2.append("cosmo 70 0 0.7\n")
fakedata2.append("xset delta -1\n")
fakedata2.append("systematic 0\n")

for i in xrange(len(pilist1)):  
   f = "data %s:%s "%(i+1,i+1)+prefix+"fake-"+Region[i]+"-"+Region[i+1]+".pi\n"
   fakedata2.insert(i+2,f)
fakedata2.append("cpd /xw\n")
fakedata2.append("plot ldata\n")
fakedata2.append("sys rm fake-data_%s.qdp fake-data_%s.pco\n"%(prefix,prefix))
fakedata2.append("iplot \n")
fakedata2.append("we fake-data_%s\n"%prefix)
fakedata2.append("exit\n")
#fakedata2.append("exit\n")
fakedata2.append("sys sed -i -e '1,3d' fake-data_%s.qdp\n"%prefix)
fakedata2.append("<<EOF\n")
output = open("fake_step2_%s.csh"%prefix,"w")
for i in xrange(len(fakedata2)):
    output.write(fakedata2[i])
output.close()



