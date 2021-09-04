#!/usr/bin/env python

from sys import*
import os

if len(argv) != 5 : 
   print ""
   print "INPUT ERROR in python";
   print "python make_fitcode.py z nH "
   print " z           - redshift "
   print " nH          - hydrogen density "
   print " sample.xcm   - sample script for xspec "
   print " mode   - 0:usual 1:crossarf mode "
   print ""
   print " output "
   print "     ./spectrum.xcm    - fit script for xspec "
   print ""
   print " e.g : python make_fitcode.py 0.182 0.0121 sample_10.xcm 0"
   exit();

z = float(argv[1])
nH = float(argv[2])
samplefile = argv[3]
mode = argv[4]

#read scal factor from log/mos1protonscale.log and log/mos2protonscale.log and log/pnprotonscale.log
file1 = open("log/mos1protonscale.log","r")
file2 = open("log/mos2protonscale.log","r")
file3 = open("log/pnprotonscale.log","r")
MOS1 = file1.readlines()
MOS2 = file2.readlines()
PN   = file3.readlines()

BL1,BL2,BL3 = [],[],[]
SL1,SL2,SL3 = [],[],[]
n1=MOS1.index("          Angle     Scaling\n") 
n2=len(MOS1)-3
for i in xrange(n2-n1):
    BL1.append(MOS1[i+n1+2].split()[1]+"\n")
    BL2.append(MOS2[i+n1+2].split()[1]+"\n")
    BL3.append(PN[i+n1+2].split()[1]+"\n")
    SL1.append(MOS1[i+n1+2].split()[2])
    SL2.append(MOS2[i+n1+2].split()[2])
    SL3.append(PN[i+n1+2].split()[2])

print "\n   Solid Angle"
print "  n mos1 mos2  pn"
for i in xrange(n2-n1):
   print " %02.0f %.2f %.2f %.2f"%(i+1,float(BL1[i].replace('\n',' ')),float(BL2[i].replace('\n',' ')),float(BL3[i].replace('\n',' ')))
print "\n Soft Proton Scaling"
print "  n mos1 mos2  pn"
for i in xrange(n2-n1):
   print " %02.0f %.2f %.2f %.2f"%(i+1,float(SL1[i].replace('\n',' ')),float(SL2[i].replace('\n',' ')),float(SL3[i].replace('\n',' ')))
print "from : log/mos1protonscale.log, log/mos2protonscale.log , log/pnprotonscale.log\n"

file1.close()
file2.close()
file3.close()

#write fit spectrum code
File = open("%s"%samplefile, "r")
specxcm = File.readlines()
RegionFile = open("SetUpFile/Region.dat", "r")
RList = RegionFile.readlines()
prefixFile = open("log/prefix.log","r")
prefixList = prefixFile.readlines()

for i in xrange(len(RList)-1):
   num = 1+3*i 
   specxcm.insert(num,"ignore %s:0.0-0.4,11.0-**\n"%(num+2))
   specxcm.insert(num,"ignore %s:0.0-0.3,11.0-**\n"%(num+1))
   specxcm.insert(num,"ignore %s:0.0-0.3,11.0-**\n"%num)

specxcm.insert(num+3,"ignore bad\n")
specxcm.insert(num+3,"statistic chi\n")
specxcm.insert(num+3,"method leven 10 0.01\n")
specxcm.insert(num+3,"abund angr\n")
specxcm.insert(num+3,"xsect bcmc\n")
specxcm.insert(num+3,"cosmo 70 0 0.7\n")
specxcm.insert(num+3,"xset delta -1\n")
specxcm.insert(num+3,"systematic 0\n")
specxcm.insert(num+3,"setpl ene\n")

#write cross arf script
if mode ==str(1):
   for i in xrange(len(RList)-1):
      if i==0 : 
         num = 1+i
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+5+len(RList)*3-1,i+3,prefixList[2].split()[0],RList[i+2].replace('\n',''),\
                                                                RList[i+3].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+5+len(RList)*3-1,i+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+2+len(RList)*3-1,i+3,prefixList[2].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+2+len(RList)*3-1,i+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+4+len(RList)*3-1,i+2,prefixList[1].split()[0],RList[i+2].replace('\n',''),\
                                                                RList[i+3].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+4+len(RList)*3-1,i+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+1+len(RList)*3-1,i+2,prefixList[1].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+1+len(RList)*3-1,i+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+3+len(RList)*3-1,i+1,prefixList[0].split()[0],RList[i+2].replace('\n',''),\
                                                                RList[i+3].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+3+len(RList)*3-1,i+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+len(RList)*3-1,i+1,prefixList[0].split()[0],RList[i+1].replace('\n',''), \
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+len(RList)*3-1,i+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
      if i==1:
         num = 1+12*i
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+10+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+10+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+7+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+7+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+9+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+9+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+6+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+6+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+8+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+8+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+5+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+5+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
      if i==2:
         num = 1+12*i
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+18+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+18+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+15+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+15+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+12+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i-2].replace('\n',''),\
                                                                RList[i-1].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+12+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+17+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+17+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+14+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+14+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+11+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i-2].replace('\n',''),\
                                                                RList[i-1].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+11+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+16+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+16+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+13+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+13+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+10+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i-2].replace('\n',''),\
                                                                RList[i-1].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+10+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
      if i==3:
         num = 7+12*i
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+29+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+29+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+26+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+26+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+23+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i-2].replace('\n',''),\
                                                                RList[i-1].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+23+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+20+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i-3].replace('\n',''),\
                                                                RList[i-2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+20+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+28+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+28+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+25+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+25+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+22+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i-2].replace('\n',''),\
                                                                RList[i-1].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+22+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+19+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i-3].replace('\n',''),\
                                                                RList[i-2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+19+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+27+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+27+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+24+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+24+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+21+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i-2].replace('\n',''),\
                                                                RList[i-1].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+21+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(i+18+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i-3].replace('\n',''),\
                                                                RList[i-2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(i+18+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
      if 3<i<len(RList)-2:
         num = 19+12*i
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%((i-2)*6+26+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%((i-2)*6+26+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%((i-2)*6+23+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%((i-2)*6+23+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%((i-2)*6+25+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%((i-2)*6+25+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%((i-2)*6+22+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%((i-2)*6+22+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%((i-2)*6+24+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i+1].replace('\n',''),\
                                                                RList[i+2].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%((i-2)*6+24+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%((i-2)*6+21+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%((i-2)*6+21+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
      if i==len(RList)-2:
         num = 19+12*i
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(3+(i-3)*6+26+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(3+(i-3)*6+26+len(RList)*3-1,i*3+3,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(3+(i-3)*6+25+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(3+(i-3)*6+25+len(RList)*3-1,i*3+2,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"arf  %i:%i %s-%s-%s-%s-%s.arf \n"%(3+(i-3)*6+24+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i-1].replace('\n',''),\
                                                                RList[i].replace('\n',''),RList[i].replace('\n',''),RList[i+1].replace('\n','')))
         specxcm.insert(num,"resp %i:%i %s-%s-%s.rmf\n"%(3+(i-3)*6+24+len(RList)*3-1,i*3+1,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))

for i in xrange(len(RList)-1):
   num = 1+3*i
   specxcm.insert(num,"response  %i:%i pn-diag.rsp.gz\n"%(num+3,num+2))
   specxcm.insert(num,"response  %i:%i mos2-diag.rsp.gz\n"%(num+2,num+1))
   specxcm.insert(num,"response  %i:%i mos1-diag.rsp.gz\n"%(num+1,num))

for i in xrange(len(RList)-1):
   num = 1+3*i
   specxcm.insert(num,"data %i:%i %s-obj-os-%s-%s-grp_cross.pi\n"%(num+2,num+2,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
   specxcm.insert(num,"data %i:%i %s-obj-%s-%s-grp_cross.pi\n"%(num+1,num+1,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
   specxcm.insert(num,"data %i:%i %s-obj-%s-%s-grp_cross.pi\n"%(num,num,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))

specxcm.insert(num+3,"data %i:%i rass.pi\n"%(num+3,num+3))


NL1 = ["A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15","A16","A17","A18","A19","A20","A21","A22","A23","A24","A25","A26","A28","A29","A30"]
NL2 = ["B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12","B13","B14","B15","B16","B17","B18","B19","B20","B21","B22","B23","B24","B25","B26","B28","B29","B30"]
NL3 = ["C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C28","C29","C30"]
NL4 = ["D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15","D16","D17","D18","D19","D20","D21","D22","D23","D24","D25","D26","D28","D29","D30"]
NL5 = ["E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13","E14","E15","E16","E17","E18","E19","E20","E21","E22","E23","E24","E25","E26","E28","E29","E30"]
NL6 = ["F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","F13","F14","F15","F16","F17","F18","F19","F20","F21","F22","F23","F24","F25","F26","F28","F29","F30"]

C1 = []
C2 = []
C3 = []
B1  = []
B2  = []
B3  = []

#write colum density and redshift
numw = specxcm.index("W\n")
specxcm.pop(numw)
specxcm.insert(numw,"%s\n"%nH)
numz = specxcm.index("Z\n")
specxcm.pop(numz)
specxcm.insert(numz,"%s\n"%z)

#write soft proton powerlaw model 
for i in xrange(n2-n1):
    C1.append(specxcm.index(NL1[i]+"\n"))
    C2.append(specxcm.index(NL2[i]+"\n"))
    C3.append(specxcm.index(NL3[i]+"\n"))
for i in xrange(n2-n1-1):
    B1.append(specxcm.index(NL4[i]+"\n"))
    B2.append(specxcm.index(NL5[i]+"\n"))
    B3.append(specxcm.index(NL6[i]+"\n"))
for i in xrange(n2-n1):
    A = specxcm.pop(C1[i])
    specxcm.insert(C1[i],BL1[i])
    B = specxcm.pop(C2[i])
    specxcm.insert(C2[i],BL2[i])
    C = specxcm.pop(C3[i])
    specxcm.insert(C3[i],BL3[i])
for i in xrange(n2-n1-1):    
    D = specxcm.pop(B1[i])
    specxcm.insert(B1[i],"="+SL1[i+1]+"*mos1_1:2\n")
    E = specxcm.pop(B2[i])
    specxcm.insert(B2[i],"="+SL2[i+1]+"*mos2_1:2\n")
    F = specxcm.pop(B3[i])
    specxcm.insert(B3[i],"="+SL3[i+1]+"*pn_1:2\n")

#write cross arf model
if mode==str(0):
   endpoint=specxcm.index("endpoint\n")
   null=specxcm.pop(endpoint)
   for i in xrange(((len(RList)-6)*6+3+33)*7):
      null=specxcm.pop(endpoint)
if mode==str(1):
   endpoint=specxcm.index("endpoint\n")
   null=specxcm.pop(endpoint)
   for i in xrange((len(RList)-6)*6+3+33):
      specxcm.insert(endpoint+8*i,"model  %s:crossmodel%s constant*constant*apec*phabs\n"%(len(RList)*3-1+i,i+1+100))

File.close()

if mode==str(0):
   output = open("spectrum.xcm","w")
if mode==str(1):
   output = open("spectrum_cross.xcm","w")

for i in xrange(len(specxcm)):
    output.write(specxcm[i])
output.close()


