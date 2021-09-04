#!/usr/bin/env python

from sys import*
import os

if len(argv) != 4 : 
   print ""
   print "INPUT ERROR in python";
   print "python add-xcm.py z nH "
   print " z           - redshift "
   print " nH          - hydrogen density "
   print " sample.xcm   - sample script for xspec "
   print ""
   print " output "
   print "     ./spectrum.xcm    - fit script for xspec "
   print ""
   print " e.g : python add-xcm.py  z nH"
   exit();

z = float(argv[1])
wabs = float(argv[2])
sample = argv[3]

file1 = open("log/mos1protonscale.log","r")
file2 = open("log/mos2protonscale.log","r")
file3 = open("log/pnprotonscale.log","r")
MOS1 = file1.readlines()
MOS2 = file2.readlines()
PN   = file3.readlines()

BL1 = []
BL2 = []
BL3 = []
SL1 = []
SL2 = []
SL3 = []

n1=MOS1.index("          Angle     Scaling\n") 
n2=len(MOS1)-3

for i in xrange(n2-n1):
    BL1.append(MOS1[i+n1+2].split()[1]+"\n")
    BL2.append(MOS2[i+n1+2].split()[1]+"\n")
    BL3.append(PN[i+n1+2].split()[1]+"\n")
    SL1.append(MOS1[i+n1+2].split()[2])
    SL2.append(MOS2[i+n1+2].split()[2])
    SL3.append(PN[i+n1+2].split()[2])
    
file1.close()
file2.close()
file3.close()

File = open("%s"%sample, "r")
specxcm = File.readlines()
RegionFile = open("SetUpFile/Region.dat", "r")
RList = RegionFile.readlines()

for i in xrange(len(RList)-1):
   num = 1+3*i
   specxcm.insert(num,"ignore %s:0.0-0.4,11.0-**\n"%(num+2))
   specxcm.insert(num,"ignore %s:0.0-0.3,11.0-**\n"%(num+1))
   specxcm.insert(num,"ignore %s:0.0-0.3,11.0-**\n"%num)
specxcm.insert(1,"setpl ene\n")
for i in xrange(len(RList)-1):
   num = 1+3*i
   specxcm.insert(num,"response  %s:%s pn-diag.rsp.gz\n"%(num+3,num+2))
   specxcm.insert(num,"response  %s:%s mos2-diag.rsp.gz\n"%(num+2,num+1))
   specxcm.insert(num,"response  %s:%s mos1-diag.rsp.gz\n"%(num+1,num))

prefixFile = open("log/prefix.log","r")
prefixList = prefixFile.readlines()

for i in xrange(len(RList)-1):
   num = 1+3*i
   specxcm.insert(num,"data %s:%s %s-obj-os-%s-%s-grp.pi\n"%(num+2,num+2,prefixList[2].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
   specxcm.insert(num,"data %s:%s %s-obj-%s-%s-grp.pi\n"%(num+1,num+1,prefixList[1].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
   specxcm.insert(num,"data %s:%s %s-obj-%s-%s-grp.pi\n"%(num,num,prefixList[0].split()[0],RList[i].replace('\n',''),RList[i+1].replace('\n','')))
specxcm.insert(num+3,"data %s:%s rass.pi\n"%(num+3,num+3))

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

numw = specxcm.index("W\n")
specxcm.pop(numw)
specxcm.insert(numw,"%s\n"%wabs)
numz = specxcm.index("Z\n")
specxcm.pop(numz)
specxcm.insert(numz,"%s\n"%z)

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

File.close()

output = open("spectrum.xcm","w")
for i in xrange(len(specxcm)):
    output.write(specxcm[i])
output.close()


