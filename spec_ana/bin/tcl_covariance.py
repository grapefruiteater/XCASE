#!/usr/bin/env python

from sys import *
import os
import commands,math,itertools
import numpy as np
import array
import pyfits
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.ticker as ptick
import pylab
from matplotlib import rcParams

fparams = {'axes.labelsize': 25,
          'axes.titlesize': 25,
          'legend.fontsize': 20,
          'xtick.labelsize': 20,
          'ytick.labelsize': 20,
          'axes.linewidth' : 2,
          'xtick.major.size' : 5,
          'xtick.major.width' : 2,
          'xtick.minor.size' : 2,
          'xtick.minor.width' : 2,
          'ytick.major.size' : 5,
          'ytick.major.width' : 2,
          'ytick.minor.size' : 2,
          'ytick.minor.width' : 2,
         'font.family': 'serif'}
rcParams.update(fparams)

#-------------------------------------------------
mos1=commands.getoutput("head -1 log/prefix.log | gawk '{print $1}'")
mos2=commands.getoutput("head -2 log/prefix.log | tail -1 | gawk '{print $1}'")
pn=commands.getoutput("tail -1 log/prefix.log  | gawk '{print $1}'")
#-------------------------------------------------

if len(argv) != 2 :
    print ""
    print "Usage:", argv[0], "clustername"
    print "";
    exit()

ClusterName = argv[1]

regionList=open("../spec_ana_crossarf/SetUpFile/Region.dat")
region=regionList.readlines()

### output xspec/fit log file for calucurate free parameter###
if os.path.exists("./xspec_fit_result.dat")==True:
    logspecfile=open("xspec_fit_result.dat","r")
    logdata=logspecfile.readlines()
    index=logdata.index("# par  comp\n")
else:
    os.system("csh bin/tclout.csh")
    logspecfile=open("xspec_fit_result.dat","r")
    logdata=logspecfile.readlines()
    index=logdata.index("# par  comp\n")

num=0
tmp=0
freetn=array.array("f")
for i in xrange(index,len(logdata)):
    if len(logdata[i].split()) >= 7:
        if logdata[i].split()[1]!="Data group: 1\n":
            if logdata[i].split()[-5]=="kT":
                freetn.append(tmp)
                print tmp,logdata[i].split()[-6],logdata[i].split()[-3],logdata[i].split()[-2],logdata[i].split()[-1]
            tmp=tmp+1
        num=num+1
error=np.zeros(num)

num=0
for i in xrange(index,len(logdata)):
    if len(logdata[i].split()) >= 7:
        error[num]=float(logdata[i].split()[-1])
        num=num+1

### -----------write num to tcl_covariance.xcm------------- ###
covxcmfile=open("/cluster491/users/miyaoka/tcl_covariance.xcm","r")
covxcm=covxcmfile.readlines()
covxcm[6]="set num %i \n"%num
os.system("rm %s-tcl_covariance.xcm"%ClusterName)
os.system("rm tcl_covariance.txt")
output = open("%s-tcl_covariance.xcm"%ClusterName,"w")
for k in xrange(len(covxcm)):
    output.write(covxcm[k])
output.close()

os.system("xspec %s-tcl_covariance.xcm"%ClusterName)
os.system("cp tcl_covariance.txt tcl_covariance.txt.log")

if os.path.exists("./tcl_covariance.txt")==True:
    print "\n####################################################"
    print "######SAVE ALL COVARIANCE to tcl_covariance.txt#####"
    print "####################################################"

savespecfile=open("tcl_covariance.txt","r")
savedata=savespecfile.readlines()
freepn=int(len(savedata)**0.5)
covariance_matrix=np.zeros((freepn,freepn))
for i in xrange(freepn):
    for j in xrange(freepn):
        covariance_matrix[i,j]=savedata[j+i*freepn]
np.savetxt("%s-covariance_matrix.txt"%ClusterName, covariance_matrix, fmt="%0.20f",delimiter=" ")

print "\n#################################"
print "######SAVE COVARIANCE MATRIX#####"
print "#################################"

######### compare between error with covariance #########
saveerrorfile=open("save-error_cross.xcm","r")
saveerror=saveerrorfile.readlines()
for i in xrange(len(saveerror)/3):
    tmp1=saveerror[3*i+2].split()[3]
    tmp2=tmp1.split(",")
    sigmamax=math.fabs(float(tmp2[0].replace("(","")))
    sigmamin=math.fabs(float(tmp2[1].replace(")","")))
    exec("sigma%i=((sigmamax**2+sigmamin**2)/2.0)**0.5"%(i+1))

##########CHANGE THIE LIST################
tempnum=[10,25,38,51,64,77]#typical number
##########################################

print "\n######################################################################"
print "###### CHECK THE ROUGH CONSISTENCY 1SIGMA_ERROR**2 WITH SIGMA**2 #####"
print "######################################################################"
print "\n (num1,num2) 1sigmaerror^2 sigma^2 sigma1*sigma2 coefficient"
####1sigmaerror is average error derived by not using XSPEC/error task.
####sigma is average 1_sigma error from XSPEC/error task.

matrix_length=int(len(tempnum))
net_covariance_matrix=np.zeros((matrix_length,matrix_length))
coefficient_matrix=np.ones((matrix_length,matrix_length))
combi=list(itertools.combinations(tempnum,2))
for i in xrange(len(combi)):
    diag1=combi[i][0]
    diag2=combi[i][1]
    index1=tempnum.index(diag1)
    index2=tempnum.index(diag2)
    print combi[i],
    exec("tmp1=sigma%i**2"%(index1+1))
    gamma=covariance_matrix[combi[i]]/math.fabs(covariance_matrix[diag1,diag1])**0.5/math.fabs(covariance_matrix[diag2,diag2])**0.5
    coefficient_matrix[index1,index2],coefficient_matrix[index2,index1]=gamma,gamma
    exec("tmp2=sigma%i*sigma%i"%(index1+1,index2+1))
    net_covariance_matrix[index1,index2],net_covariance_matrix[index2,index1]=tmp2*gamma,tmp2*gamma
    net_covariance_matrix[index1,index1]=tmp1
    print "%0.4f %0.4f %0.4f %0.4f"%(tmp1,covariance_matrix[diag1,diag1],covariance_matrix[combi[i]],gamma)

exec("net_covariance_matrix[-1,-1]=sigma%i**2"%len(tempnum))
print "\n ### Gamma coefficient matrix"
print coefficient_matrix
print "\n ### Calcurated new covariance matrix"
print net_covariance_matrix
print "\n ### Calcurated inverse covariance matrix"
print  np.linalg.inv(net_covariance_matrix)
#---OUTPUT---
np.savetxt("../result/net-covariance_matrix.dat" , net_covariance_matrix, fmt="%0.20f",delimiter=" ")
np.savetxt("../result/inverse-covariance_matrix.dat" , np.linalg.inv(net_covariance_matrix), fmt="%0.20f",delimiter=" ")

print "\n check the Gamma coefficient in the range of -1<gamma<1 "  
combi=list(itertools.combinations(range(freepn),2))
for i in xrange(len(combi)):
    diag1=combi[i][0]
    diag2=combi[i][1]
    tmp=diag1 in tempnum
    if -1.0<covariance_matrix[combi[i]]/covariance_matrix[diag1,diag1]/covariance_matrix[diag2,diag2]<1.0 and tmp==True:
        print combi[i],"%0.4f %0.4f %0.4f"%(covariance_matrix[diag1,diag1],covariance_matrix[combi[i]],covariance_matrix[combi[i]]/math.fabs(covariance_matrix[diag1,diag1])**0.5/math.fabs(covariance_matrix[diag2,diag2])**0.5)
