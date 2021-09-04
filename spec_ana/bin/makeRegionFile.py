#!/usr/bin/env python

from sys import *
import numpy as np
import pyfits
import commands

mos1=commands.getoutput("head -1 log/prefix.log | gawk '{print $1}'")

inputfits=mos1+"-bkg_region-sky.fits"
output_name=mos1+"-bkg_region-sky.reg"

t = pyfits.open(inputfits)
tbdata = t[1].data 

shape=tbdata['SHAPE'];
x=tbdata['X'];
y=tbdata['Y'];
r=tbdata['R'];
rotang=tbdata['ROTANG'];
comp=tbdata['COMPONENT'];


output=file(output_name,"w");
blnk=[""];
np.savetxt(output,["# Region file format: DS9 version 4.1"],fmt='%s');
np.savetxt(output,["global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1"],fmt='%s')
np.savetxt(output,["detector"],fmt='%s')
 
for i in np.arange(comp.size):
    shape0=shape[i]
    shape0=np.core.defchararray.replace(shape0,"!","")
    shape0=np.core.defchararray.replace(shape0," ","")
    line=np.array([str(shape0)+"("+str(x[i][0])+","+str(y[i][0])+","+str(r[i][0])+")"]);
    np.savetxt(output,line,fmt='%s',newline=" ");
    np.savetxt(output,blnk,fmt='%s');
output.close();
