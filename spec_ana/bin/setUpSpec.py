#!/usr/bin/env python


import sys,string
import numpy as np;
import commands

mos1=commands.getoutput("head -1 log/prefix.log | gawk '{print $1}'")
mos2=commands.getoutput("head -2 log/prefix.log | tail -1| gawk '{print $1}'")
pn=commands.getoutput("tail -1 log/prefix.log | gawk '{print $1}'")

r=np.loadtxt("SetUpFile/Region.dat",dtype='d',usecols=(0,),unpack=True);


mos1output_name="spspec-mos1.txt"
mos2output_name="spspec-mos2.txt"
pnoutput_name="spspec-pn.txt"

mos1output=file(mos1output_name,"w");
mos2output=file(mos2output_name,"w");
pnoutput=file(pnoutput_name,"w");
blnk=[""];

for i in np.arange(r.size-1) :
      mos11=["%s-sp-%.0lf-%.0lf.fits" % (mos1,r[i],r[i+1])]
      mos12=["%s-obj-%.0lf-%.0lf.pi" % (mos1,r[i],r[i+1])]
      mos21=["%s-sp-%.0lf-%.0lf.fits" % (mos2,r[i],r[i+1])]
      mos22=["%s-obj-%.0lf-%.0lf.pi" % (mos2,r[i],r[i+1])]
      pn1=["%s-sp-%.0lf-%.0lf.fits" % (pn,r[i],r[i+1])]
      pn2=["%s-obj-%.0lf-%.0lf.pi" % (pn,r[i],r[i+1])]
      np.savetxt(mos1output,mos11,fmt='%s',newline=" ");
      np.savetxt(mos1output,blnk,fmt='%s');
      np.savetxt(mos1output,mos12,fmt='%s',newline=" ");
      np.savetxt(mos1output,blnk,fmt='%s');
      np.savetxt(mos2output,mos21,fmt='%s',newline=" ");
      np.savetxt(mos2output,blnk,fmt='%s');
      np.savetxt(mos2output,mos22,fmt='%s',newline=" ");
      np.savetxt(mos2output,blnk,fmt='%s');
      np.savetxt(pnoutput,pn1,fmt='%s',newline=" ");
      np.savetxt(pnoutput,blnk,fmt='%s');
      np.savetxt(pnoutput,pn2,fmt='%s',newline=" ");
      np.savetxt(pnoutput,blnk,fmt='%s');

mos1output.close();
mos2output.close();
pnoutput.close();

