import os,sys,string
import math;
import numpy as np;
import scipy as scipy;
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
import pyfits

kpctocm=3.085677581e21
Mpctocm=3.085677581e24
Mp=1.6726231e-24
Msun=1.989e33
Munit=1e14

#---------------
# n_HE=(1/10)n_H
#mu=0.62
#ratio=1./1.2
#rho_factor=1.92   #rho_g= rho_factor mu M_p n_e
#--------------
# Z=0.2
mu=0.5964
ratio=1./1.174417  # ratio=n_H/n_e
ratio_HE=0.0872083 # ratio_HE=n_HE/n_H
rho_factor=1.92574264507411  #rho_g= rho_factor mu M_p n_e
n_factor=1+ratio*(1+ratio_HE) # n=n_factor*ne
#---------------------------------
H0=70.
Om0=0.3;
OL=0.7;
cos=FlatLambdaCDM(H0, Om0, Neff=0)

if len(sys.argv) != 3 : 
   print "INPUT ERROR in python";
   print "python cal_ne0.py z image"
   print " z           - redshift "
   print " image       - count_rate image"
   print ""
   print " output "
   print "     nresult/ne.dat    - electron density profile"
   print "     nresult/rho.dat   - mass density profile"
   print "     nresult/Mgas.dat  - Gas mass profile"
   print ""
   print " e.g : python cal_ne0.py 0.1832 rate-400-2300-mos1-sp.fits"
   sys.exit();

z=float(sys.argv[1]);
dataname=sys.argv[2]

fits=pyfits.open(dataname) 
header=fits[0].header
CDELT1=header['CDELT1'];
CDELT2=header['CDELT2'];
pixtoas=math.sqrt(math.fabs(CDELT1*CDELT2))*3600.;
pixtoam=pixtoas/60.
pixtodeg=pixtoas/3600.


if os.path.isdir("nresult") == False:
   os.mkdir("nresult")

output_name="nresult/ne.dat"
output_name2="nresult/Mgas.dat"
output_name3="nresult/rho.dat"
outputfile=file(output_name,"w");
outputfile2=file(output_name2,"w");
outputfile3=file(output_name3,"w");
#---------------------------------------------------------
amtokpc=cos.kpc_proper_per_arcmin(z).value
dA=cos.angular_diameter_distance(z).value
degtocm=3600.*amtokpc*kpctocm
cmtodeg=1./degtocm

ram,n,n_m,n_p,n2,n2_m,n2_p=np.loadtxt("Sxresult/ne_bestfit.dat",dtype='d',usecols=(1,2,3,4,5,6,7),unpack=True,skiprows=1);

A=(4*np.pi*(1+z)**2*(180./np.pi)**3/1e-14/pixtodeg/(dA*Mpctocm))**0.5
A_m=A
A_p=A

rMpc=np.zeros(n.size);
nbest=np.zeros(n.size); #electron density
nbest_m=np.zeros(n.size);
nbest_p=np.zeros(n.size);
rhobest=np.zeros(n.size);
rhobest_m=np.zeros(n.size);
rhobest_p=np.zeros(n.size);

for i in np.arange(ram.size) :
    rMpc[i]=ram[i]*amtokpc/1e3;
    nbest[i]=A*n[i]/math.sqrt(ratio);
    nbest_m[i]=A_m*n_m[i]/math.sqrt(ratio);
    nbest_p[i]=A_p*n_p[i]/math.sqrt(ratio);
    rhobest[i]=rho_factor*mu*Mp*nbest[i];
    rhobest_m[i]=rho_factor*mu*Mp*nbest_m[i];
    rhobest_p[i]=rho_factor*mu*Mp*nbest_p[i];

#output
line=np.array([ram,rMpc,nbest,nbest_m,nbest_p]);
label="# r[arcmin] r[Mpch_%d^{-1}] ne[cm^-3h_%d^{1/2}] ne_m ne_p" % (H0,H0)
np.savetxt(outputfile, [label], fmt='%s');
np.savetxt(outputfile, line.T, fmt='%0.8e',delimiter=' ', newline='\n');

line=np.array([ram,rMpc,rhobest,rhobest_m,rhobest_p]);
label="# r[arcmin] r[Mpch_%d^{-1}] rho[gcm^-3h_%d^{1/2}] rho_m rho_p" % (H0,H0)
np.savetxt(outputfile3, [label], fmt='%s');
np.savetxt(outputfile3, line.T, fmt='%0.8e',delimiter=' ', newline='\n');

#------------------------------------------------------------------------------
# Gas Mass
Mgas=np.zeros(n.size);
Mgas_m=np.zeros(n.size);
Mgas_p=np.zeros(n.size);

func_Mgas = interp1d(ram, 4*np.pi*rhobest*ram*ram, kind='cubic');
func_Mgas_p = interp1d(ram, 4*np.pi*rhobest_p*ram*ram, kind='cubic');
func_Mgas_m = interp1d(ram, 4*np.pi*rhobest_m*ram*ram, kind='cubic');

for i in np.arange(n.size) :
    Mgas[i]=integrate.quad(lambda x: func_Mgas(x), ram[0], ram[i])[0]*(amtokpc*kpctocm)**3/Munit/Msun
    Mgas_p[i]=integrate.quad(lambda x: func_Mgas_p(x), ram[0], ram[i])[0]*(amtokpc*kpctocm)**3/Munit/Msun
    Mgas_m[i]=integrate.quad(lambda x: func_Mgas_m(x), ram[0], ram[i])[0]*(amtokpc*kpctocm)**3/Munit/Msun


line2=np.array([ram,rMpc,Mgas,Mgas_m,Mgas_p]);
label="# r[arcmin] r[Mpch_%d^{-1}] Mgas[10^%dMsunh_%d^{-5/2}] Mgas_m Mgas_p" % (H0,math.log10(Munit),H0) 
np.savetxt(outputfile2, [label], fmt='%s');
np.savetxt(outputfile2, line2.T, fmt='%0.8e',delimiter=' ', newline='\n');





