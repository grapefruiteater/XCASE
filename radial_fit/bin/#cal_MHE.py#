import os,sys,string
import math;
import numpy as np;
import scipy as scipy;
from scipy.interpolate import UnivariateSpline
from scipy.optimize import fsolve
from scipy.optimize import bisect
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
import pyfits

kpctocm=3.085677581e21
Mpctocm=3.085677581e24
Mp=1.6726231e-24
Msun=1.989e33
Munit=1e14
G=6.6726e-8 
kB=1.60217646e-9
#keV->erg

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
n_factor=1+ratio*(1+ratio_HE) # n=n_factor*ne rho_factor=n_factor
#print rho_factor, n_factor
#--------------------------------------------------------------------------------
if len(sys.argv) != 3 : 
   print "INPUT ERROR in python";
   print "python cal_MHE.py cogname z"
   print " cogname     - cluster name "
   print " z           - redshift "
   print ""
   print " output "
   print "     MHEresult/MHE_prof.dat    - M_HE  profile"
   print "     MHEresult/Mgas_prof.dat   - M_gas profile"
   print "     MHEresult/fgas_prof.dat   - f_gas profile"
   print "     MHEresult/MHE_Delta.dat   - M_HE  at r_Delta"
   print "     MHEresult/Mgas_Delta.dat  - M_gas at r_Delta"
   print "     MHEresult/fgas_Delta.dat  - f_gas at r_Delta"
   print ""
   print " e.g : python cal_MHE.py ABELL1689 0.1832 "
   sys.exit();

cogname=sys.argv[1];
z=float(sys.argv[2]);
#--------------------------------------------------------------------------------
H0=70.
Om0=0.3; #Omega_m0
OL=0.7;  #Omega_L
cos=FlatLambdaCDM(H0, Om0, Neff=0)
amtokpc=cos.kpc_proper_per_arcmin(z).value
amtoMpc=amtokpc/1e3;
rho_cr=cos.critical_density(z).value
Omz=cos.Om(z);
#--------------------------------------------------------------------------------
if os.path.isdir("MHEresult") == False:
   os.mkdir("MHEresult")


ram1,r1,rho,rho_m,rho_p=np.loadtxt("nresult/rho.dat",dtype='d',usecols=(0,1,2,3,4),unpack=True,skiprows=1);
ram2,r2,ne,ne_m,ne_p=np.loadtxt("nresult/ne.dat",dtype='d',usecols=(0,1,2,3,4),unpack=True,skiprows=1);
ram3,r3,Mgas,Mgas_m,Mgas_p=np.loadtxt("nresult/Mgas.dat",dtype='d',usecols=(0,1,2,3,4),unpack=True,skiprows=1);
ram4,r4,T,T_m,T_p=np.loadtxt("Tresult/T3d_bestfit.dat",dtype='d',usecols=(0,1,2,3,4),unpack=True,skiprows=1);

if ram1.size != r4.size :
   print "the size of nresult/ne.dat != Tresult/T3d_bestfit.dat"
   exit;

func_rho = UnivariateSpline(r1, rho, k=3,s=0);
func_rho_m = UnivariateSpline(r1, rho_m, k=3,s=0);
func_rho_p = UnivariateSpline(r1, rho_p, k=3,s=0);
func_ne = UnivariateSpline(r2, ne, k=3,s=0);
func_ne_m = UnivariateSpline(r2, ne_m, k=3,s=0);
func_ne_p = UnivariateSpline(r2, ne_p, k=3,s=0);
func_Mgas = UnivariateSpline(r3, Mgas, k=3,s=0);
func_Mgas_m = UnivariateSpline(r3, Mgas_m, k=3,s=0);
func_Mgas_p = UnivariateSpline(r3, Mgas_p, k=3,s=0);
func_T =UnivariateSpline(r4, T, k=3,s=0);
func_T_m = UnivariateSpline(r4, T_m, k=3,s=0);
func_T_p = UnivariateSpline(r4, T_p, k=3,s=0);


func_dTdr=func_T.derivative();
func_dTdr_m=func_T_m.derivative();
func_dTdr_p=func_T_p.derivative();

func_dndr=func_ne.derivative();
func_dndr_m=func_ne_m.derivative();
func_dndr_p=func_ne_p.derivative();


#pressure
P=ne*T;
P_m=ne_m*T_m;
P_p=ne_p*T_p;
func_P=UnivariateSpline(r4, P, k=3,s=0);
func_P_m=UnivariateSpline(r4, P_m, k=3,s=0);
func_P_p=UnivariateSpline(r4, P_p, k=3,s=0);

func_dPdr=func_P.derivative();
func_dPdr_m=func_P_m.derivative();
func_dPdr_p=func_P_p.derivative();

MHE=np.zeros(r2.size);
MHE_p=np.zeros(r2.size);
MHE_m=np.zeros(r2.size);

fgas=np.zeros(r2.size);
fgas_p=np.zeros(r2.size);
fgas_m=np.zeros(r2.size);

MHE2=np.zeros(r2.size);
dMHE2=np.zeros(r2.size);
MHE2_p=np.zeros(r2.size);
MHE2_m=np.zeros(r2.size);

rvalid=r2[0];
ivalid=0;
for i in np.arange(r2.size) :
    # combined T_p and n_m --> large error
    MHE[i]=-func_T(r2[i])*r2[i]*kB*Mpctocm*(r2[i]*func_dndr(r2[i])/func_ne(r2[i])+r2[i]*func_dTdr(r2[i])/func_T(r2[i]))/G/Munit/Msun/Mp/mu;
    MHE_p[i]=-func_T_p(r2[i])*r2[i]*kB*Mpctocm*(r2[i]*func_dndr_m(r2[i])/func_ne_m(r2[i])+r2[i]*func_dTdr_p(r2[i])/func_T_p(r2[i]))/G/Munit/Msun/Mp/mu;
    MHE_m[i]=-func_T_m(r2[i])*r2[i]*kB*Mpctocm*(r2[i]*func_dndr_p(r2[i])/func_ne_p(r2[i])+r2[i]*func_dTdr_m(r2[i])/func_T_m(r2[i]))/G/Munit/Msun/Mp/mu;

    #MHE[i]=-3.68e13*func_T(r2[i])*r2[i]*(r2[i]*func_dndr(r2[i])/func_ne(r2[i])+r2[i]*func_dTdr(r2[i])/func_T(r2[i]))/Munit;
    #MHE_p[i]=-3.68e13*func_T_p(r2[i])*r2[i]*(r2[i]*func_dndr_m(r2[i])/func_ne_m(r2[i])+r2[i]*func_dTdr_p(r2[i])/func_T_p(r2[i]))/Munit;
    #MHE_m[i]=-3.68e13*func_T_m(r2[i])*r2[i]*(r2[i]*func_dndr_p(r2[i])/func_ne_p(r2[i])+r2[i]*func_dTdr_m(r2[i])/func_T_m(r2[i]))/Munit;

#    MHE[i]=-func_dPdr(r2[i])*r2[i]*r2[i]*kB*Mpctocm*rho_factor/func_rho(r2[i])/G/Munit/Msun;
#    MHE_p[i]=-func_dPdr_p(r2[i])*r2[i]*r2[i]*kB*Mpctocm*rho_factor/func_rho_p(r2[i])/G/Munit/Msun;
#    MHE_m[i]=-func_dPdr_m(r2[i])*r2[i]*r2[i]*kB*Mpctocm*rho_factor/func_rho_m(r2[i])/G/Munit/Msun;
    #MHE[i]=-func_dPdr(r2[i])*r2[i]*r2[i]*kB*Mpctocm/func_ne(r2[i])/G/Munit/Msun/mu/Mp;
    #MHE_p[i]=-func_dPdr_p(r2[i])*r2[i]*r2[i]*kB*Mpctocm/func_ne_p(r2[i])/G/Munit/Msun/mu/Mp;
   # MHE_m[i]=-func_dPdr_m(r2[i])*r2[i]*r2[i]*kB*Mpctocm/func_ne_m(r2[i])/G/Munit/Msun/mu/Mp;
#    print r2[i],func_ne(r2[i]),ne[i],func_dPdr((r2[i]+r2[i+1])*0.5),(P[i+1]-P[i])/(r2[i+1]-r2[i]);
#    print r2[i],(4.*np.pi/3.)*2500.*rho_cr*r2[i]**3*(Mpctocm**3/Munit/Msun)-MHE[i]
    #print r2[i],r2[i]*func_dndr(r2[i])/func_ne(r2[i]),r2[i]*func_dTdr(r2[i])/func_T(r2[i]);

    fgas[i]=func_Mgas(r2[i])/MHE[i];
    fgas_p[i]=func_Mgas_p(r2[i])/MHE_m[i];
    fgas_m[i]=func_Mgas_m(r2[i])/MHE_p[i];
    print MHE[i],MHE_m[i],MHE_p[i]
    if MHE[i]< 0 or MHE_m[i]<0 or MHE_p[i]<0 :
       rvalid=r2[i];       
       ivalid=i;

#rvalid=0.;

func_MHE=UnivariateSpline(r2[r2>rvalid], MHE[r2>rvalid], k=3,s=0);
func_MHE_m=UnivariateSpline(r2[r2>rvalid], MHE_m[r2>rvalid], k=3,s=0);
func_MHE_p=UnivariateSpline(r2[r2>rvalid], MHE_p[r2>rvalid], k=3,s=0);

func_Mgas=UnivariateSpline(r2[r2>rvalid], Mgas[r2>rvalid], k=3,s=0);
func_Mgas_m=UnivariateSpline(r2[r2>rvalid], Mgas_m[r2>rvalid], k=3,s=0);
func_Mgas_p=UnivariateSpline(r2[r2>rvalid], Mgas_p[r2>rvalid], k=3,s=0);

func_fgas=UnivariateSpline(r2[r2>rvalid], fgas[r2>rvalid], k=3,s=0);
func_fgas_m=UnivariateSpline(r2[r2>rvalid], fgas_m[r2>rvalid], k=3,s=0);
func_fgas_p=UnivariateSpline(r2[r2>rvalid], fgas_p[r2>rvalid], k=3,s=0);


#--------------------------------------------------------
# estimate M_delta, r_delta

Delta_180m=180.*Omz;
Delta_200m=200*Omz;
Delta_vir=18*np.pi*np.pi*(1+0.4093*(1./Omz-1.)**0.90524)*Omz
Delta=np.array([2500,1000,500,200,Delta_vir,Delta_200m,Delta_180m])
Delta_char=np.array(["2500","1000","500","200","vir","200m","180m"])

r_Delta=np.ones(Delta.size);
r_Delta_m=np.ones(Delta.size);
r_Delta_p=np.ones(Delta.size);
MHE_Delta=np.ones(Delta.size);
MHE_Delta_m=np.ones(Delta.size);
MHE_Delta_p=np.ones(Delta.size);
Mgas_Delta=np.ones(Delta.size);
Mgas_Delta_m=np.ones(Delta.size);
Mgas_Delta_p=np.ones(Delta.size);
fgas_Delta=np.ones(Delta.size);
fgas_Delta_m=np.ones(Delta.size);
fgas_Delta_p=np.ones(Delta.size);

check=np.zeros(Delta.size);
check_m=np.zeros(Delta.size);
check_p=np.zeros(Delta.size);

r_initial_guess1=0.001
r_initial_guess2=1.
r_Delta_failed=-999;
M_Delta_failed=-999;
r_lowest=1e-4;

def deviation(r,delta) :
    return (4.*np.pi/3.)*delta*rho_cr*r**3*(Mpctocm**3/Munit/Msun)-func_MHE(r)

def deviation_m(r,delta) :
    return (4.*np.pi/3.)*delta*rho_cr*r**3*(Mpctocm**3/Munit/Msun)-func_MHE_m(r)

def deviation_p(r,delta) :
    return (4.*np.pi/3.)*delta*rho_cr*r**3*(Mpctocm**3/Munit/Msun)-func_MHE_p(r)


#for i2 in np.arange(Delta.size) :
#for i in np.arange(1) :
for i in np.arange(Delta.size) :
    #i=1
    delta=Delta[i];
    Solve_MHE=lambda r_delta :   (4.*np.pi/3.)*delta*rho_cr*r_delta**3*(Mpctocm**3/Munit/Msun)-func_MHE(r_delta)
    Solve_MHE_m=lambda r_delta : (4.*np.pi/3.)*delta*rho_cr*r_delta**3*(Mpctocm**3/Munit/Msun)-func_MHE_m(r_delta)
    Solve_MHE_p=lambda r_delta : (4.*np.pi/3.)*delta*rho_cr*r_delta**3*(Mpctocm**3/Munit/Msun)-func_MHE_p(r_delta)
    
    if deviation(r_initial_guess1,Delta[i])*deviation(r_initial_guess2,Delta[i]) > 0.:
       check[i]=1;
       continue;
    else :
       r_Delta[i] = bisect(Solve_MHE, r_initial_guess1,r_initial_guess2)
    if r_Delta[i] < 1e-4 :
       check[i]=1; # failed
       r_Delta[i]=r_Delta_failed;
    else :
       check[i]=0; 
    if check[i]==0 :
       r_Delta_m[i] = bisect(Solve_MHE_m, r_Delta[i]*0.6,r_Delta[i]*1.2)
       #r_Delta_m[i] = fsolve(Solve_MHE_m, r_Delta[i])
    else :
       r_Delta_m[i] = bisect(Solve_MHE_m, r_initial_guess1,r_initial_guess2)
       #r_Delta_m[i] = fsolve(Solve_MHE_m, r_initial_guess1)
    if r_Delta_m[i] < 1e-4 :
       check_m[i]=1; # failed
       r_Delta_m[i]=r_Delta_failed;
    if check[i]==0 :
       r_Delta_p[i] = bisect(Solve_MHE_p, r_Delta[i]*0.8,r_Delta[i]*1.4)
       #r_Delta_p[i] = fsolve(Solve_MHE_p, r_Delta[i]);
    else : 
       #r_Delta_p[i] = fsolve(Solve_MHE_m, r_initial_guess1);
       r_Delta_p[i] = bisect(Solve_MHE_m, r_initial_guess1,r_initial_guess2)
    if r_Delta_p[i] < 1e-4 :
       check_p[i]=1; # failed
       r_Delta_p[i]=r_Delta_failed;

    #print r_Delta[i], (4.*np.pi/3.)*delta*rho_cr*r_Delta[i]**3*(Mpctocm**3/Munit/Msun), func_MHE(r_Delta[i])
    if check[i]==0 :
       MHE_Delta[i]=func_MHE(r_Delta[i]);
       Mgas_Delta[i]=func_Mgas(r_Delta[i]);
       fgas_Delta[i]=func_fgas(r_Delta[i]);
    else  :
       MHE_Delta[i]=M_Delta_failed
       Mgas_Delta[i]=M_Delta_failed
       fgas_Delta[i]=M_Delta_failed
    if check_m[i]==0 :
       MHE_Delta_m[i]=func_MHE_m(r_Delta_m[i]);
       Mgas_Delta_m[i]=func_Mgas_m(r_Delta_m[i]);
       fgas_Delta_m[i]=func_fgas_m(r_Delta_m[i]);    
    else  :
       MHE_Delta_m[i]=M_Delta_failed
       Mgas_Delta_m[i]=M_Delta_failed
       fgas_Delta_m[i]=M_Delta_failed
    if check_p[i]==0 :
       MHE_Delta_p[i]=func_MHE_p(r_Delta_p[i]);
       Mgas_Delta_p[i]=func_Mgas_p(r_Delta_p[i]);
       fgas_Delta_p[i]=func_fgas_p(r_Delta_p[i]);
    else  :
       MHE_Delta_p[i]=M_Delta_failed
       Mgas_Delta_p[i]=M_Delta_failed
       fgas_Delta_p[i]=M_Delta_failed
 
    if check[i]==0:
       r_initial_guess1=r_Delta[i]
       r_initial_guess2=r_Delta[i]*3.
    elif check[i]==1 and check_p[i]==0 :
       r_initial_guess1=r_Delta_p[i]
       r_initial_guess2=r_Delta_p[i]*3.
    elif check[i]==1 and check_m[i]==0 :
       r_initial_guess1=r_Delta_m[i]
       r_initial_guess2=r_Delta_m[i]*3.
    #print "ooo",check[i],r_Delta[i],r_initial_guess1,r_initial_guess2

#----------------------------------------------------------------------------
output_MHE_prof="MHEresult/MHE_profile.dat"
output_Mgas_prof="MHEresult/Mgas_profile.dat"
output_fgas_prof="MHEresult/fgas_profile.dat"
output_MHE_Delta="MHEresult/MHE_Delta.dat"
output_Mgas_Delta="MHEresult/Mgas_Delta.dat"
output_fgas_Delta="MHEresult/fgas_Delta.dat"
output_MHE_log="MHEresult/MHE.log"

file_MHE_log=file(output_MHE_log,"w");
label="Cluster   : %s" % cogname 
np.savetxt(file_MHE_log, [label], fmt='%s');
label="cosmology : H0  = %.2lf [km/s/Mpc]" % H0
np.savetxt(file_MHE_log, [label], fmt='%s');
label="cosmology : Om0 = %.2lf " % Om0
np.savetxt(file_MHE_log, [label], fmt='%s');
label="cosmology : OL = %.2lf " % OL
np.savetxt(file_MHE_log, [label], fmt='%s');
label="z         : %lf " % z
np.savetxt(file_MHE_log, [label], fmt='%s');
label="1 arcmin = %.8e [kpc]" % amtokpc
np.savetxt(file_MHE_log, [label], fmt='%s');
label="1 arcmin = %.8e [Mpc]" % amtoMpc
np.savetxt(file_MHE_log, [label], fmt='%s');
file_MHE_log.close();

file_MHE_prof=file(output_MHE_prof,"w");
file_Mgas_prof=file(output_Mgas_prof,"w");
file_fgas_prof=file(output_fgas_prof,"w");

file_MHE_Delta=file(output_MHE_Delta,"w");
file_Mgas_Delta=file(output_Mgas_Delta,"w");
file_fgas_Delta=file(output_fgas_Delta,"w");


label="# r[Mpch_%d^{-1}] MHE[10^%dMsunh_%d^{-1}] MHE_m MHE_p" % (H0,math.log10(Munit),H0) 
line2=np.array([r2[r2>rvalid], MHE[r2>rvalid], MHE_m[r2>rvalid], MHE_p[r2>rvalid]]);
np.savetxt(file_MHE_prof, [label], fmt='%s');
np.savetxt(file_MHE_prof, line2.T, fmt='%0.8e',delimiter=' ', newline='\n');
file_MHE_prof.close();

label="# r[Mpch_%d^{-1}] Mgas[10^%dMsunh_%d^{-5/2}] Mgas_m Mgas_p" % (H0,math.log10(Munit),H0) 
line2=np.array([r2[r2>rvalid], Mgas[r2>rvalid], Mgas_m[r2>rvalid], Mgas_p[r2>rvalid]]);
np.savetxt(file_Mgas_prof, [label], fmt='%s');
np.savetxt(file_Mgas_prof, line2.T, fmt='%0.8e',delimiter=' ', newline='\n');
file_Mgas_prof.close();

label="# r[Mpch_%d^{-1}] fgas[h_%d^{-3/2}] fgas_m fgas_p" % (H0,H0) 
line2=np.array([r2[r2>rvalid], fgas[r2>rvalid], fgas_m[r2>rvalid], fgas_p[r2>rvalid]]);
np.savetxt(file_fgas_prof, [label], fmt='%s');
np.savetxt(file_fgas_prof, line2.T, fmt='%0.8e',delimiter=' ', newline='\n');
file_fgas_prof.close();


label="# Delta Delta(value) r[Mpch_%d^{-1}] r_m r_p ram[arcmin] ram_m ram_p MHE[10^%dMsunh_%d^{-1}] MHE_m MHE_p" % (H0,math.log10(Munit),H0) 
np.savetxt(file_MHE_Delta, [label], fmt='%s');
blnk=[""];
for i in np.arange(Delta.size):
    if check[i]==0 :
       line1=np.array([Delta_char[i]]);
       line2=np.array([Delta[i], r_Delta[i], r_Delta_m[i], r_Delta_p[i], r_Delta[i]/amtoMpc, r_Delta_m[i]/amtoMpc, r_Delta_p[i]/amtoMpc, MHE_Delta[i], MHE_Delta_m[i], MHE_Delta_p[i]]);
       np.savetxt(file_MHE_Delta,line1,fmt='%s',newline=" ");
       np.savetxt(file_MHE_Delta,line2,fmt='%.8e',newline=" ");
       np.savetxt(file_MHE_Delta,blnk,fmt='%s');
file_MHE_Delta.close();

label="# Delta Delta(value) r[Mpch_%d^{-1}] r_m r_p ram[arcmin] ram_m ram_p Mgas[10^%dMsunh_%d^{-5/2}] Mgas_m Mgas_p" % (H0,math.log10(Munit),H0) 
np.savetxt(file_Mgas_Delta, [label], fmt='%s');
blnk=[""];
for i in np.arange(Delta.size):
    if check[i]==0 :
       line1=np.array([Delta_char[i]]);
       line2=np.array([Delta[i], r_Delta[i], r_Delta_m[i], r_Delta_p[i], r_Delta[i]/amtoMpc, r_Delta_m[i]/amtoMpc, r_Delta_p[i]/amtoMpc,  Mgas_Delta[i], Mgas_Delta_m[i], Mgas_Delta_p[i]]);
       np.savetxt(file_Mgas_Delta,line1,fmt='%s',newline=" ");
       np.savetxt(file_Mgas_Delta,line2,fmt='%.8e',newline=" ");
       np.savetxt(file_Mgas_Delta,blnk,fmt='%s');
file_Mgas_Delta.close();

label="# Delta Delta(value) r[Mpch_%d^{-1}] r_m r_p ram[arcmin] ram_m ram_p fgas[h_%d^{-3/2}] fgas_m fgas_p" % (H0,H0) 
np.savetxt(file_fgas_Delta, [label], fmt='%s');
blnk=[""];
for i in np.arange(Delta.size):
    if check[i]==0 :
       line1=np.array([Delta_char[i]]);
       line2=np.array([Delta[i], r_Delta[i], r_Delta_m[i], r_Delta_p[i], r_Delta[i]/amtoMpc, r_Delta_m[i]/amtoMpc, r_Delta_p[i]/amtoMpc, fgas_Delta[i], fgas_Delta_m[i], fgas_Delta_p[i]]);
       np.savetxt(file_fgas_Delta,line1,fmt='%s',newline=" ");
       np.savetxt(file_fgas_Delta,line2,fmt='%.8e',newline=" ");
       np.savetxt(file_fgas_Delta,blnk,fmt='%s');
file_fgas_Delta.close();

#----------------------------------------------------------------------------
