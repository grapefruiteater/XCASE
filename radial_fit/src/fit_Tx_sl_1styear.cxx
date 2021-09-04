//icpc  fit_Tx_sl_1styear.cxx  -I/usr/local/include/Minuit2- -I/usr/local/include/Math -L/usr/local/lib/ -lMinuit2  $CPGPLOT_FLAGS $GSL_FLAGS $MINUIT_FLAGS $CFITSIO_FLAGS -parallel -o "fit_Tx_sl_1styear"

#include<iostream>
#include<fstream>
#include<algorithm>
#include<cstdio>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include <vector>
#include"cpgplot.h"
#include"fitsio.h" 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/stat.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_sf_dilog.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_errno.h>
#include "spline.h"


#include<cstdio>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include <cassert>
#include<gsl/gsl_sf_gamma.h>

#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"

double PIXTOAS_PSFFIX;
double pixratio;

using namespace std;
const istreambuf_iterator<char>  EOF_MARKER;
using namespace ROOT::Minuit2;
  

#define TwoPi  2*M_PI
#define Pi M_PI

#define TINY 1.e-5
#define EPS 1.e-15
#define cv 2.99792458e10
#define ec 4.80298e-10
#define kB 1.60217646e-9
#define h  6.62606876e-27
#define Mp 1.6726231e-24
#define Me 9.1093897e-28
#define Ktoerg 1.380658e-16
#define KtokeV 8.617385e-8
#define Msun   1.989e33
#define G       6.6726e-8
#define yrtosec 3.155815e7
#define kpctocm 3.0856e21 // 1kpc=3.0856e21 cm
#define Mpctocm 3.0856e24 // 1kpc=3.0856e21 cm

#define MASK_BORDER 260
#define MASK_TH 1e4

double amtokpc;       // arcmin to kpc
double Sigma_cr;      // critical surface density
double delta_c;       // density contrast
double rho_cr;        // critical density
double zd;            // lens redshift
double zs;            // mean source redshift
double Hubble;        // Hubble constant 70 km/s/Mpc
double Om0;           // Omega Matter
double OL0;           // Omega Lambda
double Ob0;           // Omeba baryon
double E;             // E(z)=H(z)/Hubble
double Omz;           // Omega(z)=Om0*(1+zd)^3/E(z)^2
double OLz;           // Omega(z)=Om0*(1+zd)^3/E(z)^2
double Hz;            // H(z)=Hubble*E(z)
double M0;            // M0=4.*Pi*rho_cr*delta_c/3.
double Dratio;        // Distance ratio Dds/Ds      
double Ds;            // Angular Diameter Distance to source
double Dd;            // Angular Diameter Distance to lens
double Dds;           // Angular Diameter Distance from source to lens
double Dl;            // Luminosity Distance to lens
double hubble_norm;   // Hubble/100;
double sigma8;        // sigma8;

//double factor_rcool=1.; // rcool=factor_rcool*rt;
double Tmin_fixed=0;
double rcool_fixed=1e10;
double acool_fixed=0.;
double factor_rcool_fixed=0;
//double aT_fixed=0.;  
double bT_fixed=2.;  
double H0=70;

double RMAX; // maximun radius for numerical integration along the line-of-sight;

//=======================================================================
namespace ROOT {

   namespace Minuit2 {


class T3d {

public:
    T3d(double T0,   double Tmin,double factor_rcool,double rt, 
        double acool,double aT,  double bT,   double cT) :
      T0(T0), Tmin(Tmin),factor_rcool(factor_rcool), rt(rt),acool(acool), aT(aT), bT(bT),cT(cT){}
   ~T3d() {}

  double operator()(double r) const {
   
    double rnorm=r/rt;
    //double rcool=pow(10.,factor_rcool)*rt;
    //double xT=pow(r/rcool,acool);
   
    return T0*pow(rnorm,-aT)*pow(1+pow(rnorm,bT),-cT/bT);

  }

private:

  double T0;
  double Tmin;
  double factor_rcool;
  double rt;
  double acool;
  double aT;
  double bT;
  double cT;

};

class Tew {

public:
   
    Tew(double T0,   double Tmin,double factor_rcool,double rt, 
        double acool,double aT,  double bT,   double cT,const tk::spline& n2spline) :
      T0(T0), Tmin(Tmin),factor_rcool(factor_rcool), rt(rt),acool(acool), aT(aT), bT(bT),cT(cT),n2spline(n2spline){}
  ~Tew() {}


  double operator()(double r) const {


  //r/R=x = cosh(t)
    double tmin=0.; 
    double tmax=acosh(RMAX/r); //tmax=RMAX [arcmin]/r [arcmin]
    int Nline=50; 
    double dt=(tmax-tmin)/Nline;
    int iline;
    double wgt1=0.;
    double wgt2=0.;
    
    double sum1=0.;
    double sum2=0.;
    double ttmp;
    double rtmp;
    double rmin=cosh(tmin)*r;
    double rmax=cosh(tmax)*r; 

    T3d T3d(T0,Tmin,factor_rcool,rt,acool,aT,bT,cT);

    wgt1=n2spline(rmin)*pow(T3d(rmin),0.5);
    wgt2=n2spline(rmax)*pow(T3d(rmax),0.5);

    sum1=T3d(rmin)*wgt1*cosh(tmin)+T3d(rmax)*wgt2*cosh(tmax);
    sum2=wgt1*cosh(tmin)+wgt2*cosh(tmax);

    for(iline=1;iline<Nline;iline+=2){
      ttmp=tmin+iline*dt;
      rtmp=cosh(ttmp)*r;
      wgt1=n2spline(rtmp)*pow(T3d(rtmp),0.5);
      if(wgt1<0){wgt1=0;}
      sum1+=4*T3d(rtmp)*wgt1*cosh(ttmp);
      sum2+=4*wgt1*cosh(ttmp);

    }

    for(iline=2;iline<Nline-1;iline+=2){

      ttmp=tmin+iline*dt;
      rtmp=cosh(ttmp)*r;
      wgt1=n2spline(rtmp)*pow(T3d(rtmp),0.5);
      if(wgt1<0){wgt1=0;}
      sum1+=2*T3d(rtmp)*wgt1*cosh(ttmp);
      sum2+=2*wgt1*cosh(ttmp);
    }

    // sum1=2.*dt*sum1/3.;
    //  sum2=2.*dt*sum2/3.;
    
    return sum1/sum2;
    

  }

private:

  double T0;
  double Tmin;
  double factor_rcool;
  double rt;
  double acool;
  double aT;
  double bT;
  double cT;
  tk::spline n2spline;
};


class Tsl {

public:
   
    Tsl(double T0,   double Tmin,double factor_rcool,double rt, 
        double acool,double aT,  double bT,   double cT,const tk::spline& n2spline) :
      T0(T0), Tmin(Tmin),factor_rcool(factor_rcool), rt(rt),acool(acool), aT(aT), bT(bT),cT(cT),n2spline(n2spline){}
  ~Tsl() {}


  double operator()(double r) const {

    //linar grid ---> 1./sqrt(rtmp^2-r^2) "rtmp around rmin" --> infty
    //log grid is essenstial

    double tmin=0.; 
    double tmax=acosh(RMAX/r);
    int Nline=50; 
    double dt=(tmax-tmin)/Nline;
    int iline;
    double wgt1=0.;
    double wgt2=0.;
    
    double sum1=0.;
    double sum2=0.;
    double ttmp;
    double rtmp;
    double rmin=cosh(tmin)*r;
    double rmax=cosh(tmax)*r; 

    T3d T3d(T0,Tmin,factor_rcool,rt,acool,aT,bT,cT);

    wgt1=n2spline(rmin)*pow(T3d(rmin),-3./4);
    wgt2=n2spline(rmax)*pow(T3d(rmax),-3./4);

    sum1=T3d(rmin)*wgt1*cosh(tmin)+T3d(rmax)*wgt2*cosh(tmax);
    sum2=wgt1*cosh(tmin)+wgt2*cosh(tmax);

    for(iline=1;iline<Nline;iline+=2){
      ttmp=tmin+iline*dt;
      rtmp=cosh(ttmp)*r;
      wgt1=n2spline(rtmp)*pow(T3d(rtmp),-3./4);
      if(wgt1<0){wgt1=0;}
      sum1+=4*T3d(rtmp)*wgt1*cosh(ttmp);
      sum2+=4*wgt1*cosh(ttmp);

    }

    for(iline=2;iline<Nline-1;iline+=2){

      ttmp=tmin+iline*dt;
      rtmp=cosh(ttmp)*r;
      wgt1=n2spline(rtmp)*pow(T3d(rtmp),-3./4);
      if(wgt1<0){wgt1=0;}
      sum1+=2*T3d(rtmp)*wgt1*cosh(ttmp);
      sum2+=2*wgt1*cosh(ttmp);
    }

    //      sum1=2.*dt*sum1/3.;
    //  sum2=2.*dt*sum2/3.;
    
    return sum1/sum2;


  }

private:

  double T0;
  double Tmin;
  double factor_rcool;
  double rt;
  double acool;
  double aT;
  double bT;
  double cT;
  tk::spline n2spline;
};



class Chisq : public FCNBase {

public:
  Chisq(const std::vector<double>& rmid,
        const std::vector<double>& rinn,
        const std::vector<double>& rout,
        const std::vector<double>& Tx,
        const std::vector<double>& dTx,
        const tk::spline& n2spline) : 
    rmid(rmid), rinn(rinn), rout(rout), Tx(Tx), dTx(dTx),n2spline(n2spline),fErrorDef(1.) {}

  ~Chisq() {}
  virtual double Up() const {return fErrorDef;} 
  virtual double operator()(const std::vector<double>&) const;

  void SetErrorDef(double def) {fErrorDef = def;}

private:

  std::vector<double> rmid;
  std::vector<double> rinn;
  std::vector<double> rout;
  std::vector<double> Tx;
  std::vector<double> dTx;
  tk::spline n2spline;

  double fErrorDef;

};


double Chisq::operator()(const std::vector<double>& para) const {

  //Tsl gives the same answer as Tsl_area
  double chisq=0.;
  //  Tew Tew(para[0], para[1], para[2], para[3], para[4], aT_fixed, bT_fixed, para[5],n2spline);
  //  Tew Tew(para[0], Tmin_fixed, factor_rcool_fixed, para[1], acool_fixed, para[2], bT_fixed, para[3],n2spline);
  Tsl Tew(para[0], Tmin_fixed, factor_rcool_fixed, para[1], acool_fixed, para[2], bT_fixed, para[3],n2spline);



  int ibin;
  int Nbin=rmid.size();
  double Tmodel;

  for(ibin=0;ibin<Nbin;ibin++){
    Tmodel=Tew(rmid[ibin]);
       chisq+=(Tmodel-Tx[ibin])*(Tmodel-Tx[ibin])/dTx[ibin]/dTx[ibin];
 }
  return chisq;

}

//------------------------

class Chisq_3d : public FCNBase {
public:
  Chisq_3d(const std::vector<double>& rmid,
        const std::vector<double>& rinn,
        const std::vector<double>& rout,
        const std::vector<double>& Tx,
        const std::vector<double>& dTx,
        const tk::spline& n2spline) : 
    rmid(rmid), rinn(rinn), rout(rout), Tx(Tx), dTx(dTx),n2spline(n2spline),fErrorDef(1.) {}

  ~Chisq_3d() {}
  virtual double Up() const {return fErrorDef;} 
  virtual double operator()(const std::vector<double>&) const;

  void SetErrorDef(double def) {fErrorDef = def;}

private:

  std::vector<double> rmid;
  std::vector<double> rinn;
  std::vector<double> rout;
  std::vector<double> Tx;
  std::vector<double> dTx;
  tk::spline n2spline;

  double fErrorDef;

};


double Chisq_3d::operator()(const std::vector<double>& para) const {

  double chisq=0.;
  //  T3d T3d(para[0], para[1], para[2], para[3], para[4], aT_fixed, bT_fixed, para[5]);
  //T3d T3d(para[0], para[1], para[2], para[3], acool_fixed, aT_fixed, bT_fixed, para[4]);
  T3d T3d(para[0], Tmin_fixed, factor_rcool_fixed, para[1], acool_fixed, para[2], bT_fixed, para[3]); 

  int ibin;
  int Nbin=rmid.size();
  double Tmodel;

  for(ibin=0;ibin<Nbin;ibin++){
    Tmodel=T3d(rmid[ibin]);
    chisq+=(Tmodel-Tx[ibin])*(Tmodel-Tx[ibin])/dTx[ibin]/dTx[ibin];
 }
  return chisq;

}



  }  // namespace Minuit2

}  // namespace ROOT

//=======================================================================

int main(int argc, char *argv[])
{
  if(argc!=4){
    fprintf(stderr, "\n");
    fprintf(stderr, " Use : ./fit_Tx_ew temp.dat n.dat bT\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      temp.dat         - temperature data\n");
    fprintf(stderr, "      n.dat            - num density data\n");
    fprintf(stderr, "      acool            - fixed parameter in Vikhlinin model (acool=1~2 ?) : Recommendation acool=2\n");
    fprintf(stderr, "      aT               - fixed parameter in Vikhlinin model (aT=0)\n");
    fprintf(stderr, "      bT               - fixed parameter in Vikhlinin model (bT=2)\n");
    fprintf(stderr, "      model : T3d= T0*(r/rt)^-aT*(1+(r/rt)^bT)^(-cT/bT)\n");
    fprintf(stderr, "      output files in Tresult\n");
    fprintf(stderr, "      best-fit Tx profile  : Tx_bestfit.dat\n");
    fprintf(stderr, "      best-fit T3d profile : T3d_bestfit.dat\n");
    fprintf(stderr, "      log file             : Tx_bestfit.log\n");
    fprintf(stderr, "      \n");
    fprintf(stderr, "  COMMENT : \n");
    fprintf(stderr, "     Numerical integration along the line-of-sight is up to the maximam radius of n.dat (default 30 arcmin)\n");
    fprintf(stderr, "     If you analyze low-z clusters, change the maximum radius in n.dat\n");
    fprintf(stderr, "      \n");
    fprintf(stderr, "      \n");
    fprintf(stderr, " e.g. : ./fit_Tx_ew temp.dat nresult/ne.dat 2\n");
    fprintf(stderr, "\n");
    return 1;
  }

   struct stat st;
   if(stat("Tresult", &st)!=0){
     mkdir("Tresult",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   }

    char tempname[1000];
    char nename[1000];
    sprintf(tempname,"%s",argv[1]);
    sprintf(nename,"%s",argv[2]);
    bT_fixed=atof(argv[3]);
    /*
    acool_fixed=atof(argv[3]);
    aT_fixed=atof(argv[4]);
    */

    ifstream fin;
    fin.open(tempname);
   
    if(!fin){fprintf(stderr,"NO FILE : %s\n", tempname); return 1;}
    streampos current;

    //count number
    char colum1[1000];
    colum1[0]='#';
     do{
       current = fin.tellg();
       fin.getline(colum1,sizeof(colum1));
     }while(colum1[0]=='#');
    // 行数を数えて結果ファイルに書く
    fin.seekg(current);
    istreambuf_iterator<char>  inF(fin);
    int Nbin = count( inF, EOF_MARKER, '\n' );
    fin.seekg(current);
    std::vector<double> rmid(Nbin);
    std::vector<double> rinn(Nbin);
    std::vector<double> rout(Nbin);
    std::vector<double> Tx(Nbin);
    std::vector<double> Tx_m(Nbin);
    std::vector<double> Tx_p(Nbin);
    std::vector<double> dTx(Nbin);

    //setup initial parameters
    double T0_prior=0.;
    double Tmin_prior=0.;
    double rcool_prior=0.;
    double rt_prior=0.;
    double acool_prior=0.;
    double aT_prior=0.;
    double cT_prior=1.;

    int ibin;
    for(ibin=0;ibin<Nbin;ibin++){
      fin>>rmid[ibin]>>rinn[ibin]>>rout[ibin]>>Tx[ibin]>>Tx_m[ibin]>>Tx_p[ibin];
      dTx[ibin]=sqrt(0.5*((Tx_m[ibin]-Tx[ibin])*(Tx_m[ibin]-Tx[ibin])+(Tx_p[ibin]-Tx[ibin])*(Tx_p[ibin]-Tx[ibin])));

      if(T0_prior<Tx[ibin]){
	T0_prior=Tx[ibin];
        rcool_prior=rmid[ibin];
      }
    
    }
    
     rt_prior=2.;
     Tmin_prior=Tx[0]/T0_prior;
     acool_prior=log(T0_prior/Tx[0])/log(rcool_prior/1e-2);
      

    fin.close();

    //----------------------------------------------

    ifstream fin2;
    fin2.open(nename);
   
    if(!fin2){fprintf(stderr,"NO FILE : %s\n", nename); return 1;}
    streampos current2;

    //count number
    colum1[0]='#';
     do{
       current2 = fin2.tellg();
       fin2.getline(colum1,sizeof(colum1));
     }while(colum1[0]=='#');
    // 行数を数えて結果ファイルに書く
    fin2.seekg(current2);
    istreambuf_iterator<char>  inF2(fin2);
    int Nline = count( inF2, EOF_MARKER, '\n' );
    fin2.seekg(current2);

    std::vector<double> rline_pix(Nline);
    std::vector<double> rline(Nline);
    std::vector<double> rline_Mpc(Nline);
    std::vector<double> n(Nline);
    std::vector<double> n_p(Nline);
    std::vector<double> n_m(Nline);
    std::vector<double> n2(Nline);
    std::vector<double> n2_p(Nline);
    std::vector<double> n2_m(Nline);
    std::vector<double> sx(Nline);
 

    int iline;
    for(iline=0;iline<Nline;iline++){
      fin2>>rline[iline]>>rline_Mpc[iline]>>n[iline]>>n_m[iline]>>n_p[iline];

      //I found that n2_m ! = n_m^2 and n2_p ! = n_p^2 
      n2[iline]=n[iline]*n[iline];
      n2_m[iline]=n_m[iline]*n_m[iline];
      n2_p[iline]=n_p[iline]*n_p[iline];
      if(n2[iline]<0){n2[iline]=1e-10;}
      if(n2_p[iline]<0){n2_p[iline]=1e-10;}
      if(n2_m[iline]<0){n2_m[iline]=1e-10;}
    }

    double amtoMpc=rline_Mpc[0]/rline[0];
    cout<<" amtoMpc "<<amtoMpc<<endl;
    tk::spline n2spline;
    n2spline.set_points(rline,n2);  
    tk::spline n2spline_p;
    n2spline_p.set_points(rline,n2_p);  
    tk::spline n2spline_m;
    n2spline_m.set_points(rline,n2_m);  
    
    fin2.close();

    //SETUP maximum radius;
    RMAX=rline[Nline-1];

  //-----------------------------------------------------------------------
  //fitting
  //-----------------------------------------------------------------------
   int Npara=4;
   //1st fitting w/ 3D model
   //       

     double T0=T0_prior;
     double Tmin=Tmin_fixed;
     double rcool=rcool_fixed;
     double rt=rt_prior;
     double acool=acool_fixed;
     double aT=aT_prior;
     double bT=bT_fixed;
     double cT=cT_prior;
     double factor_rcool=-0.6;

     Chisq_3d fFCN0(rmid,rinn,rout,Tx,dTx,n2spline);
     
     MnUserParameters upar0;    
     upar0.Add("T0", T0, 0.1);
     upar0.Add("rt", rt, 0.2);
     upar0.Add("aT", aT, 0.05);
     upar0.Add("cT", cT, 0.2); 

     upar0.SetLimits("cT", 0.0,3.);

     MnMigrad migrad0(fFCN0, upar0);
  
     //migrad0.Fix("Tmin");
     //          migrad0.Fix("factor_rcool");
     //  migrad0.Fix("aT");
     //migrad0.Fix("ct");
     
     
     // Minimize
     FunctionMinimum min0 = migrad0();

     // output
     std::cout<<"minimum: "<<min0<<std::endl;
     T0= min0.UserState().Value("T0");
     rt= min0.UserState().Value("rt");
     aT= min0.UserState().Value("aT");
     cT= min0.UserState().Value("cT");
     
     //------------------------

     Chisq fFCN(rmid,rinn,rout,Tx,dTx,n2spline);
    
     MnUserParameters upar;
     upar.Add("T0", T0, 1.4);
     upar.Add("rt", rt, 0.2);
     upar.Add("aT", aT, 0.05);
     upar.Add("cT", cT, 0.2); 


     upar.SetLowerLimit("T0", 0.);
     upar.SetLowerLimit("rt", 0.);
     //upar.SetLimits("aT", -1.,1.);
     upar.SetLimits("cT", 0.,3);

     // create MIGRAD minimizer
     MnMigrad migrad(fFCN, upar);

     // Minimize
     FunctionMinimum min = migrad();

     // output
     std::cout<<"minimum: "<<min<<std::endl;
  

     T0= min.UserState().Value("T0");
     rt= min.UserState().Value("rt");
     aT= min.UserState().Value("aT");
     cT= min.UserState().Value("cT");
   
     Tmin=Tmin_fixed;
     rcool=rcool_fixed;
     acool=acool_fixed;
     bT=bT_fixed; 
     
     int i,j;
     ROOT::Minuit2::MnUserCovariance covar=min.UserCovariance();
     double Cov[Npara][Npara];
     double rCov[Npara][Npara];
     for(i=0;i<Npara;i++){
     for(j=0;j<Npara;j++){
     Cov[i][j]=covar(i,j); 
     rCov[i][j]=covar(i,j)/sqrt(covar(i,i)*covar(j,j));
     }
     }

     double chisq=min.Fval();
     int Ndof=Nbin-Npara;
     double redchisq=chisq/Ndof; 
     double prob=gsl_sf_gamma_inc_Q(Ndof/2.,chisq/2.);

     char *now;
     tm *newtime;
     time_t aclock;
     time(&aclock);
     newtime = localtime(&aclock);  
     now     = asctime(newtime);    
     
     
    //----------------------------------------------------------
    //fitting w/ n2_p

     Chisq fFCN_p(rmid,rinn,rout,Tx,dTx,n2spline_p);
    
     MnUserParameters upar_p;
     upar_p.Add("T0", T0, 1.4);
     upar_p.Add("rt", rt, 0.2);
     upar_p.Add("aT", aT, 0.2); 
     upar_p.Add("cT", cT, 0.2); 
    

     upar_p.SetLowerLimit("T0", 0.);
     upar_p.SetLowerLimit("rt", 0.);
     //upar_p.SetLimits("aT", -1.,1.);
     upar_p.SetLimits("cT", 0.,3);

     // create MIGRAD minimizer
     MnMigrad migrad_p(fFCN_p, upar_p);

     // Minimize
     FunctionMinimum min_p = migrad_p();

     // output
     std::cout<<"minimum: "<<min_p<<std::endl;
  

     double T0_p= min_p.UserState().Value("T0");
     double rt_p= min_p.UserState().Value("rt");
     double aT_p= min_p.UserState().Value("aT");
     double cT_p= min_p.UserState().Value("cT");
   
     double factor_rcool_p= factor_rcool_fixed;
     double Tmin_p= Tmin_fixed;   
     double rcool_p=rcool_fixed;
     double acool_p=acool_fixed;
     double bT_p=bT_fixed;

    //----------------------------------------------------------
    //fitting w/ n2_m

     Chisq fFCN_m(rmid,rinn,rout,Tx,dTx,n2spline_m);
    
     MnUserParameters upar_m;
     upar_m.Add("T0", T0, 1.4);
     upar_m.Add("rt", rt, 0.2);
     upar_m.Add("aT", aT, 0.2); 
     upar_m.Add("cT", cT, 0.2); 
    

     upar_m.SetLowerLimit("T0", 0.);
     upar_m.SetLowerLimit("rt", 0.);
     //upar_m.SetLimits("aT", -1.,1.);
     upar_m.SetLimits("cT", 0.,3);

     // create MIGRAD minimizer
     MnMigrad migrad_m(fFCN_m, upar_m);

     // Minimize
     FunctionMinimum min_m = migrad_m();

     // output
     std::cout<<"minimum: "<<min_m<<std::endl;
  

     double T0_m= min_m.UserState().Value("T0");
     double rt_m= min_m.UserState().Value("rt");
     double aT_m= min_m.UserState().Value("aT");
     double cT_m= min_m.UserState().Value("cT");

   
     double factor_rcool_m= factor_rcool_fixed;
     double Tmin_m= Tmin_fixed;   
     double rcool_m=rcool_fixed;
     double acool_m=acool_fixed;
     double bT_m=bT_fixed;



     
    //===========================================================================================================
    char log_name[1000];
    sprintf(log_name, "Tresult/Tx_bestfit.log");
    FILE *fp_log = fopen( log_name, "w" );

    fprintf(fp_log, " #------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, " # LOG %s", now);
    fprintf(fp_log, " #------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, " input temp file                                : %s\n",tempname);
    fprintf(fp_log, " input n2 file                                  : %s\n",nename);
    fprintf(fp_log, " #------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "   model T3d= T0*(r/rt)^-aT*(1+(r/rt)^bT)^(-cT/bT); \n");

    fprintf(fp_log, "   Npara        (number of parameters)          : %d \n",Npara);
    fprintf(fp_log, "   T0           [keV]                           : %.8e \n",T0);
    fprintf(fp_log, "   rt           [arcmin]                        : %.8e \n",rt);
    fprintf(fp_log, "   aT                                           : %.8e\n",aT);
    fprintf(fp_log, "   cT                                           : %.8e \n",cT);
    fprintf(fp_log, "   \n");
    // fprintf(fp_log, "   Tmin (fixed)                                 : %.8e \n",Tmin);
    //fprintf(fp_log, "   factor_rcool (fixed)                         : %.8e \n",factor_rcool);
    // fprintf(fp_log, "   acool (fixed)                                : %.8e \n",acool);
    fprintf(fp_log, "   bT (fixed)                                   : %.8e\n",bT);
    //fprintf(fp_log, "   rcool (fixed)                                : %.8e\n",rcool);
    fprintf(fp_log, " #------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "   Covariance\n");
    for(i=0;i<Npara;i++){
    for(j=0;j<Npara;j++){
      fprintf(fp_log, "%.8e ",Cov[i][j]);
    }
      fprintf(fp_log, "\n");
    }
    fprintf(fp_log, " #------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "   Coefficient\n");
    for(i=0;i<Npara;i++){
    for(j=0;j<Npara;j++){
      fprintf(fp_log, "%.8e ",rCov[i][j]);
    }
      fprintf(fp_log, "\n");
    }
    fprintf(fp_log, " #------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "   chisq                                        : %.5lf \n",chisq);
    fprintf(fp_log, "   redchisq                                     : %.5lf \n",redchisq);
    fprintf(fp_log, "   Probability Q                                : %.8e \n",prob);

    fclose(fp_log);

    //----------------------------------------------------------

    char output_name0[1000];
    sprintf(output_name0, "Tresult/Tx_profile.dat");
    FILE *fp_output0 = fopen( output_name0, "w" );

    char output_name[1000];
    sprintf(output_name, "Tresult/Tx_bestfit.dat");
    FILE *fp_output = fopen( output_name, "w" );

    char output_name2[1000];
    sprintf(output_name2, "Tresult/T3d_bestfit.dat");
    FILE *fp_output2 = fopen( output_name2, "w" );

    Tew Tew_best(T0,Tmin,factor_rcool,rt,acool,aT,bT,cT,n2spline);
    Tew Tew_p(T0_p,Tmin_p,factor_rcool_p,rt_p,acool_p,aT_p,bT_p,cT_p,n2spline_p);    
    Tew Tew_m(T0_m,Tmin_m,factor_rcool_m,rt_m,acool_m,aT_m,bT_m,cT_m,n2spline_m);    

    T3d T3d_best(T0,Tmin,factor_rcool,rt,acool,aT,bT,cT);     
    T3d T3d_p(T0_p,Tmin_p,factor_rcool_p,rt_p,acool_p,aT_p,bT_p,cT_p);    
    T3d T3d_m(T0_m,Tmin_m,factor_rcool_m,rt_m,acool_m,aT_m,bT_m,cT_m);    
 
    std::vector<double> T3dmodel(Nline);  
    std::vector<double> dT3dmodel(Nline);  
    std::vector<double> dT3dmodel1(Nline);  
    std::vector<double> dT3dmodel2(Nline);  

    int Nbin2=60;
    double drmodel=(rout[rout.size()-1]-rinn[0])/(Nbin2-1);
    std::vector<double> rmodel(Nbin2);  
    std::vector<double> Tewmodel(Nbin2);  
    std::vector<double> dTewmodel(Nbin2);  
    std::vector<double> dTewmodel1(Nbin2); 
    //std::vector<double> dTewmodel1_p(Nbin2);  
    //std::vector<double> dTewmodel1_m(Nbin2);  
    std::vector<double> dTewmodel2(Nbin2);  



    //error calculation
    gsl_matrix *cov=gsl_matrix_alloc (Npara, Npara);
     for(i=0;i<Npara;i++){
      for(j=0;j<Npara;j++){
      gsl_matrix_set (cov, i, j, Cov[i][j]);
    }}
    gsl_linalg_cholesky_decomp(cov);


    for(ibin=0;ibin<Nbin2;ibin++){
      rmodel[ibin]=rinn[0]+drmodel*ibin;
      Tewmodel[ibin]=Tew_best(rmodel[ibin]);
      dTewmodel1[ibin]=0.5*(pow(Tew_p(rmodel[ibin])-Tewmodel[ibin],2)+pow(Tew_m(rmodel[ibin])-Tewmodel[ibin],2));
      dTewmodel2[ibin]=0.;
      dTewmodel[ibin]=0.;
    }
    

   for(iline=0;iline<Nline;iline++){  
      T3dmodel[iline]=T3d_best(rline[iline]);
      dT3dmodel1[iline]=0.5*(pow(T3d_p(rline[iline])-T3dmodel[iline],2)+pow(T3d_m(rline[iline])-T3dmodel[iline],2));      
      dT3dmodel2[iline]=0.;
      dT3dmodel[iline]=0.;
   }

     int iboot;
     int Nboot=10000;
     double T0_tmp;
     double rt_tmp;
     double aT_tmp;
     double cT_tmp;

     double rcool_tmp=rcool_fixed;
     double acool_tmp=acool_fixed;
     double Tmin_tmp=Tmin_fixed;
     double factor_rcool_tmp=factor_rcool_fixed;
     double bT_tmp=bT_fixed; 

     gsl_rng *gBaseRand; 
     unsigned long randSeed;
     gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);  
     srand(time(NULL));                    // initialization for rand() 
     randSeed = rand();                    // returns a non-negative integer
     gsl_rng_set (gBaseRand, randSeed);    // seed the PRNG
     
     double gauss0,gauss1,gauss2,gauss3,gauss4,gauss5;

   
     for(iboot=0;iboot<Nboot;iboot++){
       //    if(iboot%1000==0){fprintf(stdout,"calculating error %d / %d\n",iboot,Nboot);}
    do{

    gauss0=gsl_ran_gaussian(gBaseRand,1.);
    gauss1=gsl_ran_gaussian(gBaseRand,1.);
    gauss1=gsl_ran_gaussian(gBaseRand,1.);
    gauss2=gsl_ran_gaussian(gBaseRand,1.);
    gauss3=gsl_ran_gaussian(gBaseRand,1.);
    gauss4=gsl_ran_gaussian(gBaseRand,1.);

   
    T0_tmp=gsl_matrix_get (cov, 0, 0)*gauss0+T0;
    rt_tmp  =gsl_matrix_get (cov, 1, 0)*gauss0+gsl_matrix_get (cov, 1, 1)*gauss1+rt;
    aT_tmp=gsl_matrix_get (cov, 2, 0)*gauss0+gsl_matrix_get (cov, 2, 1)*gauss1+gsl_matrix_get (cov, 2, 2)*gauss2+aT;
    cT_tmp=gsl_matrix_get (cov, 3, 0)*gauss0+gsl_matrix_get (cov, 3, 1)*gauss1+gsl_matrix_get (cov, 3, 2)*gauss2+gsl_matrix_get (cov, 3, 3)*gauss3+cT;
    
    //     }while(T0_tmp<0. || Tmin_tmp<0. || rt_tmp<0 || rcool_tmp<0 || rcool_tmp>rt_tmp);
     }while(T0_tmp<0. || rt_tmp<0 );

    Tew Tew_tmp(T0_tmp,Tmin_fixed,factor_rcool_fixed,rt_tmp,acool_fixed,aT_tmp,bT_fixed,cT_tmp,n2spline);
    T3d T3d_tmp(T0_tmp,Tmin_fixed,factor_rcool_fixed,rt_tmp,acool_fixed,aT_tmp,bT_fixed,cT_tmp);     

    
     for(ibin=0;ibin<Nbin2;ibin++){   
       dTewmodel2[ibin]+=pow(Tew_tmp(rmodel[ibin])-Tewmodel[ibin],2);
      }
    

     for(iline=0;iline<Nline;iline++){   
       dT3dmodel2[iline]+=pow(T3d_tmp(rline[iline])-T3dmodel[iline],2);
      }
	    
  }
     

     fprintf(fp_output0,"# r[arcmin] r_m r_p r[Mpch_%.0f^{-1}] r_Mpc_m r_Mpc_p Tx[keV] Tx_m Tx_p\n",H0);
    fprintf(fp_output,"# r[arcmin] r[Mpch_%.0f^{-1}] Tew[keV] Tew_m Tew_p\n",H0);
    fprintf(fp_output2,"# r[arcmin] r[Mpch_%.0f^{-1}] T3d[keV] T3d_m T3d_p\n",H0);
    
    for(ibin=0;ibin<Nbin;ibin++){  
      fprintf(fp_output0,"%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n", rmid[ibin],rinn[ibin],rout[ibin],rmid[ibin]*amtoMpc,rinn[ibin]*amtoMpc,rout[ibin]*amtoMpc,Tx[ibin],Tx_m[ibin],Tx_p[ibin]);
    }

    for(ibin=1;ibin<Nbin2;ibin++){     //ibin=0 --> rmodel=0
        dTewmodel[ibin]=sqrt(dTewmodel2[ibin]/(Nboot-1.)+dTewmodel1[ibin]);        
      //       dTewmodel[ibin]=sqrt(dTewmodel1[ibin]);        
      fprintf(fp_output,"%.8e %.8e %.8e %.8e %.8e\n",rmodel[ibin],rmodel[ibin]*amtoMpc,Tewmodel[ibin],Tewmodel[ibin]-dTewmodel[ibin],Tewmodel[ibin]+dTewmodel[ibin]);
      }
    
    for(iline=0;iline<Nline;iline++){  
      dT3dmodel[iline]=sqrt(dT3dmodel2[iline]/(Nboot-1.)+dT3dmodel1[iline]); 
      fprintf(fp_output2,"%.8e %.8e %.8e %.8e %.8e\n",rline[iline],rline[iline]*amtoMpc,T3dmodel[iline],T3dmodel[iline]-dT3dmodel[iline],T3dmodel[iline]+dT3dmodel[iline]);         
   }
     fclose(fp_output0);
     fclose(fp_output);
     fclose(fp_output2);
     

}




