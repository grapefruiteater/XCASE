//icpc  fit_sx_2beta_mos1+mos2+pn.cxx  -I/usr/local/include/Minuit2 -I/usr/local/include/Math -L/usr/local/lib/ -lMinuit2  $CPGPLOT_FLAGS $GSL_FLAGS $MINUIT_FLAGS $CFITSIO_FLAGS -parallel -o "fit_sx_2beta_mos1+mos2+pn"

#include<iostream>
#include<fstream>
#include<algorithm>
#include<cstdio>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include <vector>
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

double rate_th=0.; //do not count mask region  !!! SETUP PARAMETER !!!

//-----------------------------------------
//Vikhlinin
//-----------------------------------------
double gamma1=3.; // !!! SETUP PARAMETER !!!
double ALPHA=0;
double EPSILON=0;
double NORM2=0.;
double BETA2=0.;
double R02=0.;

//-----------------------------------------
// Afakeit parameter
//-----------------------------------------
double Amos1_norm=100.;
double Amos2_norm=100.;
double Apn_norm=300.;
double c1_A=-1e-2;
double c2_A=0;


//=======================================================================

//=======================================================================
// data map
//=======================================================================
std::vector<double> mos1_psfmap;
std::vector<double> mos1_imgmap;
std::vector<double> mos1_sigmap;

std::vector<double> mos2_psfmap;
std::vector<double> mos2_imgmap;
std::vector<double> mos2_sigmap;

std::vector<double> pn_psfmap;
std::vector<double> pn_imgmap;
std::vector<double> pn_sigmap;

//radial profile
std::vector<double> r1;
std::vector<double> r2;
std::vector<double> rmid;
std::vector<double> rmid_am;
std::vector<double> sx_mos1;
std::vector<double> dsx_mos1;
std::vector<double> area_mos1;
std::vector<double> sx_mos2;
std::vector<double> dsx_mos2;
std::vector<double> area_mos2;
std::vector<double> sx_pn;
std::vector<double> dsx_pn;
std::vector<double> area_pn;
std::vector<double> sx_ave;
std::vector<double> dsx_sys;
std::vector<double> sx_sig;
std::vector<double> sx_sig_stat;

std::vector<double> sx_mos10;
std::vector<double> dsx_mos10;
std::vector<double> sx_mos20;
std::vector<double> dsx_mos20;
std::vector<double> sx_pn0;
std::vector<double> dsx_pn0;
std::vector<double> sx_ave0;
std::vector<double> dsx_sys0;

std::vector<double> r1model;
std::vector<double> r2model;
std::vector<double> area_model;
std::vector<double> Sx_model;
std::vector<double> Sx_model_mos1;
std::vector<double> Sx_model_mos2;
std::vector<double> Sx_model_pn;
std::vector<double> Sx_model_conv_mos1;
std::vector<double> Sx_model_conv_mos2;
std::vector<double> Sx_model_conv_pn;

//=======================================================================
int i,j,ix,iy,ipsf,ixpsf,iypsf;
long fpixel[2]={1,1};
long inc[2]={1,1};
int anynull;
int status=0;

//mos1,mos2,pn
int Nximg_mos1,Nyimg_mos1;
int Nximg_mos2,Nyimg_mos2;
int Nximg_pn,Nyimg_pn;
int Nxsig_mos1,Nysig_mos1;
int Nxsig_mos2,Nysig_mos2;
int Nxsig_pn,Nysig_pn;
int Nxpsf_mos1,Nypsf_mos1;
int Nxpsf_mos2,Nypsf_mos2;
int Nxpsf_pn,Nypsf_pn;


int Nx,Ny,N;
int Nxpsf,Nypsf,Npsf;
int ixpsf_c,iypsf_c; //ixpsf_c=Nxpsf/2-1;

//pix -> arcsec
double pixtoas_mos1_img;
double pixtoas_mos2_img;
double pixtoas_pn_img;
double pixtoas_mos1_sig;
double pixtoas_mos2_sig;
double pixtoas_pn_sig;
double pixtoas_mos1_psf;
double pixtoas_mos2_psf;
double pixtoas_pn_psf;

// pixel scale for our calculation
double pixtoas;
double pixtoam;
double pixtoas_psf;
double pixtoam_psf;

std::vector<double> x; //x coordinate
std::vector<double> y; //x coordinate
std::vector<double> d;  //distance from the center

//modelmap
int Nxmodel;
int Nymodel;
int Nmodel=Nxmodel*Nymodel;
double pixtoas_model;
double pixtoam_model;
double xcmodel;
double ycmodel;
double xcmodel_psf;
double ycmodel_psf;

std::vector<double> modelmap_mos1;
std::vector<double> modelmap_mos2;
std::vector<double> modelmap_pn;
std::vector<double> xmodel;
std::vector<double> ymodel;
std::vector<double> dmodel;
std::vector<double> dmodel_rateimg;

std::vector<double> dconv;

//convolution max radius
double dconv_max_am=2.; //arcmin !!! SETUP PARAMETER !!!
double dconv_max;

//psf map
double wgt_mos1;
double wgt_mos2;
double wgt_pn;

//spline for convolution
int Nspline=50;                         // !!KEY PARAMETER for Calculating SPEED!!
int ispline;
double rspline_min=1e-5; //arcmin       // !!KEY PARAMETER for Calculating SPEED!!
double rspline_max=40.;  //arcmin       // !!KEY PARAMETER for Calculating SPEED!!
double rspline_width=(log10(rspline_max)-log10(rspline_min))/Nspline;
std::vector<double> rspline(Nspline);
std::vector<double> Sxspline_mos1(Nspline);
std::vector<double> Sxspline_mos2(Nspline);
std::vector<double> Sxspline_pn(Nspline);
tk::spline s_mos1;
tk::spline s_mos2;
tk::spline s_pn;

void resize_sx(int Nbin);
void setup()
{
    for(ix=0;ix<Nx;ix++){
      x[ix]=1.*ix+1.;
     }
    for(iy=0;iy<Ny;iy++){
      y[iy]=1.*iy+1.;
    } 
    for(ix=0;ix<Nxmodel;ix++){
      xmodel[ix]=1.*ix+1.;
    }
    for(iy=0;iy<Nymodel;iy++){
      ymodel[iy]=1.*iy+1.;
    }
  
   for(ispline=0;ispline<Nspline;ispline++){
    rspline[ispline]=rspline_min*pow(10,rspline_width*ispline)/pixtoam;    //pix  
   }
}


//=======================================================================
// Read emissivity data 
//=======================================================================

void readAfakeitdata(char dataname[],std::vector<double>& r_A,std::vector<double>& A)
{
    ifstream fin;
    fin.open(dataname);
   
    if(!fin){fprintf(stderr,"NO FILE : %s\n", dataname);}
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
    int N_A = count( inF, EOF_MARKER, '\n' );
    fin.seekg(current);

    double r_A_tmp;
    double A_tmp;
    double tmp;
    for(i=0;i<N_A;i++){
      fin>>r_A_tmp>>tmp>>tmp>>A_tmp;
      r_A.push_back(r_A_tmp);
      A.push_back(A_tmp);
    }
     fin.close();

}


//=======================================================================
// Read count_rateimage
//=======================================================================

void readimgmap(std::vector<double>& imgmap,char ratename[],int& nx,int& ny, double& pixtoas)
{
   fitsfile *rate_fits;
  
  if(fits_open_image(&rate_fits, ratename, READONLY, &status)){
    fits_report_error(stderr, status);
    exit(status);
  }
  int bitpix_rate;  // INPUT BITPIX
  int naxis_rate;   // Number of DIMENSION
  long naxes_rate[2];  //NUMBER OF PIXs (NAXIS1, NAXIS2)

  if(fits_get_img_param(rate_fits, 2, &bitpix_rate, &naxis_rate, naxes_rate, &status))
   {
      fprintf(stderr,
       "%s (line %d): Problem reading information about current HDU from file \n",
	      __FILE__, __LINE__);     
      fits_report_error(stderr, status);
      exit(status);
      }  
  
  if ( naxis_rate != 2) {
      fprintf(stderr, "Are you insane?  We expect an image here.\n");
       exit(status);
  }

  nx=naxes_rate[0];
  ny=naxes_rate[1];

  double CDELT1;
  double CDELT2;
  fits_read_key(rate_fits, TDOUBLE, "CDELT1", &CDELT1, NULL, &status ); 
  fits_read_key(rate_fits, TDOUBLE, "CDELT2", &CDELT2, NULL, &status ); 
  pixtoas=sqrt(fabs(CDELT1*CDELT2))*3600.;
  status=0;

  double *img=new double[nx*ny];
  fits_read_subset(rate_fits, TDOUBLE, fpixel, naxes_rate, inc, NULL, img, &anynull, &status);
  fits_close_file(rate_fits, &status);

  for(i=0;i<nx*ny;i++){
    imgmap.push_back(img[i]);
  }
  
  
}


//=======================================================================
// convolution
//=======================================================================
// input   : input map
// output  : convolved map
// sxmodel : convolved sx profile
// input   : input map
// output  : convolved map
// sxmodel : convolved sx profile
void conv(const std::vector<double>& input_mos1,
          const std::vector<double>& input_mos2,
          const std::vector<double>& input_pn,
          const std::vector<double>& r1model,
          const std::vector<double>& r2model,
          const std::vector<double>& area_model, 
          std::vector<double>& sxmodel_mos1,
          std::vector<double>& sxmodel_mos2,
          std::vector<double>& sxmodel_pn)
{
  std::vector<double> map_conv_mos1(Nmodel);
  std::vector<double> norm_conv_mos1(Nmodel);
  std::vector<double> map_conv_mos2(Nmodel);
  std::vector<double> norm_conv_mos2(Nmodel);
  std::vector<double> map_conv_pn(Nmodel);
  std::vector<double> norm_conv_pn(Nmodel);
  int ixnew,iynew,inew;
  int ibin;
  int Nbin=r1model.size();
 
  //initialize
  for(i=0;i<Nmodel;i++){
    map_conv_mos1[i]=0.;
    norm_conv_mos1[i]=0.;
    map_conv_mos2[i]=0.;
    norm_conv_mos2[i]=0.;
    map_conv_pn[i]=0.;
    norm_conv_pn[i]=0.;
  }
  for(ibin=0;ibin<Nbin;ibin++){	
    sxmodel_mos1[ibin]=0.;
    sxmodel_mos2[ibin]=0.;
    sxmodel_pn[ibin]=0.;
  }

  //convolution
  for(ix=0;ix<Nxmodel;ix++){
  for(iy=0;iy<Nymodel;iy++){ 
        i=ix+Nxmodel*iy;
	if(dconv[i]<dconv_max){
          for(ixpsf=0;ixpsf<Nxpsf;ixpsf++){
             for(iypsf=0;iypsf<Nypsf;iypsf++){ 
                  ipsf=ixpsf+Nxpsf*iypsf;
                  ixnew=ix+(ixpsf-ixpsf_c);
                  iynew=iy+(iypsf-iypsf_c);
                  inew=ixnew+Nxmodel*iynew;

                  if ( ixnew >= 0 && ixnew< Nxmodel && iynew >=0 && iynew<Nymodel ){
		    map_conv_mos1[inew]+=input_mos1[i]*mos1_psfmap[ipsf]; 
		    norm_conv_mos1[inew]+=mos1_psfmap[ipsf];
		    map_conv_mos2[inew]+=input_mos2[i]*mos2_psfmap[ipsf]; 
		    norm_conv_mos2[inew]+=mos2_psfmap[ipsf];
		    map_conv_pn[inew]+=input_pn[i]*pn_psfmap[ipsf]; 
		    norm_conv_pn[inew]+=pn_psfmap[ipsf];

		  }

	   }}
	     }
        else{
	  map_conv_mos1[i]+=input_mos1[i]*wgt_mos1;
          norm_conv_mos1[i]+=wgt_mos1;
	  map_conv_mos2[i]+=input_mos2[i]*wgt_mos2;
          norm_conv_mos2[i]+=wgt_mos2;
	  map_conv_pn[i]+=input_pn[i]*wgt_pn;
          norm_conv_pn[i]+=wgt_pn;

	  }
  }}

  for(ix=0;ix<Nxmodel;ix++){
  for(iy=0;iy<Nymodel;iy++){ 
    i=ix+Nxmodel*iy;
    map_conv_mos1[i]/=norm_conv_mos1[i];
    map_conv_mos2[i]/=norm_conv_mos2[i];
    map_conv_pn[i]/=norm_conv_pn[i];    

     for(ibin=0;ibin<Nbin;ibin++){
       if(dmodel[i]>r1model[ibin] &&  dmodel[i]<=r2model[ibin]){
 	sxmodel_mos1[ibin]+=map_conv_mos1[i];
 	sxmodel_mos2[ibin]+=map_conv_mos2[i];
 	sxmodel_pn[ibin]+=map_conv_pn[i];
        break;
      }
    }}
  }

  for(ibin=0;ibin<Nbin;ibin++){	
        sxmodel_mos1[ibin]/=area_model[ibin];
        sxmodel_mos2[ibin]/=area_model[ibin];
        sxmodel_pn[ibin]/=area_model[ibin];
  }

}



//=======================================================================
namespace ROOT {

   namespace Minuit2 {

     class n2_Vikhlinin {

public:

       n2_Vikhlinin(double norm1, double r01, double alpha, double beta1, double norm2, double r02, double beta2,double epsilon) :
	 norm1(norm1), r01(r01), alpha(alpha), beta1(beta1), norm2(norm2), r02(r02), beta2(beta2),epsilon(epsilon){}
  ~n2_Vikhlinin() {}


  double operator()(double r) const {
   
    return norm1*pow(r/r01,-alpha)/pow(1.+r*r/r01/r01,3.*beta1-alpha/2.)/pow(1.+pow(r/r01,gamma1),epsilon/gamma1)+norm2/pow(1.+r*r/r02/r02,3.*beta2);

  }

private:

  double norm1;
  double r01;
  double alpha;
  double beta1;
  double norm2;
  double r02;
  double beta2;
  double epsilon;

};

     class n_Vikhlinin {

public:

       n_Vikhlinin(double norm1, double r01, double alpha, double beta1, double norm2, double r02, double beta2,double epsilon) :
	 norm1(norm1), r01(r01), alpha(alpha), beta1(beta1), norm2(norm2), r02(r02), beta2(beta2),epsilon(epsilon) {}
  ~n_Vikhlinin() {}


  double operator()(double r) const {
   
    return sqrt(norm1*pow(r/r01,-alpha)/pow(1.+r*r/r01/r01,3.*beta1-alpha/2.)/pow(1.+pow(r/r01,gamma1),epsilon/gamma1)+norm2/pow(1.+r*r/r02/r02,3.*beta2));

  }

private:

  double norm1;
  double r01;
  double alpha;
  double beta1;
  double norm2;
  double r02;
  double beta2;
  double epsilon;

};


class Sx_Vikhlinin {

public:

  Sx_Vikhlinin(double norm1, double r01, double alpha, double beta1, double norm2, double r02, double beta2,double epsilon) :
    norm1(norm1), r01(r01), alpha(alpha), beta1(beta1), norm2(norm2), r02(r02), beta2(beta2),epsilon(epsilon) {}
  ~Sx_Vikhlinin() {}


  double operator()(double r) const {

    //r/R=x = cosh(t)
    double tmin=0.; 
    double tmax=acosh(20000./r); //integration along line-of-sight with unit of pix=pixtoam
    int Nline=60; 
    double dt=(tmax-tmin)/Nline;
    int iline;
    double sum=0.;
    double ttmp;
    double rtmp;
 
    n2_Vikhlinin n2_Vikhlinin(norm1, r01,alpha,beta1,norm2,r02,beta2,epsilon);
   
    sum=n2_Vikhlinin(cosh(tmin)*r)*cosh(tmin)+n2_Vikhlinin(cosh(tmax)*r)*cosh(tmax);


    for(iline=1;iline<Nline;iline+=2){
      ttmp=tmin+iline*dt;
      rtmp=cosh(ttmp)*r;
      sum+=4*n2_Vikhlinin(rtmp)*cosh(ttmp);
    }

    for(iline=2;iline<Nline-1;iline+=2){
      ttmp=tmin+iline*dt;
      rtmp=cosh(ttmp)*r;
      sum+=2*n2_Vikhlinin(rtmp)*cosh(ttmp);
    }
    
      sum=2.*dt*sum*r/3.;

   
    return sum;

  }

private:

  double norm1;
  double r01;
  double alpha;
  double beta1;
  double norm2;
  double r02;
  double beta2;
  double epsilon;
};

//----------------------------------------------------------------------------------------------------

class n2_2beta {

public:

  n2_2beta(double norm1, double r01, double beta1, double norm2, double r02, double beta2) :
    norm1(norm1), r01(r01), beta1(beta1), norm2(norm2), r02(r02), beta2(beta2){}
  ~n2_2beta() {}


  double operator()(double r) const {
   
    return norm1*pow(1+r*r/r01/r01,-3.*beta1)+norm2*pow(1+r*r/r02/r02,-3.*beta2);

  }

private:

  double norm1;
  double r01;
  double beta1;
  double norm2;
  double r02;
  double beta2;


};

class n_2beta {

public:

  n_2beta(double norm1, double r01, double beta1, double norm2, double r02, double beta2) :
    norm1(norm1), r01(r01), beta1(beta1), norm2(norm2), r02(r02), beta2(beta2){}
  ~n_2beta() {}


  double operator()(double r) const {
   
    return sqrt(norm1*pow(1+r*r/r01/r01,-3.*beta1)+norm2*pow(1+r*r/r02/r02,-3.*beta2));

  }

private:

  double norm1;
  double r01;
  double beta1;
  double norm2;
  double r02;
  double beta2;


};


class Sx_2beta {

public:

  Sx_2beta(double norm1, double r01, double beta1, double norm2, double r02, double beta2) :
          norm1(norm1), r01(r01), beta1(beta1), norm2(norm2), r02(r02), beta2(beta2){}
  ~Sx_2beta() {}


  double operator()(double r) const {
   
    return sqrt(Pi)*(gsl_sf_gamma(3*beta1-0.5)*r01*norm1*pow(1+r*r/r01/r01,0.5-3*beta1)/gsl_sf_gamma(3*beta1)
      +gsl_sf_gamma(3*beta2-0.5)*r02*norm2*pow(1+r*r/r02/r02,0.5-3*beta2)/gsl_sf_gamma(3*beta2));
    
  }

private:

  double norm1;
  double r01;
  double beta1;
  double norm2;
  double r02;
  double beta2;

};


class n2_beta {

public:

  n2_beta(double norm1, double r01, double beta1) :
    norm1(norm1), r01(r01), beta1(beta1){}
  ~n2_beta() {}


  double operator()(double r) const {
   
    return norm1*pow(1+r*r/r01/r01,-3.*beta1);

  }

private:

  double norm1;
  double r01;
  double beta1;


};


class n_beta {

public:

  n_beta(double norm1, double r01, double beta1) :
    norm1(norm1), r01(r01), beta1(beta1){}
  ~n_beta() {}


  double operator()(double r) const {
   
    return sqrt(norm1*pow(1+r*r/r01/r01,-3.*beta1));

  }

private:

  double norm1;
  double r01;
  double beta1;


};


class Sx_beta {

public:

  Sx_beta(double norm1, double r01, double beta1) :
          norm1(norm1), r01(r01), beta1(beta1){}
  ~Sx_beta() {}


  double operator()(double r) const {
   
    return sqrt(Pi)*(gsl_sf_gamma(3*beta1-0.5)*r01*norm1*pow(1+r*r/r01/r01,0.5-3*beta1)/gsl_sf_gamma(3*beta1));
    
  }

private:

  double norm1;
  double r01;
  double beta1;
  double norm2;
  double r02;
  double beta2;

};



class Chisq : public FCNBase {

public:
  Chisq(const std::vector<double>& rmid,
        const std::vector<double>& sx_mos1,
	const std::vector<double>& dsx_mos1,
        const std::vector<double>& sx_mos2,
	const std::vector<double>& dsx_mos2,
        const std::vector<double>& sx_pn,
	const std::vector<double>& dsx_pn) : 
    rmid(rmid), sx_mos1(sx_mos1), dsx_mos1(dsx_mos1),sx_mos2(sx_mos2), dsx_mos2(dsx_mos2),sx_pn(sx_pn), dsx_pn(dsx_pn),fErrorDef(1.) {}

  ~Chisq() {}
  virtual double Up() const {return fErrorDef;} //log-likelihood
  virtual double operator()(const std::vector<double>&) const;

  void SetErrorDef(double def) {fErrorDef = def;}

private:

  std::vector<double> rmid;
  std::vector<double> sx_mos1;
  std::vector<double> dsx_mos1;
  std::vector<double> sx_mos2;
  std::vector<double> dsx_mos2;
  std::vector<double> sx_pn;
  std::vector<double> dsx_pn;
  double fErrorDef;

};

//----------------------------------------------------------------------------------------------------


double Chisq::operator()(const std::vector<double>& para) const {

  double chisq=0.;
  Sx_2beta Sx_mos1(para[0], para[1],para[2],para[3],para[4],para[5]);

  double sb_mos1=para[6];
  double sb_mos2=para[7];
  double sb_pn=para[8];

  int ibin;
 
  for(ispline=0;ispline<Nspline;ispline++){
   Sxspline_mos1[ispline]=Sx_mos1(rspline[ispline]);  
  }

  s_mos1.set_points(rspline,Sxspline_mos1);  

  //!!! setup 2d map 
  //initialize
  for(ix=0;ix<Nxmodel;ix++){
  for(iy=0;iy<Nymodel;iy++){ 
    i=ix+Nxmodel*iy;
    modelmap_mos1[i]=s_mos1(dmodel_rateimg[i]);
    modelmap_mos2[i]=modelmap_mos1[i];
    modelmap_pn[i]=modelmap_mos1[i];
  }
  }

  //convolution
  conv(modelmap_mos1,modelmap_mos2,modelmap_pn,r1model,r2model,area_model,Sx_model_conv_mos1,Sx_model_conv_mos2,Sx_model_conv_pn);

  for(ibin=0;ibin<r1model.size();ibin++){
      Sx_model_conv_mos1[ibin]+=sb_mos1;
      Sx_model_conv_mos2[ibin]+=sb_mos2;
      Sx_model_conv_pn[ibin]+=sb_pn;
      chisq+=(Sx_model_conv_mos1[ibin]-sx_mos1[ibin])*(Sx_model_conv_mos1[ibin]-sx_mos1[ibin])/dsx_mos1[ibin]/dsx_mos1[ibin]
	    +(Sx_model_conv_mos2[ibin]-sx_mos2[ibin])*(Sx_model_conv_mos2[ibin]-sx_mos2[ibin])/dsx_mos2[ibin]/dsx_mos2[ibin]
            +(Sx_model_conv_pn[ibin]-sx_pn[ibin])*(Sx_model_conv_pn[ibin]-sx_pn[ibin])/dsx_pn[ibin]/dsx_pn[ibin];
  }

  return chisq;

}



class Chisq_nonconv : public FCNBase {

public:
  Chisq_nonconv(const std::vector<double>& rmid,
        const std::vector<double>& sx_mos1,
	const std::vector<double>& dsx_mos1,
        const std::vector<double>& sx_mos2,
	const std::vector<double>& dsx_mos2,
        const std::vector<double>& sx_pn,
	const std::vector<double>& dsx_pn) : 
    rmid(rmid), sx_mos1(sx_mos1), dsx_mos1(dsx_mos1),sx_mos2(sx_mos2), dsx_mos2(dsx_mos2),sx_pn(sx_pn), dsx_pn(dsx_pn),fErrorDef(1.) {}

  ~Chisq_nonconv() {}
  virtual double Up() const {return fErrorDef;} //log-likelihood
  virtual double operator()(const std::vector<double>&) const;

  void SetErrorDef(double def) {fErrorDef = def;}

private:

  std::vector<double> rmid;
  std::vector<double> sx_mos1;
  std::vector<double> dsx_mos1;
  std::vector<double> sx_mos2;
  std::vector<double> dsx_mos2;
  std::vector<double> sx_pn;
  std::vector<double> dsx_pn;
  double fErrorDef;

};


double Chisq_nonconv::operator()(const std::vector<double>& para) const {

  // original  : Sx = c Afakeit(r) * Fcl + c*bkg
  // fitting   : Sx/(c Afakeit(r)) = Fcl + bkg/Afakit(r)

  double chisq=0.;
  Sx_2beta Sx_mos1(para[0], para[1],para[2],para[3],para[4],para[5]);
  
  double Sx_tmp; 
  double sb_mos1=para[6];
  double sb_mos2=para[7];
  double sb_pn=para[8];
  int ibin;
  int Nbin=rmid.size();

  for(ibin=0;ibin<Nbin;ibin++){
      Sx_tmp=Sx_mos1(rmid[ibin]);
      Sx_model_mos1[ibin]=Sx_tmp+sb_mos1;
      Sx_model_mos2[ibin]=Sx_tmp+sb_mos2;
      Sx_model_pn[ibin]=Sx_tmp+sb_pn;
      chisq+=(Sx_model_mos1[ibin]-sx_mos1[ibin])*(Sx_model_mos1[ibin]-sx_mos1[ibin])/dsx_mos1[ibin]/dsx_mos1[ibin]
	    +(Sx_model_mos2[ibin]-sx_mos2[ibin])*(Sx_model_mos2[ibin]-sx_mos2[ibin])/dsx_mos2[ibin]/dsx_mos2[ibin]
            +(Sx_model_pn[ibin]-sx_pn[ibin])*(Sx_model_pn[ibin]-sx_pn[ibin])/dsx_pn[ibin]/dsx_pn[ibin];
  }
  
  return chisq;

}



  }  // namespace Minuit2

}  // namespace ROOT


//=======================================================================
// main func
//=======================================================================

int main(int argc, char *argv[])
{
if(argc!=24){
    fprintf(stderr, "\n");
    fprintf(stderr, "Use : ./fit_sx_vikhlinin_conv_v3 mos1_img mos1_sig mos1_psf Afakeit_mos1 con_mos1 mos2_img mos2_sig mos2_psf Afakeit_mos2 con_mos2 pn_img pn_sig pn_psf Afakeit_pn con_pn pixtoas_psfmap xc yc xc_psfmap yc_psfmap rmin rmax Nbin\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      mos1_img         - count_rate mos1 image\n");
    fprintf(stderr, "      mos1_sig         - sigma  mos1 image\n");
    fprintf(stderr, "      mos1_psf         - psfmap mos1 image\n");
    fprintf(stderr, "      Afakeit_mos1     - emissivity file for mos1 \n");
    fprintf(stderr, "      con_mos1         - con for mos1\n");
    fprintf(stderr, "      mos2_img         - count_rate mos2 image\n");
    fprintf(stderr, "      mos2_sig         - sigma  mos2 image\n");
    fprintf(stderr, "      mos2_psf         - psfmap mos2 image\n");
    fprintf(stderr, "      Afakeit_mos2     - emissivity file for mos2 \n");
    fprintf(stderr, "      con_mos2         - con for mos2\n");
    fprintf(stderr, "      pn_img           - count_rate pn image\n");
    fprintf(stderr, "      pn_sig           - sigma  pn image\n");
    fprintf(stderr, "      pn_psf           - psfmap pn image\n");
    fprintf(stderr, "      con_pn           - con for pn\n");
    fprintf(stderr, "      Afakeit_pn       - emissivity file for pn \n");
    fprintf(stderr, "      pixtoas_psfmap   - pix size for psfmap [arcsec]\n");
    fprintf(stderr, "      xc               - x center of Sx profile [pix]\n");
    fprintf(stderr, "      yc               - y center of Sx profile [pix]\n");
    fprintf(stderr, "      xc_psf           - x center to make psf map [pix]\n");
    fprintf(stderr, "      yc_psf           - y center to make psf map [pix]\n");
    fprintf(stderr, "      rmin             - minium radius of Sx profile [arcmin]\n");
    fprintf(stderr, "      rmax             - minium radius of Sx profile [arcmin]\n");
    fprintf(stderr, "      Nbin             - number of bins\n");
    fprintf(stderr, "      \n");
    fprintf(stderr, "      output files in Sxresult\n");
    fprintf(stderr, "      Sx profile          : sx_profile.dat\n");
    fprintf(stderr, "      best-fit Sx profile : sx_bestfit.dat\n");
    fprintf(stderr, "      best-fit ne profile : ne_bestfit.dat\n");
    fprintf(stderr, "      log file            : se_bestfit.log\n");
    fprintf(stderr, "      \n");
    fprintf(stderr, " e.g. : OMP_NUM_THREADS=20  ./fit_sx_vikhlinin_mos1+mos2+pn rate-400-2300-mos1-sp.fits sigma-400-2300-mos1-sp.fits mos1_psf_2bin.fits Lambda_mos1.dat 0.975813 rate-400-2300-mos2-sp.fits sigma-400-2300-mos2-sp.fits mos2_psf_2bin.fits  Lambda_mos2.dat 1.0 rate-400-2300-pn-sp.fits sigma-400-2300-pn-sp.fits psf_pn_rebin.fits Lambda_pn.dat 0.938212 2 443.050510 413.031470 443.050510 413.031470 0.1 10 28\n");
    fprintf(stderr, "\n");
    return 1;
  }

   struct stat st;
   if(stat("Sxresult", &st)!=0){
     mkdir("Sxresult",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   }

  char mos1_ratename[1000];  
  sprintf(mos1_ratename,"%s",argv[1]);  
  char mos1_signame[1000];  
  sprintf(mos1_signame,"%s",argv[2]);
  char mos1_psfname[1000];
  sprintf(mos1_psfname,"%s",argv[3]);
  char Amos1_name[1000];
  sprintf(Amos1_name,"%s",argv[4]);
  double con_mos1=atof(argv[5]);

  char mos2_ratename[1000];  
  sprintf(mos2_ratename,"%s",argv[6]);  
  char mos2_signame[1000];  
  sprintf(mos2_signame,"%s",argv[7]);
  char mos2_psfname[1000];
  sprintf(mos2_psfname,"%s",argv[8]);
  char Amos2_name[1000];
  sprintf(Amos2_name,"%s",argv[9]);
  double con_mos2=atof(argv[10]);

  char pn_ratename[1000];  
  sprintf(pn_ratename,"%s",argv[11]);  
  char pn_signame[1000];  
  sprintf(pn_signame,"%s",argv[12]);
  char pn_psfname[1000];
  sprintf(pn_psfname,"%s",argv[13]);
  char Apn_name[1000];
  sprintf(Apn_name,"%s",argv[14]);
  double con_pn=atof(argv[15]);

  double PIXTOAS_PSFFIX=atof(argv[16]);
  double xc=atof(argv[17]);
  double yc=atof(argv[18]);
  double xc_psf=atof(argv[19]); //center of psf map [pix]
  double yc_psf=atof(argv[20]); //center of psf map [pix]
  double rmin=atof(argv[21]); //arcmin
  double rmax=atof(argv[22]); //arcmin
  int    Nbin=atoi(argv[23]);


  //===================================================================================
  //read Afakeit data and fitting
  //=================================================================================== 

  //A= Afakeit
  std::vector<double> r_A_mos1;
  std::vector<double> A_mos1;
  std::vector<double> r_A_mos2;
  std::vector<double> A_mos2;
  std::vector<double> r_A_pn;
  std::vector<double> A_pn;

  readAfakeitdata(Amos1_name,r_A_mos1,A_mos1);
  readAfakeitdata(Amos2_name,r_A_mos2,A_mos2);
  readAfakeitdata(Apn_name,r_A_pn,A_pn);

    

  //===================================================================================
  //read fits data
  //===================================================================================
 
  readimgmap(mos1_psfmap,mos1_psfname, Nxpsf_mos1,Nypsf_mos1,pixtoas_mos1_psf);
  readimgmap(mos1_imgmap,mos1_ratename,Nximg_mos1,Nyimg_mos1,pixtoas_mos1_img);
  readimgmap(mos1_sigmap,mos1_signame, Nxsig_mos1,Nysig_mos1,pixtoas_mos1_sig);

  readimgmap(mos2_psfmap,mos2_psfname, Nxpsf_mos2,Nypsf_mos2,pixtoas_mos2_psf);
  readimgmap(mos2_imgmap,mos2_ratename,Nximg_mos2,Nyimg_mos2,pixtoas_mos2_img);
  readimgmap(mos2_sigmap,mos2_signame, Nxsig_mos2,Nysig_mos2,pixtoas_mos2_sig);

  readimgmap(pn_psfmap,  pn_psfname,   Nxpsf_pn,  Nypsf_pn,  pixtoas_pn_psf);
  readimgmap(pn_imgmap,  pn_ratename,  Nximg_pn,  Nyimg_pn,  pixtoas_pn_img);
  readimgmap(pn_sigmap,  pn_signame,   Nxsig_pn,  Nysig_pn,  pixtoas_pn_sig);

  if(Nximg_mos1 != Nximg_pn   || Nyimg_mos1 != Nyimg_pn   ||
     Nximg_mos2 != Nximg_pn   || Nyimg_mos2 != Nyimg_pn   ||
     Nximg_mos1 != Nximg_mos1 || Nyimg_mos1 != Nyimg_mos2){
    fprintf(stderr, "Different pix number : %s %s %s\n",mos1_ratename,mos2_ratename,pn_ratename);
    return 1;
  }
  if(Nxsig_mos1 != Nxsig_pn   || Nysig_mos1 != Nysig_pn   ||
     Nxsig_mos2 != Nxsig_pn   || Nysig_mos2 != Nysig_pn   ||
     Nxsig_mos1 != Nxsig_mos1 || Nysig_mos1 != Nysig_mos2){
    fprintf(stderr, "Different pix number : %s %s %s\n",mos1_signame,mos2_signame,pn_signame);
    return 1;
  }
  if(Nxpsf_mos1 != Nxpsf_pn   || Nypsf_mos1 != Nypsf_pn   ||
     Nxpsf_mos2 != Nxpsf_pn   || Nypsf_mos2 != Nypsf_pn   ||
     Nxpsf_mos1 != Nxpsf_mos1 || Nypsf_mos1 != Nypsf_mos2){
    fprintf(stderr, "Different pix number : %s %s %s\n",mos1_psfname,mos2_psfname,pn_psfname);
    return 1;
  }
  if(Nximg_mos1 != Nxsig_mos1 || Nyimg_mos1 != Nysig_mos1   ||
     Nximg_mos2 != Nxsig_mos2 || Nyimg_mos2 != Nysig_mos2   ||
     Nximg_pn   != Nxsig_pn   || Nyimg_pn   != Nysig_pn){
    fprintf(stderr, "Different pix number in img and sig fits: %s %s %s %s %s %s\n",
            mos1_ratename,mos2_ratename,pn_ratename,mos1_signame,mos2_signame,pn_signame);
    return 1;
  }
  
  else{
    Nx=Nximg_mos1;  
    Ny=Nyimg_mos1;    
    N=Nx*Ny;
    Nxpsf=Nxpsf_mos1;  
    Nypsf=Nypsf_mos1;
    Npsf=Nxpsf*Nypsf;
    ixpsf_c=Nxpsf_mos1/2-1;
    iypsf_c=Nypsf_mos1/2-1;

    pixtoas=pixtoas_mos1_img;
    pixtoam=pixtoas/60.;
    pixtoas_psf=PIXTOAS_PSFFIX;
    pixtoam_psf=pixtoas_psf/60.;

    //modelmap
    pixtoas_model=PIXTOAS_PSFFIX; //model pix size
    pixtoam_model=pixtoas_model/60.;
    Nxmodel=round(Nx*pixtoas/pixtoas_model);
    Nymodel=round(Ny*pixtoas/pixtoas_model);
    Nmodel=Nxmodel*Nymodel;

    xcmodel=xc*pixtoas/pixtoas_model; //xc in model map
    ycmodel=yc*pixtoas/pixtoas_model; //xc in model map

    xcmodel_psf=xc_psf*pixtoas/pixtoas_model; //psf center-x in model map
    ycmodel_psf=yc_psf*pixtoas/pixtoas_model; //psf center-y in model map

    dconv_max=dconv_max_am/pixtoam_model;  
    pixratio=pixtoas_model/pixtoas;


    //convert r [arcmin] --> r[pix]
    // c1_A*r = (c1_A*pixtoam)*(r/pixtoam)
    //c1_A*=pixtoam;
  }
    


  //========================================================================================
  // setup pix information
  //========================================================================================
  x.resize(Nx);
  y.resize(Ny);
  d.resize(N); 
  
  modelmap_mos1.resize(Nmodel);
  modelmap_mos2.resize(Nmodel);
  modelmap_pn.resize(Nmodel);
  xmodel.resize(Nxmodel);
  ymodel.resize(Nymodel);
  dmodel.resize(Nmodel);  
  dmodel_rateimg.resize(Nmodel);  
  dconv.resize(Nmodel);  

  setup();

  //-------------------------------------------
  //setup wgt for pix at r~>dgrid arcmin
  wgt_mos1=0.;
  wgt_mos2=0.;
  wgt_pn=0.;
  for(i=0;i<Npsf;i++){
    wgt_mos1+=mos1_psfmap[i];
    wgt_mos2+=mos2_psfmap[i];
    wgt_pn+=pn_psfmap[i];
  }
 


  //========================================================================================
  // setup radial profile
  //========================================================================================

   resize_sx(Nbin);

  double r_width=(log10(rmax)-log10(rmin))/Nbin;
 
  int ibin;
   for(ibin=0;ibin<Nbin;ibin++){     
     r1[ibin]=rmin*pow(10,r_width*ibin)/pixtoam;      
     r2[ibin]=rmin*pow(10,r_width*(ibin+1))/pixtoam;   
     r1model[ibin]=rmin*pow(10,r_width*ibin)/pixtoam_model;      
     r2model[ibin]=rmin*pow(10,r_width*(ibin+1))/pixtoam_model;  
     rmid[ibin]=0.5*(r1[ibin]+r2[ibin]);
     rmid_am[ibin]=rmid[ibin]*pixtoam;
     area_model[ibin]=0.; 

     area_mos1[ibin]=EPS; 
     sx_mos1[ibin]=0.; 
     dsx_mos1[ibin]=0.; 

     area_mos2[ibin]=EPS; 
     sx_mos2[ibin]=0.; 
     dsx_mos2[ibin]=0.; 

     area_pn[ibin]=EPS; 
     sx_pn[ibin]=0.; 
     dsx_pn[ibin]=0.; 

     sx_mos10[ibin]=0.; 
     dsx_mos10[ibin]=0.; 

     sx_mos20[ibin]=0.; 
     dsx_mos20[ibin]=0.; 

     sx_pn0[ibin]=0.; 
     dsx_pn0[ibin]=0.; 

     sx_ave[ibin]=0.;
     dsx_sys[ibin]=0.;

     sx_ave0[ibin]=0.;
     dsx_sys0[ibin]=0.;

   }

 //-----------------------------------------------------------------------
  //setup model map
  for(iy=0;iy<Nymodel;iy++){  
  for(ix=0;ix<Nxmodel;ix++){
      i=ix+Nxmodel*iy;
      dmodel[i]=sqrt((xmodel[ix]-xcmodel)*(xmodel[ix]-xcmodel)+(ymodel[iy]-ycmodel)*(ymodel[iy]-ycmodel)); //unit of psfmap pix size
      dmodel_rateimg[i]=dmodel[i]*pixratio; //unit of rateimage pix size

      dconv[i]=sqrt((1.*ix-xcmodel_psf)*(1.*ix-xcmodel_psf)+(1.*iy-ycmodel_psf)*(1.*iy-ycmodel_psf)); //distance from center to make PSF map 
      for(ibin=0;ibin<Nbin;ibin++){
	if(dmodel[i]>r1model[ibin] &&  dmodel[i]<=r2model[ibin]){
	  area_model[ibin]+=1.;
	}
        }
  }}
 

  //Calculating Sx from the data
  int ix,iy;
  for(iy=0;iy<Ny;iy++){  
  for(ix=0;ix<Nx;ix++){
      i=ix+Nx*iy;
      d[i]=sqrt((x[ix]-xc)*(x[ix]-xc)+(y[iy]-yc)*(y[iy]-yc)); //distance from the selected center

      if(fabs(mos1_imgmap[i])>rate_th){
      for(ibin=0;ibin<Nbin;ibin++){
	if(d[i]>r1[ibin] &&  d[i]<=r2[ibin] ){
	sx_mos1[ibin]+=mos1_imgmap[i];
	dsx_mos1[ibin]+=mos1_sigmap[i]*mos1_sigmap[i];
	area_mos1[ibin]+=1.;
        break;
      }}
      }

      if(fabs(mos2_imgmap[i])>rate_th){
      for(ibin=0;ibin<Nbin;ibin++){
	if(d[i]>r1[ibin] &&  d[i]<=r2[ibin] ){
	sx_mos2[ibin]+=mos2_imgmap[i];
	dsx_mos2[ibin]+=mos2_sigmap[i]*mos2_sigmap[i];
	area_mos2[ibin]+=1.;
        break;
      }}
      }


      if(fabs(pn_imgmap[i])>rate_th){
      for(ibin=0;ibin<Nbin;ibin++){
	if(d[i]>r1[ibin] &&  d[i]<=r2[ibin] ){
	sx_pn[ibin]+=pn_imgmap[i];
	dsx_pn[ibin]+=pn_sigmap[i]*pn_sigmap[i];
	area_pn[ibin]+=1.;
        break;
      }}
      }

   }}





  for(ibin=0;ibin<Nbin;ibin++){
    /*
    if(area_mos1[ibin]<=1){
      fprintf(stderr,"ERROR : the number of sampling pix <= 1 at %d th bin\n",ibin);
      fprintf(stderr,"        Change Nbin %s \n",mos1_ratename);

      return 0; 
    } 
    if(area_mos2[ibin]<=1){
      fprintf(stderr,"ERROR : the number of sampling pix <= 1 at %d th bin\n",ibin);
      fprintf(stderr,"        Change Nbin %s \n",mos2_ratename);

      return 0; 
    } 
    if(area_pn[ibin]<=1){
      fprintf(stderr,"ERROR : the number of sampling pix <= 1 at %d th bin\n",ibin);
      fprintf(stderr,"        Change Nbin %s \n",pn_ratename);

      return 0; 
    }
    */
       sx_mos10[ibin]=sx_mos1[ibin]/area_mos1[ibin];
       dsx_mos10[ibin]=sqrt(dsx_mos1[ibin]/area_mos1[ibin]/(area_mos1[ibin]-1.));

       sx_mos20[ibin]=sx_mos2[ibin]/area_mos2[ibin];
       dsx_mos20[ibin]=sqrt(dsx_mos2[ibin]/area_mos2[ibin]/(area_mos2[ibin]-1.));

       sx_pn0[ibin]=sx_pn[ibin]/area_pn[ibin];
       dsx_pn0[ibin]=sqrt(dsx_pn[ibin]/area_pn[ibin]/(area_pn[ibin]-1.));

       sx_mos1[ibin]=sx_mos1[ibin]/(con_mos1*A_mos1[0])/area_mos1[ibin];
       dsx_mos1[ibin]=sqrt(dsx_mos1[ibin]/area_mos1[ibin]/(area_mos1[ibin]-1.))/(con_mos1*A_mos1[0]);

       sx_mos2[ibin]=sx_mos2[ibin]/(con_mos2*A_mos2[0])/area_mos2[ibin];
       dsx_mos2[ibin]=sqrt(dsx_mos2[ibin]/area_mos2[ibin]/(area_mos2[ibin]-1.))/(con_mos2*A_mos2[0]);

       sx_pn[ibin]=sx_pn[ibin]/(con_pn*A_pn[0])/area_pn[ibin];
       dsx_pn[ibin]=sqrt(dsx_pn[ibin]/area_pn[ibin]/(area_pn[ibin]-1.))/(con_pn*A_pn[0]);

       sx_ave[ibin]=(sx_mos1[ibin]+sx_mos2[ibin]+sx_pn[ibin])/3.;
       //       sx_ave0[ibin]=(sx_mos10[ibin]+sx_mos20[ibin]+sx_pn0[ibin])/3.;

     if( sx_ave[ibin]==0){
	 r1[ibin]=0.;
	 r2[ibin]=0.;
	 rmid[ibin]=0.;
       }

  }
  

   sx_mos10.erase(std::remove(sx_mos10.begin(), sx_mos10.end(), 0), sx_mos10.end());
   sx_mos20.erase(std::remove(sx_mos20.begin(), sx_mos20.end(), 0), sx_mos20.end());
   sx_pn0.erase(std::remove(sx_pn0.begin(), sx_pn0.end(), 0), sx_pn0.end());
   sx_mos1.erase(std::remove(sx_mos1.begin(), sx_mos1.end(), 0), sx_mos1.end());
   sx_mos2.erase(std::remove(sx_mos2.begin(), sx_mos2.end(), 0), sx_mos2.end());
   sx_pn.erase(std::remove(sx_pn.begin(), sx_pn.end(), 0), sx_pn.end());
   sx_ave.erase(std::remove(sx_ave.begin(), sx_ave.end(), 0), sx_ave.end());

   dsx_mos10.erase(std::remove(dsx_mos10.begin(), dsx_mos10.end(), 0), dsx_mos10.end());
   dsx_mos20.erase(std::remove(dsx_mos20.begin(), dsx_mos20.end(), 0), dsx_mos20.end());
   dsx_pn0.erase(std::remove(dsx_pn0.begin(), dsx_pn0.end(), 0), dsx_pn0.end());
   dsx_mos1.erase(std::remove(dsx_mos1.begin(), dsx_mos1.end(), 0), dsx_mos1.end());
   dsx_mos2.erase(std::remove(dsx_mos2.begin(), dsx_mos2.end(), 0), dsx_mos2.end());
   dsx_pn.erase(std::remove(dsx_pn.begin(), dsx_pn.end(), 0), dsx_pn.end());

   r1.erase(std::remove(r1.begin(), r1.end(), 0), r1.end());
   r2.erase(std::remove(r2.begin(), r2.end(), 0), r2.end());
   rmid.erase(std::remove(rmid.begin(), rmid.end(), 0), rmid.end());

   Nbin=sx_mos10.size();  
   resize_sx(Nbin);

    char sx_name[1000];
    sprintf(sx_name, "Sxresult/sx_profile_noncor.dat");

    char sx_name2[1000];
    sprintf(sx_name2, "Sxresult/sx_profile.dat");


    FILE *fp_sx = fopen( sx_name, "w" );
    fprintf(fp_sx,"# r[pix] dr_m[pix] dr_p[pix] r[arcmin] dr_m[arcmin] dr_p[arcmin] sx_mos1 dsx_mos1 dsx_mos1_stat sx_mos2 dsx_mos2 dsx_mos2_stat sx_pn dsx_pn dsx_pn_stat\n");

    FILE *fp_sx2 = fopen( sx_name2, "w" );
    fprintf(fp_sx2,"# r[pix] dr_m[pix] dr_p[pix] r[arcmin] dr_m[arcmin] dr_p[arcmin] sx_mos1 dsx_mos1 dsx_mos1_stat sx_mos2 dsx_mos2 dsx_mos2_stat sx_pn dsx_pn dsx_pn_stat sx_all dsx_all dsx_all_stat\n");


  for(ibin=0;ibin<Nbin;ibin++){
    dsx_sys[ibin]=sqrt((pow(sx_mos1[ibin]-sx_ave[ibin],2)+pow(sx_mos2[ibin]-sx_ave[ibin],2)+pow(sx_pn[ibin]-sx_ave[ibin],2))/3./2.);
    sx_sig_stat[ibin]=sqrt((dsx_mos1[ibin]*dsx_mos1[ibin]+dsx_mos2[ibin]*dsx_mos2[ibin]+dsx_pn[ibin]*dsx_pn[ibin])/3.);
    sx_sig[ibin]=sqrt(sx_sig_stat[ibin]*sx_sig_stat[ibin]+dsx_sys[ibin]*dsx_sys[ibin]);
    dsx_mos1[ibin]=sqrt(dsx_mos1[ibin]*dsx_mos1[ibin]+dsx_sys[ibin]*dsx_sys[ibin]);
    dsx_mos2[ibin]=sqrt(dsx_mos2[ibin]*dsx_mos2[ibin]+dsx_sys[ibin]*dsx_sys[ibin]);
    dsx_pn[ibin]=sqrt(dsx_pn[ibin]*dsx_pn[ibin]+dsx_sys[ibin]*dsx_sys[ibin]);

    

    dsx_sys0[ibin]=sqrt((pow(sx_mos10[ibin]-sx_ave[ibin]*Amos1_norm*con_mos1,2)+pow(sx_mos20[ibin]-sx_ave[ibin]*Amos2_norm*con_mos2,2)+pow(sx_pn0[ibin]-sx_ave[ibin]*Apn_norm*con_pn,2))/3./2.);
    dsx_mos10[ibin]=sqrt(dsx_mos10[ibin]*dsx_mos10[ibin]+dsx_sys0[ibin]*dsx_sys0[ibin]);
    dsx_mos20[ibin]=sqrt(dsx_mos20[ibin]*dsx_mos20[ibin]+dsx_sys0[ibin]*dsx_sys0[ibin]);
    dsx_pn0[ibin]=sqrt(dsx_pn0[ibin]*dsx_pn0[ibin]+dsx_sys0[ibin]*dsx_sys0[ibin]);


      fprintf(fp_sx,"%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",rmid[ibin],r1[ibin]-rmid[ibin],r2[ibin]-rmid[ibin], rmid[ibin]*pixtoam,(r1[ibin]-rmid[ibin])*pixtoam,(r2[ibin]-rmid[ibin])*pixtoam,sx_mos10[ibin],dsx_mos10[ibin],sqrt(dsx_mos10[ibin]*dsx_mos10[ibin]-dsx_sys0[ibin]*dsx_sys0[ibin]),sx_mos20[ibin],dsx_mos20[ibin],sqrt(dsx_mos20[ibin]*dsx_mos20[ibin]-dsx_sys0[ibin]*dsx_sys0[ibin]),sx_pn0[ibin],dsx_pn0[ibin],sqrt(dsx_pn0[ibin]*dsx_pn0[ibin]-dsx_sys0[ibin]*dsx_sys0[ibin]));

      fprintf(fp_sx2,"%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",rmid[ibin],r1[ibin]-rmid[ibin],r2[ibin]-rmid[ibin], rmid[ibin]*pixtoam,(r1[ibin]-rmid[ibin])*pixtoam,(r2[ibin]-rmid[ibin])*pixtoam,sx_mos1[ibin],dsx_mos1[ibin],sqrt(dsx_mos1[ibin]*dsx_mos1[ibin]-dsx_sys[ibin]*dsx_sys[ibin]),sx_mos2[ibin],dsx_mos2[ibin],sqrt(dsx_mos2[ibin]*dsx_mos2[ibin]-dsx_sys[ibin]*dsx_sys[ibin]),sx_pn[ibin],dsx_pn[ibin],sqrt(dsx_pn[ibin]*dsx_pn[ibin]-dsx_sys[ibin]*dsx_sys[ibin]),sx_ave[ibin],sx_sig[ibin],sx_sig_stat[ibin]);

   }
     fclose(fp_sx);
     fclose(fp_sx2);
  
     std::vector<double>().swap(mos1_imgmap);
     std::vector<double>().swap(mos1_sigmap);
     std::vector<double>().swap(mos2_imgmap);
     std::vector<double>().swap(mos2_sigmap);
     std::vector<double>().swap(pn_imgmap);
     std::vector<double>().swap(pn_sigmap);

   //========================================================================================
   // fitting
   //========================================================================================
     
    Chisq_nonconv fFCN(rmid,sx_mos1,dsx_mos1,sx_mos2,dsx_mos2,sx_pn,dsx_pn);
    fprintf(stdout,"fitting without convolution...\n");
    int Npara=9;

  
    double norm1= 1.;
    double r01= 1.;
    double alpha=ALPHA;
    double beta1= 1.2;
    double norm2=0.4;
    double r02= 14.;
    double beta2= 0.6;
    double epsilon=EPSILON;    
    double sb_mos1=0.1;
    double sb_mos2=0.1;
    double sb_pn=0.1;


    MnUserParameters upar;
    upar.Add("norm1", norm1, 0.1);
    upar.Add("r01", r01, 0.1);
    upar.Add("beta1", beta1, 0.2);
    upar.Add("norm2", norm2, 0.1);
    upar.Add("r02", r02, 2.);
    upar.Add("beta2", beta2, 0.2);
    upar.Add("sb_mos1", sb_mos1, 0.01);
    upar.Add("sb_mos2", sb_mos2, 0.01);
    upar.Add("sb_pn", sb_pn, 0.01);
    
    //set lower limit
    upar.SetLowerLimit("norm1", 0.0);
    upar.SetLowerLimit("r01", 0.0);
    upar.SetLimits("beta1", 0.0,3.0);
    upar.SetLowerLimit("norm2", 0.0);
    upar.SetLimits("r02", 0.0,rmax/pixtoam);
    upar.SetLimits("beta2", 0.0,3.0);
    upar.SetLowerLimit("sb_mos1", 0.0);
    upar.SetLowerLimit("sb_mos2", 0.0);
    upar.SetLowerLimit("sb_pn", 0.0);
    
    // create MIGRAD minimizer
    MnMigrad migrad(fFCN, upar);

    // Minimize
    FunctionMinimum min = migrad();

    // output
    std::cout<<"minimum: "<<min<<std::endl;
     
    norm1= min.UserState().Value("norm1");
    r01= min.UserState().Value("r01");
    //alpha= min.UserState().Value("alpha");
    beta1= min.UserState().Value("beta1");
    norm2= min.UserState().Value("norm2");
    r02= min.UserState().Value("r02");
    beta2= min.UserState().Value("beta2");
    //epsilon= min.UserState().Value("epsilon");
    sb_mos1= min.UserState().Value("sb_mos1");
    sb_mos2= min.UserState().Value("sb_mos2");
    sb_pn= min.UserState().Value("sb_pn");

    //---------------------------------------------------------------------------------
    //second fitting with conv 
    Chisq fFCN2(rmid,sx_mos1,dsx_mos1,sx_mos2,dsx_mos2,sx_pn,dsx_pn);
    fprintf(stdout,"fitting with convolution...\n");
    
    MnUserParameters upar2;
    upar2.Add("norm1", norm1, 200.);
    upar2.Add("r01", r01, 0.1);
    //upar2.Add("alpha", alpha, 0.1);
    upar2.Add("beta1", beta1, 0.4);
    upar2.Add("norm2", norm2, 0.1);
    upar2.Add("r02", r02, 2.);
    upar2.Add("beta2", beta2, 0.4);
    //upar2.Add("epsilon", epsilon, 0.4);
    upar2.Add("sb_mos1", sb_mos1, 0.1);
    upar2.Add("sb_mos2", sb_mos2, 0.1);
    upar2.Add("sb_pn", sb_pn, 0.1);    

    //set lower limit
    upar2.SetLowerLimit("norm1", 0.0);
    upar2.SetLowerLimit("r01", 0.0);
    //upar2.SetLowerLimit("alpha", 0.);
    upar2.SetLowerLimit("beta1", 0.0);
    //upar2.SetLimits("norm2", 0.0,1.0);
    upar2.SetLowerLimit("norm2", 0.0);
    upar2.SetLowerLimit("r02", 0.0);
    upar2.SetLowerLimit("beta2", 0.0);
    //  upar2.SetUpperLimit("epsilon", 5.0);
    //upar2.SetLimits("epsilon", -1.0, 5.0);
    upar2.SetLowerLimit("sb_mos1", 0.0);
    upar2.SetLowerLimit("sb_mos2", 0.0);
    upar2.SetLowerLimit("sb_pn", 0.0);

    // create MIGRAD minimizer
    MnMigrad migrad2(fFCN2, upar2);
    
    // Minimize
    FunctionMinimum min2 = migrad2();

    // output
    std::cout<<"minimum: "<<min2<<std::endl;


    migrad2.RemoveLimits("norm1");
    migrad2.RemoveLimits("r01");
    //  migrad2.RemoveLimits("alpha");
    migrad2.RemoveLimits("beta1");
    migrad2.RemoveLimits("norm2");
    migrad2.RemoveLimits("r02");
    migrad2.RemoveLimits("beta2");
    migrad2.RemoveLimits("sb_mos1");
    migrad2.RemoveLimits("sb_mos2");
    migrad2.RemoveLimits("sb_pn");

    // Minimize
    FunctionMinimum min3 = migrad2();

    // output
    std::cout<<"minimum: "<<min3<<std::endl;


    norm1= min3.UserState().Value("norm1");
    r01= min3.UserState().Value("r01");
    //alpha= min3.UserState().Value("alpha");
    beta1= min3.UserState().Value("beta1");
    norm2= min3.UserState().Value("norm2");
    r02= min3.UserState().Value("r02");
    beta2= min3.UserState().Value("beta2");
    //epsilon= min3.UserState().Value("epsilon");
    sb_mos1= min3.UserState().Value("sb_mos1");
    sb_mos2= min3.UserState().Value("sb_mos2");
    sb_pn= min3.UserState().Value("sb_pn");


    ROOT::Minuit2::MnUserCovariance covar=min3.UserCovariance();
    double Cov[Npara][Npara];
    double rCov[Npara][Npara];
    for(i=0;i<Npara;i++){
    for(j=0;j<Npara;j++){
    Cov[i][j]=covar(i,j); 
    rCov[i][j]=covar(i,j)/sqrt(covar(i,i)*covar(j,j));
    }
    }

    double chisq=min3.Fval();
    int Ndof=3*Nbin-Npara;
    double redchisq=chisq/Ndof; 
    double prob=gsl_sf_gamma_inc_Q(Ndof/2.,chisq/2.);

    //=================================================================================================================
    char *now;
    tm *newtime;
    time_t aclock;
    time(&aclock);
    newtime = localtime(&aclock);  
    now     = asctime(newtime);    
    
    char log_name[1000];
    sprintf(log_name, "Sxresult/sx_bestfit.log");
    FILE *fp_log = fopen( log_name, "w" );

    fprintf(fp_log, "#------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "# LOG %s", now);
    fprintf(fp_log, "#------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "   xc           [pix]                           : %.8e \n",xc);
    fprintf(fp_log, "   yc           [pix]                           : %.8e \n",yc);
    fprintf(fp_log, "   rmin         [arcmin]                        : %.8e \n",rmin);
    fprintf(fp_log, "   rmax         [arcmin]                        : %.8e \n",rmax);
    fprintf(fp_log, "   Nbin                                         : %.d \n",Nbin);
    fprintf(fp_log, "#------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "# fits image \n");
    fprintf(fp_log, "#------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, " mos1 rate image : %s\n", mos1_ratename);
    fprintf(fp_log, " mos1 sig  image : %s\n", mos1_signame);
    fprintf(fp_log, " mos1 psf  image : %s\n", mos1_psfname);
    fprintf(fp_log, " mos2 rate image : %s\n", mos2_ratename);
    fprintf(fp_log, " mos2 sig  image : %s\n", mos2_signame);
    fprintf(fp_log, " mos2 psf  image : %s\n", mos2_psfname);
    fprintf(fp_log, " pn   rate image : %s\n", pn_ratename);
    fprintf(fp_log, " pn   sig  image : %s\n", pn_signame);
    fprintf(fp_log, " pn   psf  image : %s\n", pn_psfname);
    fprintf(fp_log, "#------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "# Emissivity \n");
    fprintf(fp_log, "#------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "   mos1 file  : %s \n",Amos1_name);
    fprintf(fp_log, "   mos2 file  : %s \n",Amos2_name);
    fprintf(fp_log, "   pn   file  : %s \n",Apn_name);
    fprintf(fp_log, "   polynomial function \n");
    fprintf(fp_log, "   Lambda(r)=Lamba_0*(1+c1_A*r) r=pix unit\n");
    fprintf(fp_log, "   Lamba_mos1_0                                    : %.8e \n",Amos1_norm);
    fprintf(fp_log, "   Lamba_mos2_0                                    : %.8e \n",Amos2_norm);
    fprintf(fp_log, "   Lamba_pn_0                                      : %.8e \n",Apn_norm);
    fprintf(fp_log, "   c1_A                                            : %.8e \n",c1_A);
    //  fprintf(fp_log, "   c2_A                                            : %.8e \n",c2_A);

    fprintf(fp_log, "#------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "   model nenp= norm1*(r/r01)^-alpha/((1+(r/r01)^2)^{3beta1-alpha/2}/(1.+(r/r01)**gamma)**(epsilon/gamma) + norm2*(1+(r/r02)^2)^{3beta2})\n");
    fprintf(fp_log, "   Npara        (number of parameters)             : %d \n",Npara);
    fprintf(fp_log, "  0:  norm1        =n01^2 for mos1                 : %.8e \n",norm1);
    fprintf(fp_log, "  1:  r01          [pix]                           : %.8e \n",r01);
    fprintf(fp_log, "  1:  r01          [arcmin]                        : %.8e \n",r01*pixtoam);
    fprintf(fp_log, "      alpha(fix)                                   : %.8e \n",alpha);
    fprintf(fp_log, "  2:  beta1                                        : %.8e \n",beta1);
    fprintf(fp_log, "  3:  norm2        =n02^2                          : %.8e \n",norm2);
    fprintf(fp_log, "  4:  r02          [pix]                           : %.8e \n",r02);
    fprintf(fp_log, "  4:  r02          [arcmin]                        : %.8e \n",r02*pixtoam);
    fprintf(fp_log, "  5:  beta2                                        : %.8e \n",beta2);
    fprintf(fp_log, "      gamma(fix)                                   : %.8e \n",gamma1);
    fprintf(fp_log, "      epsilon(fix)                                 : %.8e \n",epsilon);
    fprintf(fp_log, "  6:  sb_mos1                                      : %.8e \n",sb_mos1);
    fprintf(fp_log, "  7:  sb_mos2                                      : %.8e \n",sb_mos2);
    fprintf(fp_log, "  8:  sb_pn                                        : %.8e \n",sb_pn);
    fprintf(fp_log, "#------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "   Covariance\n");
    for(i=0;i<Npara;i++){
    for(j=0;j<Npara;j++){
      fprintf(fp_log, "%.8e ",Cov[i][j]);
    }
      fprintf(fp_log, "\n");
    }
    fprintf(fp_log, "#------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "   Coefficient\n");
    for(i=0;i<Npara;i++){
    for(j=0;j<Npara;j++){
      fprintf(fp_log, "%.8e ",rCov[i][j]);
    }
      fprintf(fp_log, "\n");
    }
    fprintf(fp_log, "#------------------------------------------------------------------------------------------\n");
    fprintf(fp_log, "   chisq                                        : %.5lf \n",chisq);
    fprintf(fp_log, "   d.o.f                                        : %d \n",Ndof);
    fprintf(fp_log, "   redchisq                                     : %.5lf \n",redchisq);
    fprintf(fp_log, "   Probability Q                                : %.8e \n",prob);
    

    //=================================================================================================================
    //
    // best fit
    //=================================================================================================================
    Sx_2beta Sx_mos1(norm1,r01,beta1, norm2,r02,beta2);
    n_2beta n_Vikhlinin0(norm1,r01,beta1, norm2,r02,beta2);
    n2_2beta n2_Vikhlinin0(norm1,r01,beta1, norm2,r02,beta2);
  
  for(ispline=0;ispline<Nspline;ispline++){
    Sxspline_mos1[ispline]=Sx_mos1(rspline[ispline]);  
  }

  s_mos1.set_points(rspline,Sxspline_mos1);  

  //!!! setup 2d map 
  //initialize
  for(ix=0;ix<Nxmodel;ix++){
  for(iy=0;iy<Nymodel;iy++){ 
    i=ix+Nxmodel*iy;
    modelmap_mos1[i]=s_mos1(dmodel_rateimg[i]);
    modelmap_mos2[i]=modelmap_mos1[i];
    modelmap_pn[i]=modelmap_mos1[i];
  }
  }

  //convolution
  conv(modelmap_mos1,modelmap_mos2,modelmap_pn,r1model,r2model,area_model,Sx_model_conv_mos1,Sx_model_conv_mos2,Sx_model_conv_pn);

    //---------------------------------------------------
    // output
    char output_name[1000];
    sprintf(output_name, "Sxresult/sx_bestfit.dat");
  
    FILE *fp_output = fopen( output_name, "w" );
    fprintf(fp_output,"# r[pix] r[arcmin] sx_conv+sb(mos1) sx_conv+sb(mos2) sx_conv+sb(pn) sb_mos1 sb_mos2 sb_pn\n");

    char output_name0[1000];
    sprintf(output_name0, "Sxresult/sx_bestfit_noncor.dat");
  
    FILE *fp_output0 = fopen( output_name0, "w" );
    fprintf(fp_output0,"# r[pix] r[arcmin] sx_conv+sb(mos1) sx_conv+sb(mos2) sx_conv+sb(pn)\n");

    char output_name1[1000];
    sprintf(output_name1, "Sxresult/sx_bestfit_bkgsub.dat");
  
    FILE *fp_output1 = fopen( output_name1, "w" );
    fprintf(fp_output1,"# r[pix] r[arcmin] sx_conv(mos1) sx_conv(mos2) sx_conv(pn)\n");

    char sx_name3[1000];
    sprintf(sx_name3, "Sxresult/sx_profile_bkgsub.dat");

    FILE *fp_sx3 = fopen( sx_name3, "w" );
    fprintf(fp_sx3,"# r[pix] dr_m[pix] dr_p[pix] r[arcmin] dr_m[arcmin] dr_p[arcmin] sx_mos1 dsx_mos1 dsx_mos1_stat sx_mos2 dsx_mos2 dsx_mos2_stat sx_pn dsx_pn dsx_pn_stat sx_all dsx_all dsx_all_stat\n");

  
   for(ibin=0;ibin<Nbin;ibin++){

     //bkgsub	
      fprintf(fp_output1,"%.8e %.8e %.8e %.8e %.8e\n",
                          rmid[ibin],rmid[ibin]*pixtoam,
                          Sx_model_conv_mos1[ibin],
                          Sx_model_conv_mos2[ibin],
                          Sx_model_conv_pn[ibin]);
   
      //bkgsub data profile
      fprintf(fp_sx3,"%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
                     rmid[ibin],r1[ibin]-rmid[ibin],r2[ibin]-rmid[ibin], 
                     rmid[ibin]*pixtoam,(r1[ibin]-rmid[ibin])*pixtoam,(r2[ibin]-rmid[ibin])*pixtoam,
	             sx_mos1[ibin]-sb_mos1,dsx_mos1[ibin],sqrt(dsx_mos1[ibin]*dsx_mos1[ibin]-dsx_sys[ibin]*dsx_sys[ibin]),
	             sx_mos2[ibin]-sb_mos2,dsx_mos2[ibin],sqrt(dsx_mos2[ibin]*dsx_mos2[ibin]-dsx_sys[ibin]*dsx_sys[ibin]),
	             sx_pn[ibin]-sb_pn,dsx_pn[ibin],sqrt(dsx_pn[ibin]*dsx_pn[ibin]-dsx_sys[ibin]*dsx_sys[ibin]),
                     sx_ave[ibin]-(sb_mos1+sb_mos2+sb_pn)/3.,sx_sig[ibin],sx_sig_stat[ibin]);

      //noncor
      fprintf(fp_output0,"%.8e %.8e %.8e %.8e %.8e\n",
                          rmid[ibin],rmid[ibin]*pixtoam,
	                  (Sx_model_conv_mos1[ibin]+sb_mos1)*con_mos1*Amos1_norm,
	                  (Sx_model_conv_mos2[ibin]+sb_mos2)*con_mos2*Amos2_norm,
	                  (Sx_model_conv_pn[ibin]+sb_pn)*con_pn*Apn_norm);

      //corr
     Sx_model_conv_mos1[ibin]+=sb_mos1;
     Sx_model_conv_mos2[ibin]+=sb_mos2;
     Sx_model_conv_pn[ibin]+=sb_pn;
      
      fprintf(fp_output,"%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
                        rmid[ibin],rmid[ibin]*pixtoam,
                        Sx_model_conv_mos1[ibin],
                        Sx_model_conv_mos2[ibin],
                        Sx_model_conv_pn[ibin],
                        sb_mos1,
                        sb_mos2,
	                sb_pn);
      
    }
   fclose(fp_output);
   fclose(fp_output0);
   fclose(fp_output1);
   fclose(fp_sx3);
    //---------------------------------------------------
    // output
    //---------------------------------------------------
    char output_name2[1000];
    sprintf(output_name2, "Sxresult/ne_bestfit.dat");

    FILE *fp_output2 = fopen( output_name2, "w" );
    fprintf(fp_output2,"# r[pix] r[arcmin] n=(nenp)**0.5 n_m n_p  n2=(nenp) n2_m n2_p\n");


    int Nline=500;
    int iline;
    double rline_min=0.001; //arcmin
    double rline_max=30.;  //arcmin
    double rline_width=(log(rline_max)-log(rline_min))/Nline;
    std::vector<double> rline(Nline);  
    std::vector<double> nline(Nline);  
    std::vector<double> dnline(Nline);  
    std::vector<double> n2line(Nline);  
    std::vector<double> dn2line(Nline);  


    //error calculation
    gsl_matrix *cov=gsl_matrix_alloc (Npara, Npara);
    for(i=0;i<Npara;i++){
     for(j=0;j<Npara;j++){
     gsl_matrix_set (cov, i, j, Cov[i][j]);
    }}
    gsl_linalg_cholesky_decomp(cov);


    int iboot;
    int Nboot=10000;
    double norm1_tmp;
    double r01_tmp;
    double alpha_tmp;
    double beta1_tmp;
    double norm2_tmp;
    double r02_tmp;
    double beta2_tmp;
    double epsilon_tmp;


    for(iline=0;iline<Nline;iline++){     
     rline[iline]=rline_min*exp(rline_width*iline)/pixtoam;      
     nline[iline]=n_Vikhlinin0(rline[iline]);
     dnline[iline]=0;
     n2line[iline]=n2_Vikhlinin0(rline[iline]);
     dn2line[iline]=0;
   }

      gsl_rng *gBaseRand; 
      unsigned long randSeed;
      gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);
      srand(time(NULL));                    /* initialization for rand() */
      randSeed = rand();                    /* returns a non-negative integer */
      gsl_rng_set (gBaseRand, randSeed);    /* seed the PRNG */

      double gauss0,gauss1,gauss2,gauss3,gauss4,gauss5,gauss6,gauss7,gauss8;
   
  for(iboot=0;iboot<Nboot;iboot++){
    do{

    gauss0=gsl_ran_gaussian(gBaseRand,1.);
    gauss1=gsl_ran_gaussian(gBaseRand,1.);
    gauss2=gsl_ran_gaussian(gBaseRand,1.);
    gauss3=gsl_ran_gaussian(gBaseRand,1.);
    gauss4=gsl_ran_gaussian(gBaseRand,1.);
    gauss5=gsl_ran_gaussian(gBaseRand,1.);
    gauss6=gsl_ran_gaussian(gBaseRand,1.);
    gauss7=gsl_ran_gaussian(gBaseRand,1.);
    

    norm1_tmp=gsl_matrix_get (cov, 0, 0)*gauss0+norm1;
    r01_tmp  =gsl_matrix_get (cov, 1, 0)*gauss0+gsl_matrix_get (cov, 1, 1)*gauss1+r01;
    beta1_tmp=gsl_matrix_get (cov, 2, 0)*gauss0+gsl_matrix_get (cov, 2, 1)*gauss1+gsl_matrix_get (cov, 2, 2)*gauss2+beta1;
    norm2_tmp=gsl_matrix_get (cov, 3, 0)*gauss0+gsl_matrix_get (cov, 3, 1)*gauss1+gsl_matrix_get (cov, 3, 2)*gauss2+gsl_matrix_get (cov, 3, 3)*gauss3+norm2;
    r02_tmp  =gsl_matrix_get (cov, 4, 0)*gauss0+gsl_matrix_get (cov, 4, 1)*gauss1+gsl_matrix_get (cov, 4, 2)*gauss2+gsl_matrix_get (cov, 4, 3)*gauss3+gsl_matrix_get (cov, 4, 4)*gauss4+r02;
    beta2_tmp=gsl_matrix_get (cov, 5, 0)*gauss0+gsl_matrix_get (cov, 5, 1)*gauss1+gsl_matrix_get (cov, 5, 2)*gauss2+gsl_matrix_get (cov, 5, 3)*gauss3+gsl_matrix_get (cov, 5, 4)*gauss4+gsl_matrix_get (cov, 5, 5)*gauss5+beta2;
    
    }while(r01_tmp<0. || r02_tmp<0. || norm1_tmp<0 || norm2_tmp<0);

    
    n_2beta n_Vikhlinin1(norm1_tmp,r01_tmp,beta1_tmp, norm2_tmp,r02_tmp,beta2_tmp);
    n2_2beta n2_Vikhlinin1(norm1_tmp,r01_tmp,beta1_tmp, norm2_tmp,r02_tmp,beta2_tmp);

     for(iline=0;iline<Nline;iline++){   
       dnline[iline]+=pow(n_Vikhlinin1(rline[iline])-nline[iline],2);
       dn2line[iline]+=pow(n2_Vikhlinin1(rline[iline])-n2line[iline],2);
     }
  }

   for(iline=0;iline<Nline;iline++){     
     dnline[iline]=sqrt(dnline[iline]/(Nboot-1.)); 
     dn2line[iline]=sqrt(dn2line[iline]/(Nboot-1.)); 
     fprintf(fp_output2,"%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",rline[iline],rline[iline]*pixtoam,
                                                               	    nline[iline],nline[iline]-dnline[iline],nline[iline]+dnline[iline],
	     n2line[iline],n2line[iline]-dn2line[iline],n2line[iline]+dn2line[iline]);
   }

   fclose(fp_output2);


}




void resize_sx(int Nbin)
{
  r1.resize(Nbin);
  r2.resize(Nbin);
  rmid.resize(Nbin);
  rmid_am.resize(Nbin);
  sx_mos1.resize(Nbin); 
  dsx_mos1.resize(Nbin); 
  area_mos1.resize(Nbin);
  sx_mos2.resize(Nbin); 
  dsx_mos2.resize(Nbin); 
  area_mos2.resize(Nbin);
  sx_pn.resize(Nbin); 
  dsx_pn.resize(Nbin); 
  area_pn.resize(Nbin);
  sx_ave.resize(Nbin); 
  dsx_sys.resize(Nbin);
  sx_sig.resize(Nbin); 
  sx_sig_stat.resize(Nbin); 

  sx_mos10.resize(Nbin); 
  dsx_mos10.resize(Nbin); 
  sx_mos20.resize(Nbin); 
  dsx_mos20.resize(Nbin); 
  sx_pn0.resize(Nbin); 
  dsx_pn0.resize(Nbin); 
  sx_ave0.resize(Nbin); 
  dsx_sys0.resize(Nbin);

  r1model.resize(Nbin);
  r2model.resize(Nbin);
  area_model.resize(Nbin); 
  Sx_model_mos1.resize(Nbin); 
  Sx_model_mos2.resize(Nbin); 
  Sx_model_pn.resize(Nbin); 
  Sx_model_conv_mos1.resize(Nbin); 
  Sx_model_conv_mos2.resize(Nbin); 
  Sx_model_conv_pn.resize(Nbin); 

}

