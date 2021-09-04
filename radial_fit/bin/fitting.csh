#! /bin/csh -f

#----------------------------------------------------
#vikhlinin model fixed parameter
set acool=1
set at=0
set bt=2
# gamma is a free parameter
#---------------------------------------------------


if ( $#argv != 5 ) then
cat <<EOF

bin/fit_sx.csh name z rinn rout Nbin

   name : output cluster name
   z    : cluster redshift
   rinn : innermost radius [arcmin]
   rout : outermost radius [arcmin]
   Nbin : num of bins
   
   comoslogy : H0=70 km/s/Mpc Om0=0.3 OL=0.7 

Use:  bin/fit_sx.cs ABELL1689 0.1832 0.1 10 22


EOF

exit

endif
#-------------------------------------------------------------------
set name=$1
set z=$2
set rinn=$3 #innermost radius [arcmin]
set rout=$4  #outermost radius [arcmin]
set Nbin=$5  #num of bins
#set rinn=0.1 #innermost radius [arcmin]
#set rout=10  #outermost radius [arcmin]
#set Nbin=22  #num of bins
#--------------------------------------------------------------------

setenv OMP_NUM_THREADS 4

set xc=`gawk '{if($1=="PixX,PixY") print $3}' ../spec_ana/log/Center.log` # x center [pixel]
set yc=`gawk '{if($1=="PixX,PixY") print $4}' ../spec_ana/log/Center.log` # y center [pixel]

set bindir=${PWD}/bin
set imgdir="../image_ana"

#setup prefix
set prefixlog="../spec_ana/log/prefix.log"
set mos1prefix=`head -1 $prefixlog | gawk '{print $1}'`
set mos2prefix=`head -2 $prefixlog | tail -1 |  gawk '{print $1}'`
set pnprefix=`tail -1 $prefixlog | gawk '{print $1}'`
#------------------------------------------------------------------------------

set img_mos1=${imgdir}/rate-400-2300-mos1-sp.fits.mask
set sig_mos1=${imgdir}/sigma-400-2300-mos1-sp.fits.mask
set psf_mos1=${imgdir}/mos1_psf_2bin.fits
set Afakeit_mos1=../result/Lambda_${mos1prefix}.dat
set con_mos1=`tail -1 ../result/constant.dat | gawk '{print $1}'`

set img_mos2=${imgdir}/rate-400-2300-mos2-sp.fits.mask
set sig_mos2=${imgdir}/sigma-400-2300-mos2-sp.fits.mask
set psf_mos2=${imgdir}/mos2_psf_2bin.fits
set Afakeit_mos2=../result/Lambda_${mos2prefix}.dat
set con_mos2=1.00000000

set img_pn=${imgdir}/rate-400-2300-pn-sp.fits.mask
set sig_pn=${imgdir}/sigma-400-2300-pn-sp.fits.mask
set psf_pn=${imgdir}/pn_psf_2bin.fits
set Afakeit_pn=../result/Lambda_${pnprefix}.dat
set con_pn=`tail -1 ../result/constant.dat | gawk '{print $3}'`
set pixtoas_PSFmap=2 # pixsel size [arcsec] of PSFmap
 
set xcpsf=$xc  # x center for making PSFmap
set ycpsf=$yc  # x center for making PSFmap
 
set temp=../result/temp.dat
set ne=nresult/ne.dat

#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
echo " # surface brightness fitting"

echo $img_mos1 $sig_mos1 $psf_mos1 $Afakeit_mos1 $con_mos1 $img_mos2 $sig_mos2 $psf_mos2 $Afakeit_mos2 $con_mos2 $img_pn $sig_pn $psf_pn $Afakeit_pn $con_pn $pixtoas_PSFmap $xc $yc $xcpsf $ycpsf $rinn $rout $Nbin
${bindir}/fit_sx_vikhlinin_mos1+mos2+pn $img_mos1 $sig_mos1 $psf_mos1 $Afakeit_mos1 $con_mos1 $img_mos2 $sig_mos2 $psf_mos2 $Afakeit_mos2 $con_mos2 $img_pn $sig_pn $psf_pn $Afakeit_pn $con_pn $pixtoas_PSFmap $xc $yc $xcpsf $ycpsf $rinn $rout $Nbin

#beta model
#${bindir}/fit_sx_beta_mos1+mos2+pn $img_mos1 $sig_mos1 $psf_mos1 $Afakeit_mos1 $con_mos1 $img_mos2 $sig_mos2 $psf_mos2 $Afakeit_mos2 $con_mos2 $img_pn $sig_pn $psf_pn $Afakeit_pn $con_pn $pixtoas_PSFmap $xc $yc $xcpsf $ycpsf $rinn $rout $Nbin

#double beta model
#${bindir}/fit_sx_2beta_mos1+mos2+pn $img_mos1 $sig_mos1 $psf_mos1 $Afakeit_mos1 $con_mos1 $img_mos2 $sig_mos2 $psf_mos2 $Afakeit_mos2 $con_mos2 $img_pn $sig_pn $psf_pn $Afakeit_pn $con_pn $pixtoas_PSFmap $xc $yc $xcpsf $ycpsf $rinn $rout $Nbin

echo " # convert to ne"
python ${bindir}/cal_ne0.py $z $img_mos1

echo " # Temeprature fitting"
#${bindir}/fit_Tx_ew $temp $ne $acool $at $bt
${bindir}/fit_Tx_sl_1styear $temp $ne $acool $at $bt

echo " # M_HE calculation"
python ${bindir}/cal_MHE.py $name $z

#----------------------------
#save in result at the parent dir
cp -r MHEresult ../result/
cp -r Sxresult ../result/
cp -r Tresult ../result/
cp -r nresult ../result/
