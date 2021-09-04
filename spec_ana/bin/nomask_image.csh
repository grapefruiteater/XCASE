#!/bin/csh -f

if ( ! -d ../image_ana ) then
     mkdir ../image_ana
endif
#---------------------------------------------------------
#setup prefix
set prefixlog="log/prefix.log"
set mos1prefix=`head -1 $prefixlog | gawk '{print $1}'`
set mos2prefix=`head -2 $prefixlog | tail -1 |  gawk '{print $1}'`
set pnprefix=`tail -1 $prefixlog | gawk '{print $1}'`

set mos1input=`head -1 $prefixlog | gawk '{print $2}'`
set mos2input=`head -2 $prefixlog | tail -1 |  gawk '{print $2}'`
set pninput=`tail -1 $prefixlog | gawk '{print $2}'`
set mosinput=`echo "'"$mos1input" "$mos2input"'"`

#setup energy range
set eRangeFile="SetUpFile/eRange.dat"
set LOWENE=`head -1 $eRangeFile | gawk '{print $1}'`
set HIGHENE=`head -2 $eRangeFile | gawk 'NR==2 {print $1}'`
#----------------------------------------------------------------

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`

echo "Running mos-spectra/pn-spectra..."
mos-spectra prefix=$mos1input caldb=$SAS_ESAS_CALDB region=nonfile mask=0 elow=$LOWENE ehigh=$HIGHENE ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/08mos-spectra_nomask.log
mos-spectra prefix=$mos2input caldb=$SAS_ESAS_CALDB region=nonfile mask=0 elow=$LOWENE ehigh=$HIGHENE ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/09mos-spectra_nomask.log
pn-spectra prefix=$pninput caldb=$SAS_ESAS_CALDB region=nonfile mask=0 elow=$LOWENE ehigh=$HIGHENE quad1=1 quad2=1 quad3=1 quad4=1 >& log/10pn-spectra_nomask.log

echo "Running mos_back/pn_back..."
mos_back prefix=$mos1input caldb=$SAS_ESAS_CALDB diag=0 elow=$LOWENE ehigh=$HIGHENE ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/11mos_back_nomask.log
mos_back prefix=$mos2input caldb=$SAS_ESAS_CALDB diag=0 elow=$LOWENE ehigh=$HIGHENE ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/12mos_back_nomask.log
pn_back prefix=$pninput caldb=$SAS_ESAS_CALDB diag=0 elow=$LOWENE ehigh=$HIGHENE quad1=1 quad2=1 quad3=1 quad4=1 >& log/13pn_back_nomask.log

echo "Running rot-im-det-sky..."
rot-im-det-sky prefix=$mos1input mask=0 elow=$LOWENE ehigh=$HIGHENE mode=1 >& log/14rot-im-det-sky_nomask.log
rot-im-det-sky prefix=$mos2input mask=0 elow=$LOWENE ehigh=$HIGHENE mode=1 >& log/15rot-im-det-sky_nomask.log
rot-im-det-sky prefix=$pninput mask=0 elow=$LOWENE ehigh=$HIGHENE mode=1 >& log/16rot-im-det-sky_nomask.log

echo "Running comb mos1+mos2+pn..."
#MOS1&MOS2&PN
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=0 prefixlist="$mos1input $mos2input $pninput" >& log/17comb_nomask.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=0 withswcxcontrol=0 clobber=1 >& log/18adapt_nomask.log
bin_image thresholdmasking=0.02 detector=0 binning=1 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& log/19bin_image_nomask.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-all_nomask.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-all_nomask.fits
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-all_nomask.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-all_nomask.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-all_nomask.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-all_nomask.fits

#MOS1&MOS2
echo "Running comb mos1+mos2..."
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=0 prefixlist="$mos1input $mos2input" >& log/17comb-mos_nomask.log
bin_image thresholdmasking=0.02 detector=1 binning=1 prefix="$mos1input $mos2input" elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& log/19bin_image-mos_nomask.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos_nomask.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos_nomask.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos_nomask.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos_nomask.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos_nomask.fits

#MOS1
echo "Running comb mos1..."
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=0 prefixlist="$mos1input" >& log/17comb-mos1_nomask.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=0 withswcxcontrol=0 prefix="$mos1input" clobber=1 >& log/18adapt-mos1_nomask.log
bin_image thresholdmasking=0.02 detector=1 binning=1 prefix="$mos1input" elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& log/19bin_image-mos1_nomask.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos1_nomask.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos1_nomask.fits
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-mos1_nomask.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos1_nomask.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos1_nomask.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos1_nomask.fits
#MOS2
echo "Running comb mos2..."
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=0 prefixlist="$mos2input" >& log/17comb-mos2_nomask.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=0 withswcxcontrol=0 prefix="$mos2input" clobber=1 >& log/18adapt-mos2_nomask.log
bin_image thresholdmasking=0.02 detector=1 binning=1 prefix="$mos2input" elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& log/19bin_image-mos2_nomask.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos2_nomask.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos2_nomask.fits
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-mos2_nomask.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos2_nomask.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos2_nomask.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos2_nomask.fits
#PN
echo "Running comb pn..."
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=0 prefixlist=$pninput >& log/17comb-pn_nomask.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=0 withswcxcontrol=0 prefix=$pninput clobber=1 >& log/18adapt-pn_nomask.log
bin_image thresholdmasking=0.02 detector=2 binning=1 prefix=$pninput elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& log/19bin_image-pn_nomask.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-pn_nomask.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-pn_nomask.fits
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-pn_nomask.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-pn_nomask.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-pn_nomask.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-pn_nomask.fits





