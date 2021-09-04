#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`

#-----------------------------------------------------------------------------------------
#setup prefix
set prefixlog="../spec_ana/log/prefix.log"
set mos1prefix=`head -1 $prefixlog | gawk '{print $1}'`
set mos2prefix=`head -2 $prefixlog | tail -1 |  gawk '{print $1}'`
set pnprefix=`tail -1 $prefixlog | gawk '{print $1}'`

set mos1input=`head -1 $prefixlog | gawk '{print $2}'`
set mos2input=`head -2 $prefixlog | tail -1 |  gawk '{print $2}'`
set pninput=`tail -1 $prefixlog | gawk '{print $2}'`
set mosinput=`echo "'"$mos1input" "$mos2input"'"`

#setup energy range
set eRangeFile="../spec_ana/SetUpFile/eRange.dat"
set LOWENE=`head -1 $eRangeFile | gawk '{print $1}'`
set HIGHENE=`head -2 $eRangeFile | gawk 'NR==2 {print $1}'`
#-----------------------------------------------------------------------------------------
###---make image---

#MOS1&MOS2&PN
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="$mos1input $mos2input $pninput" >& log/21comb.log
bin_image thresholdmasking=0.02 detector=0 prefix="$mos1input $mos2input $pninput" binning=1 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=no >& log/22bin_image.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=0 >& log/23adapt.log
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-all-sp.fits
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-all-sp.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-all-sp.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-all-sp.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-all-sp.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-all-sp.fits
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-all-sp.fits
#MOS1&MOS2
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="$mos1input $mos2input" >& log/21comb_mos.log
bin_image thresholdmasking=0.02 detector=0 prefix="$mos1input $mos2input" binning=1 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=no >& log/23bin_image_mos.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos-sp.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos-sp.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos-sp.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos-sp.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos-sp.fits
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos-sp.fits
#MOS1
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="$mos1input" >& log/21comb_mos1.log
bin_image thresholdmasking=0.02 detector=1 prefix="$mos1input" binning=1 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=no >& log/23bin_image_mos1.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos1-sp.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos1-sp.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos1-sp.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos1-sp.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos1-sp.fits
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos1-sp.fits
#MOS2
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="$mos2input" >& log/21comb_mos2.log
bin_image thresholdmasking=0.02 detector=1 prefix="$mos2input" binning=1 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=no >& log/23bin_image_mos2.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos2-sp.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos2-sp.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos2-sp.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos2-sp.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos2-sp.fits
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos2-sp.fits
#PN
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="$pninput" >& log/21comb_pn.log
bin_image thresholdmasking=0.02 detector=2 binning=1 prefix="$pninput" elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=no >& log/23bin_image_pn.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-pn-sp.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-pn-sp.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-pn-sp.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-pn-sp.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-pn-sp.fits
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-pn-sp.fits

#Remask

foreach prefix ( "mos1" "mos2" "pn" )
if ( -f rate-${LOWENE}-${HIGHENE}-${prefix}-sp.fits.mask) then
rm rate-${LOWENE}-${HIGHENE}-${prefix}-sp.fits.mask
endif
if ( -f sigma-${LOWENE}-${HIGHENE}-${prefix}-sp.fits.mask) then
rm sigma-${LOWENE}-${HIGHENE}-${prefix}-sp.fits.mask
endif
end

rm rate-$LOWENE-$HIGHENE-mos1-sp.fits.mask rate-$LOWENE-$HIGHENE-mos2-sp.fits.mask rate-$LOWENE-$HIGHENE-pn-sp.fits.mask
rm sigma-$LOWENE-$HIGHENE-mos1-sp.fits.mask sigma-$LOWENE-$HIGHENE-mos2-sp.fits.mask sigma-$LOWENE-$HIGHENE-pn-sp.fits.mask
farith rate-$LOWENE-$HIGHENE-mos1-sp.fits ${mos1prefix}-cheese.fits rate-$LOWENE-$HIGHENE-mos1-sp.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-mos1-sp.fits ${mos1prefix}-cheese.fits sigma-$LOWENE-$HIGHENE-mos1-sp.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-mos2-sp.fits ${mos2prefix}-cheese.fits rate-$LOWENE-$HIGHENE-mos2-sp.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-mos2-sp.fits ${mos2prefix}-cheese.fits sigma-$LOWENE-$HIGHENE-mos2-sp.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-pn-sp.fits ${pnprefix}-cheese.fits rate-$LOWENE-$HIGHENE-pn-sp.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-pn-sp.fits ${pnprefix}-cheese.fits sigma-$LOWENE-$HIGHENE-pn-sp.fits.mask MUL


#---------------------------------------------------------------------------------------------
foreach prefix ( "mos1" "mos2" "pn" )

if (! -f rate-${LOWENE}-${HIGHENE}-${prefix}-sp.fits.mask) then
echo "ERROR : no file of "rate-${LOWENE}-${HIGHENE}-${prefix}-sp.fits.mask
exit 
endif
if (! -f sigma-${LOWENE}-${HIGHENE}-${prefix}-sp.fits.mask) then
echo "ERROR : no file of "sigma-${LOWENE}-${HIGHENE}-${prefix}-sp.fits.mask
exit 
endif
end
#---------------------------------------------------------------------------------


#copy file
cp rate-**.fits.mask sigma-**.fits.mask ../data/


echo "======================================================="
echo "Success  : all files are completed"
echo "Next run : "
echo "bin/psfgen.csh"
echo "======================================================="

