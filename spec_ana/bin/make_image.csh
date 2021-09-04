#!/bin/csh -f

if ( ! -d ../image_ana ) then
     mkdir ../image_ana
endif

if ($#argv != 1) then
    echo "#----------------usage----------------#"
    echo "csh make_image.csh binning"
    echo ""
    echo "binning=1 or 2 or 3"
    echo ""
    echo "Example :"
    echo "bin/make_image.csh 1"
    echo "#-------------------------------------#"

exit

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
set maskswitch=1
source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`
set binning=$argv[1]

echo "Binning : $binning"

echo "Running mos-spectra/pn-spectra..."
cat log/06mos-filter.log |grep "Limit"|head -6|egrep "\*\*\*\*" > log/remove${mos1input}_06mos-filter.log
cat log/06mos-filter.log |grep "Limit"|tail -6|egrep "\*\*\*\*" > log/remove${mos2input}_06mos-filter.log

set mos1ccd1 = 1 ;set mos1ccd2 = 1 ;set mos1ccd3 = 1 ;set mos1ccd4 = 1 ;set mos1ccd5 = 1 ;set mos1ccd6 = 1 ;set mos1ccd7 = 1;
set mos2ccd1 = 1 ;set mos2ccd2 = 1 ;set mos2ccd3 = 1 ;set mos2ccd4 = 1 ;set mos2ccd5 = 1 ;set mos2ccd6 = 1 ;set mos2ccd7 = 1;
foreach num (`cat log/remove${mos1input}_06mos-filter.log|awk '{print $2}'|awk -F'[=]' '{print $2}'`)
set mos1ccd${num} = 0
end

foreach num (`cat log/remove${mos2input}_06mos-filter.log|awk '{print $2}'|awk -F'[=]' '{print $2}'`)
set mos2ccd${num} = 0
end

mos-spectra prefix=$mos1input caldb=$SAS_ESAS_CALDB region=nonfile mask=$maskswitch elow=$LOWENE ehigh=$HIGHENE ccd1=$mos1ccd1 ccd2=$mos1ccd2 ccd3=$mos1ccd3 ccd4=$mos1ccd4 ccd5=$mos1ccd5 ccd6=$mos1ccd6 ccd7=$mos1ccd7 >& log/08mos-spectra.log
mos-spectra prefix=$mos2input caldb=$SAS_ESAS_CALDB region=nonfile mask=$maskswitch elow=$LOWENE ehigh=$HIGHENE ccd1=$mos2ccd1 ccd2=$mos2ccd2 ccd3=$mos2ccd3 ccd4=$mos2ccd4 ccd5=$mos2ccd5 ccd6=$mos2ccd6 ccd7=$mos2ccd7 >& log/09mos-spectra.log
pn-spectra prefix=$pninput caldb=$SAS_ESAS_CALDB region=nonfile mask=$maskswitch elow=$LOWENE ehigh=$HIGHENE quad1=1 quad2=1 quad3=1 quad4=1 >& log/10pn-spectra.log

echo "Running mos_back/pn_back..."
mos_back prefix=$mos1input caldb=$SAS_ESAS_CALDB diag=0 elow=$LOWENE ehigh=$HIGHENE ccd1=$mos1ccd1 ccd2=$mos1ccd2 ccd3=$mos1ccd3 ccd4=$mos1ccd4 ccd5=$mos1ccd5 ccd6=$mos1ccd6 ccd7=$mos1ccd7 >& log/11mos_back.log
mos_back prefix=$mos2input caldb=$SAS_ESAS_CALDB diag=0 elow=$LOWENE ehigh=$HIGHENE ccd1=$mos2ccd1 ccd2=$mos2ccd2 ccd3=$mos2ccd3 ccd4=$mos2ccd4 ccd5=$mos2ccd5 ccd6=$mos2ccd6 ccd7=$mos2ccd7 >& log/12mos_back.log
pn_back prefix=$pninput caldb=$SAS_ESAS_CALDB diag=0 elow=$LOWENE ehigh=$HIGHENE quad1=1 quad2=1 quad3=1 quad4=1 >& log/13pn_back.log

echo "Running rot-im-det-sky..."
rot-im-det-sky prefix=$mos1input mask=$maskswitch elow=$LOWENE ehigh=$HIGHENE mode=1 >& log/14rot-im-det-sky.log
rot-im-det-sky prefix=$mos2input mask=$maskswitch elow=$LOWENE ehigh=$HIGHENE mode=1 >& log/15rot-im-det-sky.log
rot-im-det-sky prefix=$pninput mask=$maskswitch elow=$LOWENE ehigh=$HIGHENE mode=1 >& log/16rot-im-det-sky.log

echo "Running comb mos1+mos2+pn..."
#MOS1&MOS2&PN
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=$maskswitch prefixlist="$mos1input $mos2input $pninput" >& log/17comb.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=0 withswcxcontrol=0 clobber=1 >& log/18adapt.log
bin_image thresholdmasking=0.02 detector=0 binning=$binning elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& log/19bin_image.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-all.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-all.fits
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-all.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-all.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-all.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-all.fits

#MOS1&MOS2
echo "Running comb mos1+mos2..."
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=$maskswitch prefixlist="$mos1input $mos2input" >& log/17comb-mos.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=0 withswcxcontrol=0 prefix="$mos1input $mos2input" clobber=1 >& log/18adapt-mos.log
bin_image thresholdmasking=0.02 detector=1 binning=$binning prefix="$mos1input $mos2input" elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& log/19bin_image-mos.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos.fits
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-mos.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos.fits

#MOS1
echo "Running comb mos1..."
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=$maskswitch prefixlist="$mos1input" >& log/17comb-mos1.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=0 withswcxcontrol=0 prefix="$mos1input" clobber=1 >& log/18adapt-mos1.log
bin_image thresholdmasking=0.02 detector=1 binning=$binning prefix="$mos1input" elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& log/19bin_image-mos1.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos1.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos1.fits
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-mos1.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos1.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos1.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos1.fits
#MOS2
echo "Running comb mos2..."
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=$maskswitch prefixlist="$mos2input" >& log/17comb-mos2.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=0 withswcxcontrol=0 prefix="$mos2input" clobber=1 >& log/18adapt-mos2.log
bin_image thresholdmasking=0.02 detector=1 binning=$binning prefix="$mos2input" elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& log/19bin_image-mos2.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos2.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos2.fits
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-mos2.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos2.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos2.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos2.fits
#PN
echo "Running comb pn..."
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=0 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=$maskswitch prefixlist=$pninput >& log/17comb-pn.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=0 withswcxcontrol=0 prefix=$pninput clobber=1 >& log/18adapt-pn.log
bin_image thresholdmasking=0.02 detector=2 binning=$binning prefix=$pninput elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=no withswcxcontrol=no >& log/19bin_image-pn.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-pn.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-pn.fits
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-pn.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-pn.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-pn.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-pn.fits


mv ${mos1prefix}-obj.pi ${mos1prefix}-obj-full.pi
mv ${mos1prefix}.rmf ${mos1prefix}-full.rmf
mv ${mos1prefix}.arf ${mos1prefix}-full.arf
mv ${mos1prefix}-back.pi ${mos1prefix}-back-full.pi
mv ${mos1prefix}-obj-im-sp-det.fits ${mos1prefix}-sp-full.fits
mv ${mos2prefix}-obj.pi ${mos2prefix}-obj-full.pi
mv ${mos2prefix}.rmf ${mos2prefix}-full.rmf
mv ${mos2prefix}.arf ${mos2prefix}-full.arf
mv ${mos2prefix}-back.pi ${mos2prefix}-back-full.pi
mv ${mos2prefix}-obj-im-sp-det.fits ${mos2prefix}-sp-full.fits
mv ${pnprefix}-obj.pi ${pnprefix}-obj-full.pi
mv ${pnprefix}.rmf ${pnprefix}-full.rmf
mv ${pnprefix}.arf ${pnprefix}-full.arf
mv ${pnprefix}-back.pi ${pnprefix}-back-full.pi
mv ${pnprefix}-obj-im-sp-det.fits ${pnprefix}-sp-full.fits


if( -f ${mos1prefix}-obj-full-grp.pi ) then
rm -f ${mos1prefix}-obj-full-grp.pi
endif
if( -f ${mos2prefix}-obj-full-grp.pi ) then
rm -f ${mos2prefix}-obj-full-grp.pi
endif
if( -f ${pnprefix}-obj-full-grp.pi ) then
rm -f ${pnprefix}-obj-full-grp.pi
endif


grppha ${mos1prefix}-obj-full.pi ${mos1prefix}-obj-full-grp.pi <<EOF
chkey BACKFILE ${mos1prefix}-back-full.pi
chkey RESPFILE ${mos1prefix}-full.rmf
chkey ANCRFILE ${mos1prefix}-full.arf
group min 100 
exit
EOF
grppha ${mos2prefix}-obj-full.pi ${mos2prefix}-obj-full-grp.pi <<EOF
chkey BACKFILE ${mos2prefix}-back-full.pi
chkey RESPFILE ${mos2prefix}-full.rmf
chkey ANCRFILE ${mos2prefix}-full.arf
group min 100 
exit
EOF
grppha ${pnprefix}-obj-full.pi ${pnprefix}-obj-full-grp.pi <<EOF
chkey BACKFILE ${pnprefix}-back-full.pi
chkey RESPFILE ${pnprefix}-full.rmf
chkey ANCRFILE ${pnprefix}-full.arf
group min 100 
exit
EOF


rm rate-$LOWENE-$HIGHENE-mos1.fits.mask rate-$LOWENE-$HIGHENE-mos2.fits.mask rate-$LOWENE-$HIGHENE-pn.fits.mask
rm sigma-$LOWENE-$HIGHENE-mos1.fits.mask sigma-$LOWENE-$HIGHENE-mos2.fits.mask sigma-$LOWENE-$HIGHENE-pn.fits.mask
farith rate-$LOWENE-$HIGHENE-mos1.fits ${mos1prefix}-cheese.fits rate-$LOWENE-$HIGHENE-mos1.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-mos1.fits ${mos1prefix}-cheese.fits sigma-$LOWENE-$HIGHENE-mos1.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-mos2.fits ${mos2prefix}-cheese.fits rate-$LOWENE-$HIGHENE-mos2.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-mos2.fits ${mos2prefix}-cheese.fits sigma-$LOWENE-$HIGHENE-mos2.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-pn.fits ${pnprefix}-cheese.fits rate-$LOWENE-$HIGHENE-pn.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-pn.fits ${pnprefix}-cheese.fits sigma-$LOWENE-$HIGHENE-pn.fits.mask MUL

#ファイルをイメージ解析用ディレクトリにコピーする。
cp * ../image_ana/



#--------------------------------------------------------------------------
# Check whether output files exist or not
#-------------------------------------------------------------------------
set errlog="log/make_image_errfiles.log"
echo "# These files were not created" > $errlog 


echo "Checking whether outputs exist or not...."
foreach prefix ( "all" "mos" "mos1" "mos2" "pn" )

set rate=rate-$LOWENE-$HIGHENE-${prefix}.fits
set sigma=sigma-$LOWENE-$HIGHENE-${prefix}.fits
set adapt=adapt-$LOWENE-$HIGHENE-${prefix}.fits
set comb=comb-obj-im-$LOWENE-$HIGHENE-${prefix}.fits
set combback=comb-back-im-sky-$LOWENE-$HIGHENE-${prefix}.fits
set combexp=comb-exp-im-$LOWENE-$HIGHENE-${prefix}.fits

if (! -f $rate ) then
echo $rate >> $errlog
endif
if (! -f $sigma ) then
echo $sigma >> $errlog
endif
if (! -f $adapt && $prefix != "mos") then
echo $adapt >> $errlog
endif
if (! -f $comb) then
echo $comb  >> $errlog
endif
if (! -f $combback ) then
echo $combback >> $errlog
endif

end


if ( `more $errlog | wc -l` > 1 ) then;
echo "======================================================="
echo "WARNING : Some files were not producted"
echo "see $errlog"
echo "======================================================="
more $errlog
echo "#------------------------------------------------------"
else 
echo "======================================================="
echo "Success  : all files are completed"
echo "Next run : "
echo "bin/SearchCentroid.py"
echo "======================================================="
endif

