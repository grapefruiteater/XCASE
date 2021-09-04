#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`
#======================================================================================================
if ($#argv != 1) then
    echo "#----------------usage----------------#"
    echo "make_spec.csh grpbin"
    echo ""
    echo "grpbin : bin for grppha"
    echo ""
    echo "Example :"
    echo "bin/make_spec.csh 50"
    echo "#-------------------------------------#"
exit
endif

set grpbin=$1

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

#setup region file
set RegionFile="SetUpFile/Region.dat"

#setup energy range
set eRangeFile="SetUpFile/eRange.dat"
set LOWENE=`head -1 $eRangeFile | gawk '{print $1}'`
set HIGHENE=`head -2 $eRangeFile | gawk 'NR==2 {print $1}'`

#setup group bin 
set GrpBinFile="SetUpFile/GrpBin.dat"
echo $grpbin > $GrpBinFile
#----------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
if ( ! -f $RegionFile ) then
     echo ""
     echo "ERROR : NO "$RegionFile
     echo ""
     exit
endif
bin/setUpSpec.py
#----------------------------------------------------------------------------------------------------------------------------------------
if ( ! -f spspec-mos1.txt ) then
     echo ""
     echo "ERROR : NO spspec-mos1.txt"
     echo ""
     exit
endif

if ( ! -f spspec-mos2.txt ) then
     echo ""
     echo "ERROR : NO spspec-mos2.txt"
     echo ""
     exit
endif

if ( ! -f spspec-pn.txt ) then
     echo ""
     echo "ERROR : NO spspec-pn.txt"
     echo ""
     exit
endif
if ( ! -d ../data ) then
    mkdir ../data
endif

#----------------------------------------------------------------------------------------------------------------------------------------
#銀河団中心から同心円円環領域のスペクトルとレスポンスを作成

rm -fv *-obj-im.fits *obj-im-sp-det.fits *exp-im.fits *mask-im.fits *obj.pi *obj-oot.pi 

echo "------------------------------"
echo "# Radii for regions"
set RediusList = ()
set COUNT = `more $RegionFile | wc -l`
@ i = 1
while($i <= $COUNT)
    set Radius = "`cat $RegionFile | head -$i |tail -1`"
    echo "$Radius"
    set RediusList = ($RediusList $Radius)
    @ i++
end
echo "------------------------------"

set i=0
@ max = ${#RediusList} - 1

while($i<$max)
@ i = $i + 1
set inner = $RediusList[$i]
@ j = $i + 1
set outer = $RediusList[$j]
set regionname = $inner-$outer
echo $regionname" / mos1"
mos-spectra prefix=$mos1input caldb=$SAS_ESAS_CALDB region=${mos1prefix}_${regionname}secR.reg mask=1 elow=0 ehigh=0 ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/20mos1-spectra-${regionname}.log
mos_back prefix=$mos1input caldb=$SAS_ESAS_CALDB diag=2 elow=0 ehigh=0 ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/21mos1_back-${regionname}.log

mv ${mos1prefix}-obj.pi ${mos1prefix}-obj-${regionname}.pi
mv ${mos1prefix}-back.pi ${mos1prefix}-back-${regionname}.pi
mv ${mos1prefix}.rmf ${mos1prefix}-${regionname}.rmf
mv ${mos1prefix}.arf ${mos1prefix}-${regionname}.arf
mv ${mos1prefix}-obj-im-sp-det.fits ${mos1prefix}-sp-${regionname}.fits
end

set i=0
while($i<$max)
@ i = $i + 1 
set inner = $RediusList[$i]
@ j = $i + 1
set outer = $RediusList[$j]
set regionname = $inner-$outer
echo $regionname" / mos2"
mos-spectra prefix=$mos2input caldb=$SAS_ESAS_CALDB region=${mos2prefix}_${regionname}secR.reg mask=1 elow=0 ehigh=0 ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/20mos2-spectra-${regionname}.log
mos_back prefix=$mos2input caldb=$SAS_ESAS_CALDB diag=2 elow=0 ehigh=0 ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/21mos2_back-${regionname}.log

mv ${mos2prefix}-obj.pi ${mos2prefix}-obj-${regionname}.pi
mv ${mos2prefix}-back.pi ${mos2prefix}-back-${regionname}.pi
mv ${mos2prefix}.rmf ${mos2prefix}-${regionname}.rmf
mv ${mos2prefix}.arf ${mos2prefix}-${regionname}.arf
mv ${mos2prefix}-obj-im-sp-det.fits ${mos2prefix}-sp-${regionname}.fits
end

set i=0
while($i<$max)
@ i = $i + 1 
set inner = $RediusList[$i]
@ j = $i + 1
set outer = $RediusList[$j]
set regionname = $inner-$outer
echo $regionname" / pn"
pn-spectra prefix=$pninput caldb=$SAS_ESAS_CALDB region=${pnprefix}_${regionname}secR.reg mask=1 elow=0 ehigh=0 quad1=1 quad2=1 quad3=1 quad4=1 >& log/20pn-spectra-${regionname}.log
pn_back prefix=$pninput caldb=$SAS_ESAS_CALDB diag=0 elow=0 ehigh=0 quad1=1 quad2=1 quad3=1 quad4=1 >& log/21pn_back-${regionname}.log 

mv ${pnprefix}-obj-os.pi ${pnprefix}-obj-os-${regionname}.pi
mv ${pnprefix}-obj.pi ${pnprefix}-obj-${regionname}.pi
mv ${pnprefix}-obj-oot.pi ${pnprefix}-obj-oot-${regionname}.pi
mv ${pnprefix}-back.pi ${pnprefix}-back-${regionname}.pi
mv ${pnprefix}.rmf ${pnprefix}-${regionname}.rmf
mv ${pnprefix}.arf ${pnprefix}-${regionname}.arf
mv ${pnprefix}-obj-im-sp-det.fits ${pnprefix}-sp-${regionname}.fits
end

#各円環のSolid AngleとSoft Proton normalization scale factor を計算する。
echo "Running proton_scale...."
proton_scale caldb=$SAS_ESAS_CALDB mode=2 detector=1 spfile=spspec-mos1.txt >& log/mos1protonscale.log
proton_scale caldb=$SAS_ESAS_CALDB mode=2 detector=2 spfile=spspec-mos2.txt >& log/mos2protonscale.log
proton_scale caldb=$SAS_ESAS_CALDB mode=2 detector=3 spfile=spspec-pn.txt >& log/pnprotonscale.log


bin/run_grppha.csh $grpbin

#--------------------------------------------------------------------------
# Check whether output files exist or not
#-------------------------------------------------------------------------
set errlog="log/make_spec_errfiles.log"
echo "# These files were not created" > $errlog 


echo "Checking whether output pi/rmf/arf exist or not...."
foreach prefix ( $mos1prefix $mos2prefix $pnprefix )
set i=0
while($i<$max)
@ i = $i + 1
set inner = $RediusList[$i]
@ j = $i + 1
set outer = $RediusList[$j]
set regionname = $inner-$outer

if (! -f ${prefix}-obj-${regionname}.pi) then
echo ${prefix}-obj-${regionname}.pi >> $errlog
endif
if (! -f ${prefix}-back-${regionname}.pi) then
echo ${prefix}-back-${regionname}.pi >> $errlog
endif
if (! -f ${prefix}-${regionname}.rmf) then
echo ${prefix}-${regionname}.rmf >> $errlog
endif
if (! -f ${prefix}-${regionname}.arf) then
echo ${prefix}-${regionname}.arf  >> $errlog
endif
if (! -f ${prefix}-sp-${regionname}.fits) then
echo ${prefix}-sp-${regionname}.fits >> $errlog
endif

end
end

if ( `more $errlog | wc -l` > 1 ) then;
echo "======================================================="
echo "WARNING : Some files were not producted"
echo "see $errlog"
echo "======================================================="
more $errlog
echo "======================================================="
echo "Next run : "
echo "bin/remake_spec.csh det rinn rout"
echo "======================================================="

echo "-------------------------------------------------------"
else 
echo "======================================================="
echo "Success  : all files are completed"
echo "Next run : "
echo "bin/make_xspecscript.csh z sample.xcm"
echo "======================================================="
endif





