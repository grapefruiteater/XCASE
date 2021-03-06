#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`
#======================================================================================================
if ($#argv != 3) then
    echo "#----------------usage----------------#"
    echo "remake_spec.csh det rinn rout"
    echo ""
    echo "When you failed to make some spec (pi/back) files in make_spec.csh,"
    echo "you would better run this script"
    echo ""
    echo "det  : mos1/mos2/pn"
    echo "rinn : rinn"
    echo "rinn : rout"
    echo ""
    echo "Example :"
    echo "bin/remake_spec.csh pn 240 360"
    echo "#-------------------------------------#"
exit
endif

set prefix0=$1
set rinn=$2
set rout=$3
set regionname=${rinn}-${rout}

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

#----------------------------------------------------------------
set check=`gawk 'BEGIN{check=0}{if($1=='${rinn}'){check=1;}}END{print check}' $RegionFile`
if ( $check == 0 ) then
echo "ERROR : enter the proper rinn"
more $RegionFile
exit
endif
set check=`gawk 'BEGIN{check=0}{if($1=='${rout}'){check=1;}}END{print check}' $RegionFile`
if ( $check == 0 ) then
echo "ERROR : enter the proper rout"
more $RegionFile
exit
endif
#--------------------------------
if ( $prefix0 == "mos1") then
set prefix=$mos1prefix

echo $regionname" / mos1"
mos-spectra prefix=$mos1input caldb=$SAS_ESAS_CALDB region=${mos1prefix}_${regionname}secR.reg mask=1 elow=0 ehigh=0 ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/20mos1-spectra-${regionname}.log
mos_back prefix=$mos1input caldb=$SAS_ESAS_CALDB diag=2 elow=0 ehigh=0 ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/21mos1_back-${regionname}.log

mv ${mos1prefix}-obj.pi ${mos1prefix}-obj-${regionname}.pi
mv ${mos1prefix}-back.pi ${mos1prefix}-back-${regionname}.pi
mv ${mos1prefix}.rmf ${mos1prefix}-${regionname}.rmf
mv ${mos1prefix}.arf ${mos1prefix}-${regionname}.arf
mv ${mos1prefix}-obj-im-sp-det.fits ${mos1prefix}-sp-${regionname}.fits


else if ( $prefix0 == "mos2") then
set prefix=$mos2prefix

echo $regionname" / mos2"
mos-spectra prefix=$mos2input caldb=$SAS_ESAS_CALDB region=${mos2prefix}_${regionname}secR.reg mask=1 elow=0 ehigh=0 ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/20mos2-spectra-${regionname}.log
mos_back prefix=$mos2input caldb=$SAS_ESAS_CALDB diag=2 elow=0 ehigh=0 ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 >& log/21mos2_back-${regionname}.log

mv ${mos2prefix}-obj.pi ${mos2prefix}-obj-${regionname}.pi
mv ${mos2prefix}-back.pi ${mos2prefix}-back-${regionname}.pi
mv ${mos2prefix}.rmf ${mos2prefix}-${regionname}.rmf
mv ${mos2prefix}.arf ${mos2prefix}-${regionname}.arf
mv ${mos2prefix}-obj-im-sp-det.fits ${mos2prefix}-sp-${regionname}.fits

else if ( $prefix0 == "pn") then
set prefix=$pnprefix

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


else
echo "Please enter mos1/mos2/pn"
exit 
endif






#--------------------------------------------------------------------------
# Check whether output files exist or not
#-------------------------------------------------------------------------
set errlog="log/remake_spec_errfiles.log"
echo "# These files were not created" > $errlog 

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

if ( `more $errlog | wc -l` > 1 ) then;
echo "======================================================="
echo "WARNING : Some files were not producted"
echo "see $errlog"
echo "======================================================="
more $errlog
echo "-------------------------------------------------------"
else 
echo "======================================================="
echo "Success  : all files are completed"
echo "Next run : "
echo "bin/run_grppha.csh Ngrp"
echo "======================================================="
endif





