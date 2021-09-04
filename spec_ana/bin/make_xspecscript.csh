#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`
#======================================================================================================
if ($#argv != 2) then
    echo "#----------------usage----------------#"
    echo "make_xspecscript.csh z sample.xcm"
    echo ""
    echo "z          : cluster redshift"
    echo "sample.xcm : sample of models"
    echo ""
    echo "#1 : auto-check nH"

    echo "output     : spectrum.xcm"
    echo "Example :"
    echo "make_xspecscript.csh 0.1832 sample.xcm"
    echo "#-------------------------------------#"
exit
endif

set z=$1
set sample=$2

set centerlog="log/Center.log"
set ra=`head -1 $centerlog | gawk '{print $3}'`
set dec=`head -1 $centerlog | gawk '{print $4}'`
#echo $ra,$dec
#=======================================================================================================
echo "Checking nH...."
nh 2000 $ra $dec > log/nh.dat

set nH=`tail -1 log/nh.dat | gawk '{print $7/1e22}'`
echo "Weighted nH [10^22 cm^-2] : "$nH

echo "Dowloading ROSAT data..."
set radec=`python bin/convert_radec.py $ra $dec`
echo $radec

bin/rosat.rb $radec

if ( -f rass.pi) then
rm rass.pi
endif

grppha rass_orig.pi rass.pi<<EOF
chkey RESPFILE pspcc.rsp
exit
EOF


echo "Making XspecSript..."
python bin/add-xcm.py $z $nH $sample


if ( -z spectrum.xcm ) then
echo "Failed to make Xspec script" 
exit
endif

if ( -z rass.pi ) then
echo "ERROR : no ROSAT data" 
exit
endif


echo "===================================================="
echo " Next run : @spectrum.xcm in xspec"
echo "===================================================="



