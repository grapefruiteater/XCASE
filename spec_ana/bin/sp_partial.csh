#!/bin/csh -f

if ( $#argv != 1 ) then
cat <<EOF

bin/sp_partial.csh bestfitfile errorfile

   bestfitfile : bestfit output of XSPEC

   output : 
          log/sp_partial-mos1.log  
          log/sp_partial-mos2.log  
          log/sp_partial-pn.log  

Use:  bin/sp_partial.csh save-all.xcm


EOF

exit

endif

set xspecoutput=$1
#set xspecoutputErr=$2

echo "# input xpsec output file" > log/sp_partial-inputfile.log
echo "Best : "$xspecoutput >> log/sp_partial-inputfile.log
#echo "Err  : "$xspecoutputErr >> log/sp_partial-inputfile.log


source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`


set RediusList = ()
set RegionFile = 'SetUpFile/Region.dat'
set COUNT = `wc $RegionFile | awk '{print $1}'`
@ i = 1
while($i <= $COUNT)
    set Radius = "`cat $RegionFile | head -$i |tail -1`"
    #echo "$Radius"
    set RediusList = ($RediusList $Radius)
    @ i++
end
echo "Radii in arcsec"
echo $RediusList
#---------------------------------------------------------
#setup prefix
set prefixlog="log/prefix.log"
set mos1prefix=`head -1 $prefixlog | gawk '{print $1}'`
set mos2prefix=`head -2 $prefixlog | tail -1 |  gawk '{print $1}'`
set pnprefix=`tail -1 $prefixlog | gawk '{print $1}'`

#setup energy range
set eRangeFile="SetUpFile/eRange.dat"
set LOWENE=`head -1 $eRangeFile | gawk '{print $1}'`
set HIGHENE=`head -2 $eRangeFile | gawk 'NR==2 {print $1}'`
#----------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
#Soft proton background image作成のために必要な Scaled Normalization を計算。
#スペクトルフィットから得られるソフトプロトンのpowerlawのnormのベストフィットパラメータを入れる。

set Nmos1=`gawk '{if($1=="model" &&  $2=="2:mos1_1" && $3=="powerlaw") print NR+2}' $xspecoutput`
set Nmos2=`gawk '{if($1=="model" &&  $2=="3:mos2_1" && $3=="powerlaw") print NR+2}' $xspecoutput`
set Npn=`gawk '{if($1=="model" &&  $2=="4:pn_1" && $3=="powerlaw") print NR+2}' $xspecoutput`


set mos1rnorm = `head -$Nmos1 $xspecoutput | tail -1 | gawk '{print $1}'`
set mos2rnorm = `head -$Nmos2 $xspecoutput | tail -1 | gawk '{print $1}'`
set pnrnorm = `head -$Npn $xspecoutput | tail -1 | gawk '{print $1}'`

echo "#--------------------------------------"
echo "# Best-fit normalization for soft proton"
echo "mos1 :"$mos1rnorm
echo "mos2 :"$mos2rnorm
echo "pn   :"$pnrnorm
echo "#--------------------------------------"

#----------------------------------------------------------------------------------
set inner = $RediusList[1]
set outer = $RediusList[2]
set regionname = $inner-$outer


foreach prefix ( $mos1prefix $mos2prefix $pnprefix )

if (! -f ${prefix}-sp-full.fits) then
echo "ERROR : no file of "${prefix}-sp-full.fits 
exit 
endif
if (-z ${prefix}-sp-full.fits) then
echo "ERROR : "${prefix}-sp-full.fits" exists but zero size" 
exit 
endif

if (! -f ${prefix}-obj-full.pi) then
echo "ERROR : no file of "${prefix}-obj-full.pi
exit 
endif
if (-z ${prefix}-obj-full.pi) then
echo "ERROR : "${prefix}-obj-full.pi" exists but zero size" 
exit 
endif

if (! -f ${prefix}-sp-${regionname}.fits) then
echo "ERROR : no file of "${prefix}-sp-${regionname}.fits
exit 
endif
if (-z ${prefix}-sp-${regionname}.fits) then
echo "ERROR : "${prefix}-sp-${regionname}.fits" exists but zero size" 
exit 
endif

if (! -f ${prefix}-obj-${regionname}.pi) then
echo "ERROR : no file of "${prefix}-obj-${regionname}.pi
exit 
endif
if (-z ${prefix}-obj-${regionname}.pi) then
echo "ERROR : "${prefix}-obj-${regionname}.pi" exists but zero size" 
exit 
endif


end


sp_partial caldb=$SAS_ESAS_CALDB detector=1 fullimage=${mos1prefix}-sp-full.fits fullspec=${mos1prefix}-obj-full.pi regionimage=${mos1prefix}-sp-${regionname}.fits regionspec=${mos1prefix}-obj-${regionname}.pi rnorm=$mos1rnorm  >& log/sp_partial-mos1.log
sp_partial caldb=$SAS_ESAS_CALDB detector=2 fullimage=${mos2prefix}-sp-full.fits fullspec=${mos2prefix}-obj-full.pi regionimage=${mos2prefix}-sp-${regionname}.fits regionspec=${mos2prefix}-obj-${regionname}.pi rnorm=$mos2rnorm >& log/sp_partial-mos2.log
sp_partial caldb=$SAS_ESAS_CALDB detector=3 fullimage=${pnprefix}-sp-full.fits fullspec=${pnprefix}-obj-full.pi regionimage=${pnprefix}-sp-${regionname}.fits regionspec=${pnprefix}-obj-${regionname}.pi rnorm=$pnrnorm >& log/sp_partial-pn.log

#----------------------------------------------------------------------------------
# check output
gawk '{if($1=="Scaled" && $2=="Normalization:") print $3}' log/sp_partial-mos1.log > tmp
if (-z tmp ) then
echo "ERROR : failed to calculate the normalization of mos1"
exit
endif
gawk '{if($1=="Scaled" && $2=="Normalization:") print $3}' log/sp_partial-mos2.log > tmp
if (-z tmp ) then
echo "ERROR : failed to calculate the normalization of mos2"
exit
endif
gawk '{if($1=="Scaled" && $2=="Normalization:") print $3}' log/sp_partial-pn.log > tmp
if (-z tmp ) then
echo "ERROR : failed to calculate the normalization of pn"
exit
endif


echo "======================================================="
echo "Success  : all files are completed"
echo "Next run : "
echo "cd ../image_ana"
echo "bin/proton.csh"     
echo "======================================================="

