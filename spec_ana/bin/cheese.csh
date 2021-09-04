#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`

if ($#argv != 3) then
    echo "#----------------usage----------------#"
    echo "csh cheese.csh scale rate dist"
    echo ""
    echo "scale=0.25"
    echo "rate=0.001"
    echo "dist=20"
    echo ""
    echo "Example :"
    echo "bin/cheese.csh 0.15 0.001 20"
    echo "#-------------------------------------#"

exit

endif
 
set scale = $argv[1]
set rate = $argv[2]
set dist = $argv[3]

#--------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------
echo "    Running cheese...."
cheese prefixm="$mos1input $mos2input" prefixp=$pninput scale=$scale rate=$rate dist=$dist clobber=1 elow=$LOWENE ehigh=$HIGHENE >& log/07cheese.log

echo "    Running ModifyRegionR.py ...."
bin/ModifyRegionR.py 

cp ${mos1prefix}-bkg_region-det.fits ${mos1prefix}-bkg_region-det.fits.log
cp ${mos1prefix}-bkg_region-sky.fits ${mos1prefix}-bkg_region-sky.fits.log
cp ${mos2prefix}-bkg_region-det.fits ${mos2prefix}-bkg_region-det.fits.log
cp ${mos2prefix}-bkg_region-sky.fits ${mos2prefix}-bkg_region-sky.fits.log
cp ${pnprefix}-bkg_region-det.fits ${pnprefix}-bkg_region-det.fits.log
cp ${pnprefix}-bkg_region-sky.fits ${pnprefix}-bkg_region-sky.fits.log

cp ${mos1prefix}-cheese.fits ${mos1prefix}-cheese.fits.log
cp ${mos2prefix}-cheese.fits ${mos2prefix}-cheese.fits.log
cp ${pnprefix}-cheese.fits ${pnprefix}-cheese.fits.log

echo "   Running makeRegionFile.py.... "
echo "     output : ${mos1prefix}-bkg_region-sky.reg"
#making mos1S001-bkg_region-sky.reg
bin/makeRegionFile.py

if ( -z ${mos1prefix}-bkg_region-sky.reg ) then;
echo "ERROR : Failed to make region file"
exit  
endif

echo "===================================================="
echo " NEXT TO DO "
echo "===================================================="
echo "1: Remove fake point sources"
echo "ds9 ${mos1prefix}-obj-image-sky.fits &"
echo "    open ${mos1prefix}-bkg_region-sky.reg"
echo " "
echo "2:bin/make_mask.csh savedregfile" 
#>ds9 mos1S001-obj-image-sky.fits &
#mos1S001-bkg_region-sky.reg
#で目でみてフェイクリージョンを省く

