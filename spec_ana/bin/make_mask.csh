#!/bin/csh -f

if ( $#argv != 1 ) then
cat <<EOF

bin/make_mask.csh regionfile

   regionfile : region file (fk5) removing fake point sources

Use:  bin/make_mask.csh mos1S001-bkg_region-sky_v2.reg


EOF

exit

endif

set reg=$1
if (! -f $reg) then
echo "ERROR no $reg file"
exit
endif

set format=`head -3 $reg | tail -1`
if ( $format != "fk5" ) then
echo "ERROR : region file $reg is not fk5 format"
exit
endif

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`

#----------------------------------------------------------------------------------------
#setup prefix
set prefixlog="log/prefix.log"
set mos1prefix=`head -1 $prefixlog | gawk '{print $1}'`
set mos2prefix=`head -2 $prefixlog | tail -1 |  gawk '{print $1}'`
set pnprefix=`tail -1 $prefixlog | gawk '{print $1}'`

#setup energy range
set eRangeFile="SetUpFile/eRange.dat"
set LOWENE=`head -1 $eRangeFile | gawk '{print $1}'`
set HIGHENE=`head -2 $eRangeFile | gawk 'NR==2 {print $1}'`
#----------------------------------------------------------------------------------------

#the above program overwrites mos1S001-bkg_region-det.fits and so on.
#This procedure is for remake mask files.
cp ${mos1prefix}-bkg_region-det.fits.log ${mos1prefix}-bkg_region-det.fits
cp ${mos1prefix}-bkg_region-sky.fits.log ${mos1prefix}-bkg_region-sky.fits
cp ${mos2prefix}-bkg_region-det.fits.log ${mos2prefix}-bkg_region-det.fits
cp ${mos2prefix}-bkg_region-sky.fits.log ${mos2prefix}-bkg_region-sky.fits
cp ${pnprefix}-bkg_region-det.fits.log ${pnprefix}-bkg_region-det.fits
cp ${pnprefix}-bkg_region-sky.fits.log ${pnprefix}-bkg_region-sky.fits


#remake region fits table
bin/RemakeRegionFile.py $reg


#make_maskで再度マスクイメージを作成する。
echo "   Making mask images...."
make_mask inimage=${mos1prefix}-obj-im.fits inmask=${mos1prefix}-mask-im.fits outmask=${mos1prefix}-cheese.fits reglist=${mos1prefix}-bkg_region-sky.fits
make_mask inimage=${mos2prefix}-obj-im.fits inmask=${mos2prefix}-mask-im.fits outmask=${mos2prefix}-cheese.fits reglist=${mos2prefix}-bkg_region-sky.fits
make_mask inimage=${pnprefix}-obj-im.fits inmask=${pnprefix}-mask-im.fits outmask=${pnprefix}-cheese.fits reglist=${pnprefix}-bkg_region-sky.fits

echo "===================================================="
echo " check output "
echo "===================================================="
echo "ds9 ${mos1prefix}-cheese.fits &"
echo " "
echo "===================================================="
echo " Next run : bin/make_image.csh"
echo "===================================================="
