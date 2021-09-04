#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`
#---------------------------------------------------------
#setup prefix
set prefixlog="../spec_ana/log/prefix.log"
set mos1prefix=`head -1 $prefixlog | gawk '{print $1}'`
set mos2prefix=`head -2 $prefixlog | tail -1 |  gawk '{print $1}'`
set pnprefix=`tail -1 $prefixlog | gawk '{print $1}'`
#---------------------------------------------------------
#setup energy range
set eRangeFile="../spec_ana/SetUpFile/eRange.dat"
set LOWENE=`head -1 $eRangeFile | gawk '{print $1}'`
set HIGHENE=`head -2 $eRangeFile | gawk 'NR==2 {print $1}'`
#----------------------------------------------------------------
@ lowene=`expr $LOWENE`
@ highene=`expr $HIGHENE`

@ CenterENE = ($lowene + $highene)/ 2

echo $CenterENE

#DET座標の銀河団中心位置を入れる。

set centerlog="../spec_ana/log/Center.log"

set CentX_MOS1=`head -3 $centerlog | tail -1 | gawk '{print $5}'`
set CentY_MOS1=`head -3 $centerlog | tail -1 | gawk '{print $6}'`

set CentX_MOS2=`head -4 $centerlog | tail -1 | gawk '{print $5}'`
set CentY_MOS2=`head -4 $centerlog | tail -1 | gawk '{print $6}'`

set CentX_PN=`head -5 $centerlog | tail -1 | gawk '{print $5}'`
set CentY_PN=`head -5 $centerlog | tail -1 | gawk '{print $6}'`

echo $CentX_MOS1
echo $CentX_MOS2
echo "ok"
echo $CentX_PN
echo $CentY_PN
echo "ok"

foreach prefix ( "mos1" "mos2" "pn" )
if ( -f ${prefix}_psf.fits) then
rm ${prefix}_psf.fits
endif
end

psfgen region="(DETX, DETY) IN box(${CentX_MOS1},${CentY_MOS1},2000,2000)" instrument=M1 output=mos1_psf.fits level=ELLBETA energy=${CenterENE} coortype=XY image=${mos1prefix}-clean.fits
psfgen region="(DETX, DETY) IN box(${CentX_MOS2},${CentY_MOS2},2000,2000)" instrument=M2 output=mos2_psf.fits level=ELLBETA energy=${CenterENE} coortype=XY image=${mos2prefix}-clean.fits
psfgen region="(DETX, DETY) IN box(${CentX_PN},${CentY_PN},2000,2000)" instrument=PN output=pn_psf.fits level=ELLBETA energy=${CenterENE} coortype=XY image=${pnprefix}-clean.fits

foreach prefix ( "mos1" "mos2" "pn" )
if ( -f ${prefix}_psf_2bin.fits) then
rm ${prefix}_psf_2bin.fits
endif
end

fimgbin mos1_psf.fits mos1_psf_2bin.fits 2
fimgbin mos2_psf.fits mos2_psf_2bin.fits 2
fimgbin pn_psf.fits pn_psf_2bin.fits 2

#---------------------------------------------------------------------------------
#check output files
foreach prefix ( "mos1" "mos2" "pn" )
if (! -f ${prefix}_psf.fits) then
echo "ERROR : no output file " ${prefix}_psf.fits
exit 
endif
if (! -f ${prefix}_psf_2bin.fits) then
echo "ERROR : no output file " ${prefix}_psf_2bin.fits
exit 
endif
end


cp *psf* ../data/

echo "======================================================="
echo "Success  : all files are completed"
echo "Next run : "
echo "cd ../cal_EI"
echo "bin/Make_fake.csh z"     
echo "======================================================="

