#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`
#======================================================================================================

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

#----------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
if ( ! -f $RegionFile ) then
     echo ""
     echo "ERROR : NO "$RegionFile
     echo ""
     exit
endif

echo "# Radii for regions"                                                                                                                                        
set RadiusList = ()                                                                                                                                               
set COUNT = `more $RegionFile | wc -l`                                                                                                                            
@ i = 1                                                                                                                                                           
while($i <= $COUNT)                                                                                                                                               
    set Radius = "`cat $RegionFile | head -$i |tail -1`"                                                                                                          
    echo "$Radius"                                                                                                                                                
    set RadiusList = ($RadiusList $Radius)                                                                                                                        
    @ i++                                                                                                                                                         
end                                                                                                                                                               
echo "------------------------------"                                                                                                                              

#Making Cross ARF
#First region
rm detmapmos2_1.ds
 evselect table=mos2S002-clean.fits:EVENTS ignorelegallimits=yes imageset=detmapmos2_1.ds \
   xcolumn=DETX ximagebinsize=50 ycolumn=DETY yimagebinsize=50 squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-mos2-$RadiusList[1]-$RadiusList[2].log

rm mos2S002-$RadiusList[1]-$RadiusList[2]-$RadiusList[1]-$RadiusList[2].arf
 arfgen arfset=mos2S002-$RadiusList[1]-$RadiusList[2]-$RadiusList[1]-$RadiusList[2].arf \
   spectrumset=mos2S002-obj-$RadiusList[1]-$RadiusList[2].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[1]-$RadiusList[2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[1]-$RadiusList[2].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_1.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no >& log/arfgen-mos2_1.log
rm mos2S002-$RadiusList[1]-$RadiusList[2]-$RadiusList[2]-$RadiusList[3].arf
 arfgen arfset=mos2S002-$RadiusList[1]-$RadiusList[2]-$RadiusList[2]-$RadiusList[3].arf \
   spectrumset=mos2S002-obj-$RadiusList[2]-$RadiusList[3].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[1]-$RadiusList[2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[2]-$RadiusList[3].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_1.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no >& log/arfgen-mos2_2.log
rm mos2S002-$RadiusList[1]-$RadiusList[2]-$RadiusList[3]-$RadiusList[4].arf
 arfgen arfset=mos2S002-$RadiusList[1]-$RadiusList[2]-$RadiusList[3]-$RadiusList[4].arf \
   spectrumset=mos2S002-obj-$RadiusList[3]-$RadiusList[4].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[1]-$RadiusList[2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[3]-$RadiusList[4].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_1.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no >& log/arfgen-mos2_3.log
rm mos2S002-$RadiusList[2]-$RadiusList[3]-$RadiusList[2]-$RadiusList[3].arf
 arfgen arfset=mos2S002-$RadiusList[2]-$RadiusList[3]-$RadiusList[2]-$RadiusList[3].arf \
   spectrumset=mos2S002-obj-$RadiusList[2]-$RadiusList[3].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[2]-$RadiusList[3].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[2]-$RadiusList[3].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_1.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no >& log/arfgen-mos2_4.log
rm mos2S002-$RadiusList[2]-$RadiusList[3]-$RadiusList[1]-$RadiusList[2].arf
 arfgen arfset=mos2S002-$RadiusList[2]-$RadiusList[3]-$RadiusList[1]-$RadiusList[2].arf \
   spectrumset=mos2S002-obj-$RadiusList[1]-$RadiusList[2].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[2]-$RadiusList[3].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[1]-$RadiusList[2].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_1.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no >& log/arfgen-mos2_5.log
rm mos2S002-$RadiusList[2]-$RadiusList[3]-$RadiusList[3]-$RadiusList[4].arf
 arfgen arfset=mos2S002-$RadiusList[2]-$RadiusList[3]-$RadiusList[3]-$RadiusList[4].arf \
   spectrumset=mos2S002-obj-$RadiusList[3]-$RadiusList[4].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[2]-$RadiusList[3].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[3]-$RadiusList[4].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_1.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no >& log/arfgen-mos2_6.log

#Second region
 rm detmapmos2_2.ds
 evselect table=mos2S002-clean.fits:EVENTS ignorelegallimits=yes imageset=detmapmos2_2.ds \
   xcolumn=DETX ximagebinsize=70 ycolumn=DETY yimagebinsize=70 squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-mos2-$RadiusList[2]-$RadiusList[3].log

rm mos2S002-$RadiusList[1]-$RadiusList[2]-$RadiusList[4]-$RadiusList[5].arf
 arfgen arfset=mos2S002-$RadiusList[1]-$RadiusList[2]-$RadiusList[4]-$RadiusList[5].arf \
   spectrumset=mos2S002-obj-$RadiusList[4]-$RadiusList[5].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[1]-$RadiusList[2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[4]-$RadiusList[5].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_2.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no  >& log/arfgen-mos2_7.log
rm mos2S002-$RadiusList[2]-$RadiusList[3]-$RadiusList[4]-$RadiusList[5].arf
 arfgen arfset=mos2S002-$RadiusList[2]-$RadiusList[3]-$RadiusList[4]-$RadiusList[5].arf \
   spectrumset=mos2S002-obj-$RadiusList[4]-$RadiusList[5].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[2]-$RadiusList[3].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[4]-$RadiusList[5].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_2.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no >& log/arfgen-mos2_8.log
rm mos2S002-$RadiusList[3]-$RadiusList[4]-$RadiusList[3]-$RadiusList[4].arf
 arfgen arfset=mos2S002-$RadiusList[3]-$RadiusList[4]-$RadiusList[3]-$RadiusList[4].arf \
   spectrumset=mos2S002-obj-$RadiusList[3]-$RadiusList[4].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[3]-$RadiusList[4].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[3]-$RadiusList[4].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_2.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no  >& log/arfgen-mos2_9.log
rm mos2S002-$RadiusList[3]-$RadiusList[4]-$RadiusList[1]-$RadiusList[2].arf
 arfgen arfset=mos2S002-$RadiusList[3]-$RadiusList[4]-$RadiusList[1]-$RadiusList[2].arf \
   spectrumset=mos2S002-obj-$RadiusList[1]-$RadiusList[2].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[3]-$RadiusList[4].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[1]-$RadiusList[2].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_2.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no >& log/arfgen-mos2_10.log
rm mos2S002-$RadiusList[3]-$RadiusList[4]-$RadiusList[2]-$RadiusList[3].arf
 arfgen arfset=mos2S002-$RadiusList[3]-$RadiusList[4]-$RadiusList[2]-$RadiusList[3].arf \
   spectrumset=mos2S002-obj-$RadiusList[2]-$RadiusList[3].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[3]-$RadiusList[4].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[2]-$RadiusList[3].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_2.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no  >& log/arfgen-mos2_11.log
rm mos2S002-$RadiusList[3]-$RadiusList[4]-$RadiusList[4]-$RadiusList[5].arf
 arfgen arfset=mos2S002-$RadiusList[3]-$RadiusList[4]-$RadiusList[4]-$RadiusList[5].arf \
   spectrumset=mos2S002-obj-$RadiusList[4]-$RadiusList[5].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[3]-$RadiusList[4].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[4]-$RadiusList[5].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_2.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no  >& log/arfgen-mos2_12.log


#Third region
rm detmapmos2_3.ds
evselect table=mos2S002-clean.fits:EVENTS ignorelegallimits=yes imageset=detmapmos2_3.ds \
   xcolumn=DETX ximagebinsize=70 ycolumn=DETY yimagebinsize=70 squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-mos2-$RadiusList[3]-$RadiusList[4].log

rm mos2S002-$RadiusList[1]-$RadiusList[2]-$RadiusList[5]-$RadiusList[6].arf
 arfgen arfset=mos2S002-$RadiusList[1]-$RadiusList[2]-$RadiusList[5]-$RadiusList[6].arf \
   spectrumset=mos2S002-obj-$RadiusList[5]-$RadiusList[6].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[1]-$RadiusList[2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[5]-$RadiusList[6].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_3.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no  >& log/arfgen-mos2_13.log
rm mos2S002-$RadiusList[2]-$RadiusList[3]-$RadiusList[5]-$RadiusList[6].arf
 arfgen arfset=mos2S002-$RadiusList[2]-$RadiusList[3]-$RadiusList[5]-$RadiusList[6].arf \
   spectrumset=mos2S002-obj-$RadiusList[5]-$RadiusList[6].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[2]-$RadiusList[3].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[5]-$RadiusList[6].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_3.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no >& log/arfgen-mos2_14.log
rm mos2S002-$RadiusList[4]-$RadiusList[5]-$RadiusList[4]-$RadiusList[5].arf
 arfgen arfset=mos2S002-$RadiusList[4]-$RadiusList[5]-$RadiusList[4]-$RadiusList[5].arf \
   spectrumset=mos2S002-obj-$RadiusList[4]-$RadiusList[5].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[4]-$RadiusList[5].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[4]-$RadiusList[5].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_3.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no  >& log/arfgen-mos2_15.log
rm mos2S002-$RadiusList[4]-$RadiusList[5]-$RadiusList[1]-$RadiusList[2].arf
 arfgen arfset=mos2S002-$RadiusList[4]-$RadiusList[5]-$RadiusList[1]-$RadiusList[2].arf \
   spectrumset=mos2S002-obj-$RadiusList[1]-$RadiusList[2].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[4]-$RadiusList[5].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[1]-$RadiusList[2].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_3.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no >& log/arfgen-mos2_16.log
rm mos2S002-$RadiusList[4]-$RadiusList[5]-$RadiusList[2]-$RadiusList[3].arf
 arfgen arfset=mos2S002-$RadiusList[4]-$RadiusList[5]-$RadiusList[2]-$RadiusList[3].arf \
   spectrumset=mos2S002-obj-$RadiusList[2]-$RadiusList[3].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[4]-$RadiusList[5].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[2]-$RadiusList[3].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_3.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no  >& log/arfgen-mos2_17.log
rm mos2S002-$RadiusList[3]-$RadiusList[4]-$RadiusList[5]-$RadiusList[6].arf
 arfgen arfset=mos2S002-$RadiusList[3]-$RadiusList[4]-$RadiusList[5]-$RadiusList[6].arf \
   spectrumset=mos2S002-obj-$RadiusList[5]-$RadiusList[6].pi \
   crossreg_spectrumset=mos2S002-obj-$RadiusList[3]-$RadiusList[4].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos2S002-$RadiusList[5]-$RadiusList[6].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos2_3.ds badpixlocation=mos2S002-clean.fits \
   modelootcorr=no  >& log/arfgen-mos2_18.log
