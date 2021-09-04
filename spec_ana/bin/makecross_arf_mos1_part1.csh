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
rm detmapmos1_1.ds
 evselect table=mos1S001-clean.fits:EVENTS ignorelegallimits=yes imageset=detmapmos1_1.ds \
   xcolumn=DETX ximagebinsize=50 ycolumn=DETY yimagebinsize=50 squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-mos1-$RadiusList[1]-$RadiusList[2].log

rm mos1S001-$RadiusList[1]-$RadiusList[2]-$RadiusList[1]-$RadiusList[2].arf
 arfgen arfset=mos1S001-$RadiusList[1]-$RadiusList[2]-$RadiusList[1]-$RadiusList[2].arf \
   spectrumset=mos1S001-obj-$RadiusList[1]-$RadiusList[2].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[1]-$RadiusList[2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[1]-$RadiusList[2].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_1.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_1.log
rm mos1S001-$RadiusList[1]-$RadiusList[2]-$RadiusList[2]-$RadiusList[3].arf
 arfgen arfset=mos1S001-$RadiusList[1]-$RadiusList[2]-$RadiusList[2]-$RadiusList[3].arf \
   spectrumset=mos1S001-obj-$RadiusList[2]-$RadiusList[3].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[1]-$RadiusList[2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[2]-$RadiusList[3].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_1.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_2.log
rm mos1S001-$RadiusList[1]-$RadiusList[2]-$RadiusList[3]-$RadiusList[4].arf
 arfgen arfset=mos1S001-$RadiusList[1]-$RadiusList[2]-$RadiusList[3]-$RadiusList[4].arf \
   spectrumset=mos1S001-obj-$RadiusList[3]-$RadiusList[4].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[1]-$RadiusList[2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[3]-$RadiusList[4].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_1.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_3.log
rm mos1S001-$RadiusList[2]-$RadiusList[3]-$RadiusList[2]-$RadiusList[3].arf
 arfgen arfset=mos1S001-$RadiusList[2]-$RadiusList[3]-$RadiusList[2]-$RadiusList[3].arf \
   spectrumset=mos1S001-obj-$RadiusList[2]-$RadiusList[3].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[2]-$RadiusList[3].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[2]-$RadiusList[3].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_1.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_4.log
rm mos1S001-$RadiusList[2]-$RadiusList[3]-$RadiusList[1]-$RadiusList[2].arf
 arfgen arfset=mos1S001-$RadiusList[2]-$RadiusList[3]-$RadiusList[1]-$RadiusList[2].arf \
   spectrumset=mos1S001-obj-$RadiusList[1]-$RadiusList[2].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[2]-$RadiusList[3].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[1]-$RadiusList[2].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_1.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_5.log
rm mos1S001-$RadiusList[2]-$RadiusList[3]-$RadiusList[3]-$RadiusList[4].arf
 arfgen arfset=mos1S001-$RadiusList[2]-$RadiusList[3]-$RadiusList[3]-$RadiusList[4].arf \
   spectrumset=mos1S001-obj-$RadiusList[3]-$RadiusList[4].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[2]-$RadiusList[3].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[3]-$RadiusList[4].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_1.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_6.log

#Second region
 rm detmapmos1_2.ds
 evselect table=mos1S001-clean.fits:EVENTS ignorelegallimits=yes imageset=detmapmos1_2.ds \
   xcolumn=DETX ximagebinsize=70 ycolumn=DETY yimagebinsize=70 squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-mos1-$RadiusList[2]-$RadiusList[3].log

rm mos1S001-$RadiusList[1]-$RadiusList[2]-$RadiusList[4]-$RadiusList[5].arf
 arfgen arfset=mos1S001-$RadiusList[1]-$RadiusList[2]-$RadiusList[4]-$RadiusList[5].arf \
   spectrumset=mos1S001-obj-$RadiusList[4]-$RadiusList[5].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[1]-$RadiusList[2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[4]-$RadiusList[5].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_2.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no  >& log/arfgen-mos1_7.log
rm mos1S001-$RadiusList[2]-$RadiusList[3]-$RadiusList[4]-$RadiusList[5].arf
 arfgen arfset=mos1S001-$RadiusList[2]-$RadiusList[3]-$RadiusList[4]-$RadiusList[5].arf \
   spectrumset=mos1S001-obj-$RadiusList[4]-$RadiusList[5].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[2]-$RadiusList[3].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[4]-$RadiusList[5].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_2.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_8.log
rm mos1S001-$RadiusList[3]-$RadiusList[4]-$RadiusList[3]-$RadiusList[4].arf
 arfgen arfset=mos1S001-$RadiusList[3]-$RadiusList[4]-$RadiusList[3]-$RadiusList[4].arf \
   spectrumset=mos1S001-obj-$RadiusList[3]-$RadiusList[4].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[3]-$RadiusList[4].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[3]-$RadiusList[4].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_2.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no  >& log/arfgen-mos1_9.log
rm mos1S001-$RadiusList[3]-$RadiusList[4]-$RadiusList[1]-$RadiusList[2].arf
 arfgen arfset=mos1S001-$RadiusList[3]-$RadiusList[4]-$RadiusList[1]-$RadiusList[2].arf \
   spectrumset=mos1S001-obj-$RadiusList[1]-$RadiusList[2].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[3]-$RadiusList[4].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[1]-$RadiusList[2].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_2.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_10.log
rm mos1S001-$RadiusList[3]-$RadiusList[4]-$RadiusList[2]-$RadiusList[3].arf
 arfgen arfset=mos1S001-$RadiusList[3]-$RadiusList[4]-$RadiusList[2]-$RadiusList[3].arf \
   spectrumset=mos1S001-obj-$RadiusList[2]-$RadiusList[3].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[3]-$RadiusList[4].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[2]-$RadiusList[3].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_2.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no  >& log/arfgen-mos1_11.log
rm mos1S001-$RadiusList[3]-$RadiusList[4]-$RadiusList[4]-$RadiusList[5].arf
 arfgen arfset=mos1S001-$RadiusList[3]-$RadiusList[4]-$RadiusList[4]-$RadiusList[5].arf \
   spectrumset=mos1S001-obj-$RadiusList[4]-$RadiusList[5].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[3]-$RadiusList[4].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[4]-$RadiusList[5].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_2.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no  >& log/arfgen-mos1_12.log


#Third region
rm detmapmos1_3.ds
evselect table=mos1S001-clean.fits:EVENTS ignorelegallimits=yes imageset=detmapmos1_3.ds \
   xcolumn=DETX ximagebinsize=70 ycolumn=DETY yimagebinsize=70 squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-mos1-$RadiusList[3]-$RadiusList[4].log

rm mos1S001-$RadiusList[1]-$RadiusList[2]-$RadiusList[5]-$RadiusList[6].arf
 arfgen arfset=mos1S001-$RadiusList[1]-$RadiusList[2]-$RadiusList[5]-$RadiusList[6].arf \
   spectrumset=mos1S001-obj-$RadiusList[5]-$RadiusList[6].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[1]-$RadiusList[2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[5]-$RadiusList[6].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_3.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no  >& log/arfgen-mos1_13.log
rm mos1S001-$RadiusList[2]-$RadiusList[3]-$RadiusList[5]-$RadiusList[6].arf
 arfgen arfset=mos1S001-$RadiusList[2]-$RadiusList[3]-$RadiusList[5]-$RadiusList[6].arf \
   spectrumset=mos1S001-obj-$RadiusList[5]-$RadiusList[6].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[2]-$RadiusList[3].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[5]-$RadiusList[6].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_3.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_14.log
rm mos1S001-$RadiusList[4]-$RadiusList[5]-$RadiusList[4]-$RadiusList[5].arf
 arfgen arfset=mos1S001-$RadiusList[4]-$RadiusList[5]-$RadiusList[4]-$RadiusList[5].arf \
   spectrumset=mos1S001-obj-$RadiusList[4]-$RadiusList[5].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[4]-$RadiusList[5].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[4]-$RadiusList[5].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_3.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no  >& log/arfgen-mos1_15.log
rm mos1S001-$RadiusList[4]-$RadiusList[5]-$RadiusList[1]-$RadiusList[2].arf
 arfgen arfset=mos1S001-$RadiusList[4]-$RadiusList[5]-$RadiusList[1]-$RadiusList[2].arf \
   spectrumset=mos1S001-obj-$RadiusList[1]-$RadiusList[2].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[4]-$RadiusList[5].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[1]-$RadiusList[2].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_3.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_16.log
rm mos1S001-$RadiusList[4]-$RadiusList[5]-$RadiusList[2]-$RadiusList[3].arf
 arfgen arfset=mos1S001-$RadiusList[4]-$RadiusList[5]-$RadiusList[2]-$RadiusList[3].arf \
   spectrumset=mos1S001-obj-$RadiusList[2]-$RadiusList[3].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[4]-$RadiusList[5].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[2]-$RadiusList[3].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_3.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no  >& log/arfgen-mos1_17.log
rm mos1S001-$RadiusList[3]-$RadiusList[4]-$RadiusList[5]-$RadiusList[6].arf
 arfgen arfset=mos1S001-$RadiusList[3]-$RadiusList[4]-$RadiusList[5]-$RadiusList[6].arf \
   spectrumset=mos1S001-obj-$RadiusList[5]-$RadiusList[6].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[3]-$RadiusList[4].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[5]-$RadiusList[6].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_3.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no  >& log/arfgen-mos1_18.log
