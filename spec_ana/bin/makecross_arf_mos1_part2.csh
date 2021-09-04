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

#Intermediate annulus
@ n = 3
while ($n <= ${#RadiusList} - 4)
@ n1 = $n + 1
@ n2 = $n + 2
@ n3 = $n + 3
echo $n1,$n2,$n3

rm detmapmos1_$n.ds
evselect table=mos1S001-clean.fits:EVENTS ignorelegallimits=yes imageset=detmapmos1_$n.ds \
   xcolumn=DETX ximagebinsize=50 ycolumn=DETY yimagebinsize=50 squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-mos1-$RadiusList[$n]-$RadiusList[$n1].log

rm mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n1]-$RadiusList[$n2].arf
arfgen arfset=mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n1]-$RadiusList[$n2].arf \
   spectrumset=mos1S001-obj-$RadiusList[$n1]-$RadiusList[$n2].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[$n1]-$RadiusList[$n2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[$n1]-$RadiusList[$n2].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_$n.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_19.log

rm mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n]-$RadiusList[$n1].arf
arfgen arfset=mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n]-$RadiusList[$n1].arf \
   spectrumset=mos1S001-obj-$RadiusList[$n]-$RadiusList[$n1].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[$n1]-$RadiusList[$n2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[$n]-$RadiusList[$n1].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_$n.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_20.log

rm mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n2]-$RadiusList[$n3].arf
arfgen arfset=mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n2]-$RadiusList[$n3].arf \
   spectrumset=mos1S001-obj-$RadiusList[$n2]-$RadiusList[$n3].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[$n1]-$RadiusList[$n2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[$n2]-$RadiusList[$n3].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_$n.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_21.log
   @ n++ 
end

