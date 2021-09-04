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
@ m =  ${#RadiusList}
@ n =  $m - 3
while ($n <= ${#RadiusList} - 2)
@ n1 = $n + 1
@ n2 = $n + 2
@ n3 = $n + 3
echo $n1,$n2,$n3
rm detmapmos1_$n.ds
 evselect table=mos1S001-clean.fits:EVENTS ignorelegallimits=yes imageset=detmapmos1_$n.ds \
   xcolumn=DETX ximagebinsize=100 ycolumn=DETY yimagebinsize=100 squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-mos1-$RadiusList[$n]-$RadiusList[$n1].log

rm mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n1]-$RadiusList[$n2].arf
 arfgen arfset=mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n1]-$RadiusList[$n2].arf \
   spectrumset=mos1S001-obj-$RadiusList[$n1]-$RadiusList[$n2].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[$n1]-$RadiusList[$n2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[$n1]-$RadiusList[$n2].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_$n.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_22.log

rm mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n]-$RadiusList[$n1].arf
arfgen arfset=mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n]-$RadiusList[$n1].arf \
   spectrumset=mos1S001-obj-$RadiusList[$n]-$RadiusList[$n1].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[$n1]-$RadiusList[$n2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[$n]-$RadiusList[$n1].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_$n.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_23.log

rm mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n2]-$RadiusList[$n3].arf
 arfgen arfset=mos1S001-$RadiusList[$n1]-$RadiusList[$n2]-$RadiusList[$n2]-$RadiusList[$n3].arf \
   spectrumset=mos1S001-obj-$RadiusList[$n2]-$RadiusList[$n3].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[$n1]-$RadiusList[$n2].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[$n2]-$RadiusList[$n3].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_$n.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_24.log
   @ n++ 
end

#Final annulus                                                                                                    
echo ${#RadiusList}
@ m = ${#RadiusList}
@ n = $m - 1
@ l = $m - 2
rm detmapmos1_last.ds
 evselect table=mos1S001-clean.fits:EVENTS ignorelegallimits=yes imageset=detmapmos1_last.ds \
   xcolumn=DETX ximagebinsize=200 ycolumn=DETY yimagebinsize=200 squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-mos1-$RadiusList[$n]-$RadiusList[$m].log

rm mos1S001-$RadiusList[$n]-$RadiusList[$m]-$RadiusList[$n]-$RadiusList[$m].arf
 arfgen arfset=mos1S001-$RadiusList[$n]-$RadiusList[$m]-$RadiusList[$n]-$RadiusList[$m].arf \
   spectrumset=mos1S001-obj-$RadiusList[$n]-$RadiusList[$m].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[$n]-$RadiusList[$m].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[$n]-$RadiusList[$m].rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmapmos1_last.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_25.log

rm mos1S001-$RadiusList[$n]-$RadiusList[$m]-$RadiusList[$l]-$RadiusList[$n].arf
 arfgen arfset=mos1S001-$RadiusList[$n]-$RadiusList[$m]-$RadiusList[$l]-$RadiusList[$n].arf \
   spectrumset=mos1S001-obj-$RadiusList[$l]-$RadiusList[$n].pi \
   crossreg_spectrumset=mos1S001-obj-$RadiusList[$n]-$RadiusList[$m].pi crossregionarf=yes withrmfset=yes \
   rmfset=mos1S001-$RadiusList[$l]-$RadiusList[$n].rmf extendedsource=yes modelee=no \
   withbadpixcorr=yes detmaptype=dataset detmaparray=detmapmos1_last.ds badpixlocation=mos1S001-clean.fits \
   modelootcorr=no >& log/arfgen-mos1_26.log

