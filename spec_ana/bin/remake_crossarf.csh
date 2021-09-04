#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`
#======================================================================================================
if ($#argv != 3) then
    echo "#----------------usage----------------#"
    echo "csh remake_crossarf.csh prefix imagebin rinn1-rout1-rinn2-rout2 "
    echo ""
    echo "prefix   : mos1 or mos2 or pn"
    echo "imagebin : bin for count image"
    echo "rinn1-rout1-rinn2-rout"
    echo "    rinn1 : rint for leaking out "
    echo "    rinn1 : rout for leaking out"
    echo "    rinn2 : rinn for leaking in"
    echo "    rinn2 : rout for leaking in"
    echo ""
    echo "Example :"
    echo "bin/remake_crossarf.csh mos1 50 0-40-40-60"
    echo "#-------------------------------------#"
exit
endif
set prefix0=$1
set imagebin=$2
set makeregion=$3

set regionnameout=`echo $makeregion | cut -d "-" -f 1-2`
set regionnamein=`echo $makeregion | cut -d "-" -f 3-4`
echo $regionnameout
echo $regionnamein
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
echo "  "                                                                                                                              
#Making Cross ARF

if ( $prefix0 == "mos1") then
set prefix=$mos1prefix 
echo "  remaking : "${prefix}-${regionnameout}-${regionnamein}.arf
if ( -f detmap${prefix}_remake.ds ) then
    rm detmap${prefix}_remake.ds
endif
 evselect table=$prefix-clean.fits:EVENTS ignorelegallimits=yes imageset=detmap${prefix0}_remake.ds \
   xcolumn=DETX ximagebinsize=$imagebin ycolumn=DETY yimagebinsize=$imagebin squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-$prefix0-${regionnamein}_remake.log
if ( -f ${prefix}-${regionnameout}-${regionnamein}.arf ) then
    rm ${prefix}-${regionnameout}-${regionnamein}.arf
endif
 arfgen arfset=${prefix}-${regionnameout}-${regionnamein}.arf \
   spectrumset=${prefix}-obj-${regionnamein}.pi \
   crossreg_spectrumset=${prefix}-obj-${regionnameout}.pi crossregionarf=yes withrmfset=yes \
   rmfset=${prefix}-${regionnamein}.rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmap${prefix0}_remake.ds badpixlocation=${prefix}-clean.fits \
   modelootcorr=no >& log/arfgen-${prefix0}_remake.log
endif

if ( $prefix0 == "mos2") then
set prefix=$mos2prefix 
echo "  remaking : "${prefix}-${regionnameout}-${regionnamein}.arf
if ( -f detmap${prefix}_remake.ds ) then
    rm detmap${prefix}_remake.ds
endif
 evselect table=$prefix-clean.fits:EVENTS ignorelegallimits=yes imageset=detmap${prefix0}_remake.ds \
   xcolumn=DETX ximagebinsize=$imagebin ycolumn=DETY yimagebinsize=$imagebin squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-$prefix0-${regionnamein}_remake.log
if ( -f ${prefix}-${regionnameout}-${regionnamein}.arf ) then
    rm ${prefix}-${regionnameout}-${regionnamein}.arf
endif
 arfgen arfset=${prefix}-${regionnameout}-${regionnamein}.arf \
   spectrumset=${prefix}-obj-${regionnamein}.pi \
   crossreg_spectrumset=${prefix}-obj-${regionnameout}.pi crossregionarf=yes withrmfset=yes \
   rmfset=${prefix}-${regionnamein}.rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmap${prefix0}_remake.ds badpixlocation=${prefix}-clean.fits \
   modelootcorr=no >& log/arfgen-${prefix0}_remake.log
endif

if ( $prefix0 == "pn") then
set prefix=$pnprefix 
echo "  remaking : "${prefix}-${regionnameout}-${regionnamein}.arf
if ( -f detmap${prefix}_remake.ds ) then
    rm detmap${prefix}_remake.ds
endif
 evselect table=$prefix-clean.fits:EVENTS ignorelegallimits=yes imageset=detmap${prefix0}_remake.ds \
   xcolumn=DETX ximagebinsize=$imagebin ycolumn=DETY yimagebinsize=$imagebin squarepixels=yes \
   imagebinning=binSize withimageset=yes >& log/evselect-$prefix0-${regionnamein}_remake.log
if ( -f ${prefix}-${regionnameout}-${regionnamein}.arf ) then
    rm ${prefix}-${regionnameout}-${regionnamein}.arf
endif
 arfgen arfset=${prefix}-${regionnameout}-${regionnamein}.arf \
   spectrumset=${prefix}-obj-${regionnamein}.pi \
   crossreg_spectrumset=${prefix}-obj-${regionnameout}.pi crossregionarf=yes withrmfset=yes \
   rmfset=${prefix}-${regionnamein}.rmf extendedsource=yes modelee=no withbadpixcorr=yes \
   detmaptype=dataset detmaparray=detmap${prefix0}_remake.ds badpixlocation=${prefix}-clean.fits \
   modelootcorr=no >& log/arfgen-${prefix0}_remake.log
endif
