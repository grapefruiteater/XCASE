#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`

if ($#argv != 1) then
    echo "#----------------usage----------------#"
    echo "run_grppha.csh grpbin"
    echo ""
    echo "grppha and copy outputfiles to ../data"
    echo "grpbin : bin for grppha"
    echo ""
    echo "Example :"
    echo "bin/run_grppha.csh 50"
    echo "#-------------------------------------#"

exit
endif

set grpbin=$1

set RegionFile='SetUpFile/Region.dat'
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

#setup energy range
set eRangeFile="SetUpFile/eRange.dat"
set LOWENE=`head -1 $eRangeFile | gawk '{print $1}'`
set HIGHENE=`head -2 $eRangeFile | gawk 'NR==2 {print $1}'`

#setup group bin 
set GrpBinFile="SetUpFile/GrpBin.dat"
echo $grpbin > $GrpBinFile
#----------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
#銀河団中心から同心円円環領域のスペクトルとレスポンスを作成

set RediusList = ()
set COUNT = `more $RegionFile | wc -l`
@ i = 1
while($i <= $COUNT)
    set Radius = "`cat $RegionFile | head -$i |tail -1`"
    set RediusList = ($RediusList $Radius)
    @ i++
end

set i=0
@ max = ${#RediusList} - 1

while($i<$max)
@ i = $i + 1
set inner = $RediusList[$i]
@ j = $i + 1
set outer = $RediusList[$j]
set regionname = $inner-$outer
echo $regionname
if (-f ${mos1prefix}-obj-${regionname}-grp.pi) then
rm ${mos1prefix}-obj-${regionname}-grp.pi
endif
grppha ${mos1prefix}-obj-${regionname}.pi ${mos1prefix}-obj-${regionname}-grp.pi <<EOF
chkey BACKFILE ${mos1prefix}-back-${regionname}.pi 
chkey RESPFILE ${mos1prefix}-${regionname}.rmf 
chkey ANCRFILE ${mos1prefix}-${regionname}.arf 
group min ${grpbin}
show group
exit
EOF
end

set i=0
while($i<$max)
@ i = $i + 1 
set inner = $RediusList[$i]
@ j = $i + 1
set outer = $RediusList[$j]
set regionname = $inner-$outer
echo $regionname
if ( -f ${mos2prefix}-obj-${regionname}-grp.pi ) then
rm ${mos2prefix}-obj-${regionname}-grp.pi
endif
grppha ${mos2prefix}-obj-${regionname}.pi ${mos2prefix}-obj-${regionname}-grp.pi <<EOF
chkey BACKFILE ${mos2prefix}-back-${regionname}.pi 
chkey RESPFILE ${mos2prefix}-${regionname}.rmf 
chkey ANCRFILE ${mos2prefix}-${regionname}.arf 
group min ${grpbin}
exit
EOF
end

set i=0
while($i<$max)
@ i = $i + 1 
set inner = $RediusList[$i]
@ j = $i + 1
set outer = $RediusList[$j]
set regionname = $inner-$outer
echo $regionname
if( -f ${pnprefix}-obj-os-${regionname}-grp.pi ) then
rm ${pnprefix}-obj-os-${regionname}-grp.pi
endif
grppha ${pnprefix}-obj-${regionname}.pi ${pnprefix}-obj-os-${regionname}-grp.pi <<EOF
chkey BACKFILE ${pnprefix}-back-${regionname}.pi 
chkey RESPFILE ${pnprefix}-${regionname}.rmf 
chkey ANCRFILE ${pnprefix}-${regionname}.arf 
group min ${grpbin}
exit
EOF
end

#copy spec file

cp *obj**.pi *.arf *.rmf ../data/
