#!/bin/csh -f

if ( ! -d log ) then
    mkdir log
endif

if ( ! -d SetUpFile ) then
    mkdir SetUpFile
endif

if ( ! -d ../data ) then
    mkdir ../data
endif

if ( ! -d ../image_ana ) then
    mkdir ../image_ana
endif

if ( ! -d ../image_ana/log ) then
    mkdir ../image_ana/log
endif

if ( ! -d ../cal_EI ) then
    mkdir ../cal_EI
endif

if ( ! -d ../result ) then
    mkdir ../result
endif

if ( ! -d ../radial_fit ) then
    mkdir ../radial_fit
endif


source ${SAS_DIR}/setsas.csh >& /dev/null
#--------------------------
#setup odf
set specpath=$PWD
setenv SAS_ODF ${specpath}/../odf
#---------------------------


cifbuild >& log/00cifbuild.log

setenv SAS_CCF ccf.cif 

odfingest >& log/01odfingest.log

setenv SAS_ODF `ls -1 *SUM.SAS`
setenv SAS_RAND_SEED 0

epchain >& log/02epchain.log
epchain withoutoftime=true >& log/03epchain.log
emchain >& log/04emchain.log
pn-filter >& log/05pn-filter.log
mos-filter >& log/06mos-filter.log

echo "---- Setup Prefix -----"
#setup prefix
set mos1=`more log/06mos-filter.log | grep eventset | gawk '{if($1=="espfilt") print $0}' | cut -d"." -f1 | sed 's/mos/ /g' | gawk '{print $3}' | head -1`
set mos2=`more log/06mos-filter.log | grep eventset | gawk '{if($1=="espfilt") print $0}' | cut -d"." -f1 | sed 's/mos/ /g' | gawk '{print $3}' | tail -1`
set pn=`more log/05pn-filter.log | grep eventset | gawk '{if($1=="espfilt") print $0}' | cut -d"." -f1 | sed 's/pn/ /g' | gawk '{print $3}' | tail -1`

set prefixlog="log/prefix.log"
echo "mos"$mos1" "$mos1 > $prefixlog
echo "mos"$mos2" "$mos2 >> $prefixlog
echo "pn"$pn" "$pn >> $prefixlog



echo "======================================================="
echo "Next run : "
echo "bin/cheese.csh scale rate dist"
echo "======================================================="

