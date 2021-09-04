#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`

if ( ! -d log ) then
    mkdir log
endif


#-----------------------------------------------------------------------------------------
#setup prefix
set prefixlog="../spec_ana/log/prefix.log"
set mos1prefix=`head -1 $prefixlog | gawk '{print $1}'`
set mos2prefix=`head -2 $prefixlog | tail -1 |  gawk '{print $1}'`
set pnprefix=`tail -1 $prefixlog | gawk '{print $1}'`

set mos1input=`head -1 $prefixlog | gawk '{print $2}'`
set mos2input=`head -2 $prefixlog | tail -1 |  gawk '{print $2}'`
set pninput=`tail -1 $prefixlog | gawk '{print $2}'`
set mosinput=`echo "'"$mos1input" "$mos2input"'"`

#setup energy range
set eRangeFile="../spec_ana/SetUpFile/eRange.dat"
set LOWENE=`head -1 $eRangeFile | gawk '{print $1}'`
set HIGHENE=`head -2 $eRangeFile | gawk 'NR==2 {print $1}'`
#-------------------------------------------------------------------------------------------
#sp_partialコマンドの出力 log/sp_partial-mos1.log から、”Scaled Normalization” の値を入れる。
set mos1pnorm = `gawk '{if($1=="Scaled" && $2=="Normalization:") print $3}' ../spec_ana/log/sp_partial-mos1.log`
set mos2pnorm = `gawk '{if($1=="Scaled" && $2=="Normalization:") print $3}' ../spec_ana/log/sp_partial-mos2.log`
set pnpnorm   = `gawk '{if($1=="Scaled" && $2=="Normalization:") print $3}' ../spec_ana/log/sp_partial-pn.log`

echo "#--------------------------------------"
echo "# Scaled normalization for soft prton"
echo "mos1 :"$mos1pnorm
echo "mos2 :"$mos2pnorm
echo "pn   :"$pnpnorm
echo "#--------------------------------------"

set inputlog="log/20sp_partial_input.log"
echo "# Scaled normalization for soft prton" > $inputlog
echo "mos1 :"$mos1pnorm >> $inputlog
echo "mos2 :"$mos2pnorm >> $inputlog
echo "pn   :"$pnpnorm >> $inputlog 
echo "#----------------------------------" >> $inputlog
#-------------------------------------------------------------------------------------------
#スペクトルフィットから得られるソフトプロトンのpowerlawのphotoindexのベストフィットパラメータを入れる。
set xspecoutput=`gawk '{if($1=="Best") print "../spec_ana/"$3}' ../spec_ana/log/sp_partial-inputfile.log`

set Nmos1=`gawk '{if($1=="model" &&  $2=="2:mos1_1" && $3=="powerlaw") print NR+1}' $xspecoutput`
set Nmos2=`gawk '{if($1=="model" &&  $2=="3:mos2_1" && $3=="powerlaw") print NR+1}' $xspecoutput`
set Npn=`gawk '{if($1=="model" &&  $2=="4:pn_1" && $3=="powerlaw") print NR+1}' $xspecoutput`

set mos1pindex=`head -$Nmos1 $xspecoutput | tail -1 | gawk '{print $1}'` 
set mos2pindex=$mos1pindex
set pnpindex=`head -$Npn $xspecoutput | tail -1 | gawk '{print $1}'`

#---------------------------------------
echo "# PowerLaw for soft proton" >> $inputlog
echo "mos1 :"$mos1pindex >> $inputlog
echo "mos2 :"$mos2pindex >> $inputlog
echo "pn   :"$pnpindex >> $inputlog
echo "#--------------------------------------" >> $inputlog

#-------------------------------------------------------------------------------------------------------

foreach prefix ( $mos1prefix $mos2prefix $pnprefix )

if (! -f ${prefix}-obj-full-grp.pi) then
echo "ERROR : no file of "${prefix}-obj-full-grp.pi
exit 
endif
if (-z ${prefix}-obj-full-grp.pi) then
echo "ERROR : "${prefix}-obj-full-grp.pi" exists but zero size" 
exit 
endif

end

proton prefix=$mos1input caldb=$SAS_ESAS_CALDB specname=${mos1prefix}-obj-full-grp.pi ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 elow=$LOWENE ehigh=$HIGHENE spectrumcontrol=1 pindex=$mos1pindex pnorm=$mos1pnorm clobber=1 >& log/20sp_partial_mos1.log
proton prefix=$mos2input caldb=$SAS_ESAS_CALDB specname=${mos2prefix}-obj-full-grp.pi ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 elow=$LOWENE ehigh=$HIGHENE spectrumcontrol=1 pindex=$mos2pindex pnorm=$mos2pnorm clobber=1 >& log/20sp_partial_mos2.log
proton prefix=$pninput caldb=$SAS_ESAS_CALDB specname=${pnprefix}-obj-full-grp.pi  elow=$LOWENE ehigh=$HIGHENE spectrumcontrol=1 pindex=$pnpindex pnorm=$pnpnorm clobber=1 >& log/20sp_partial_pn.log


rot-im-det-sky prefix=$mos1input mask=1 elow=$LOWENE ehigh=$HIGHENE mode=2 clobber=1
rot-im-det-sky prefix=$mos2input mask=1 elow=$LOWENE ehigh=$HIGHENE mode=2 clobber=1
rot-im-det-sky prefix=$pninput mask=1 elow=$LOWENE ehigh=$HIGHENE mode=2 clobber=1


echo ""
echo "========= INPUT Parameters ========="
more $inputlog


#------------------------------------------------------------------------
# check files
if (! -f ${mos1prefix}-prot-im-sky-${LOWENE}-${HIGHENE}.fits) then
echo "ERROR : no output file "${mos1prefix}-prot-im-sky-${LOWENE}-${HIGHENE}.fits
exit 
endif
if (! -f ${mos2prefix}-prot-im-sky-${LOWENE}-${HIGHENE}.fits) then
echo "ERROR : no output file "${mos2prefix}-prot-im-sky-${LOWENE}-${HIGHENE}.fits
exit
endif
if (! -f ${pnprefix}-prot-im-sky-${LOWENE}-${HIGHENE}.fits) then
echo "ERROR : no output file "${pnprefix}-prot-im-sky-${LOWENE}-${HIGHENE}.fits
exit
endif
#-----------------------------------------------------------------------


echo "======================================================="
echo "if 6 parameters are found, you have successed."
echo "Next run : "
echo "bin/remake_image.csh"
echo "======================================================="

