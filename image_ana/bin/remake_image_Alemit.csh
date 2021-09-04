#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`

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
set LOWENE="1400"
set HIGHENE="1600"

#-----------------------------------------------------------------------------------------
###---make image---

#MOS1&MOS2&PN
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="$mos1input $mos2input $pninput" >& log/21comb.log
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-all.fits
#MOS1&MOS2
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="$mos1input $mos2input" >& log/21comb_mos.log
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos.fits
#MOS1
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="$mos1input" >& log/21comb_mos1.log
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos1.fits
#MOS2
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="$mos2input" >& log/21comb_mos2.log
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos2.fits
#PN
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="$pninput" >& log/21comb_pn.log
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-pn.fits

#resetup energy range  
set eRangeFile="../spec_ana/SetUpFile/eRange.dat"
set LOWENE=`head -1 $eRangeFile | gawk '{print $1}'`
set HIGHENE=`head -2 $eRangeFile | gawk 'NR==2 {print $1}'`

rm -f tmp1.fits tmp2.fits tmp3.fits rate-$LOWENE-$HIGHENE-mos1-sp_self_v1.fits rate-$LOWENE-$HIGHENE-mos1-sp_self_v2.fits rate-$LOWENE-$HIGHENE-mos1-sp_selfscal_v1.fits rate-$LOWENE-$HIGHENE-mos1-sp_selfscal_v2.fits rate-$LOWENE-$HIGHENE-mos1-sp_selfscal_v1.fits.mask rate-$LOWENE-$HIGHENE-mos1-sp_selfscal_v2.fits.mask
farith comb-obj-im-$LOWENE-$HIGHENE-mos1.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos1-sp.fits tmp1.fits SUB 
farith tmp1.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos1-sp.fits tmp2.fits SUB 
farith tmp2.fits comb-exp-im-$LOWENE-$HIGHENE-mos1.fits rate-$LOWENE-$HIGHENE-mos1-sp_self_v1.fits DIV
farith tmp2.fits mos${mos1input}-exp-im-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos1-sp_self_v2.fits DIV
fcarith tmp2.fits 2073600.0 tmp3.fits MUL
farith tmp3.fits comb-exp-im-$LOWENE-$HIGHENE-mos1.fits rate-$LOWENE-$HIGHENE-mos1-sp_selfscal_v1.fits DIV
farith tmp3.fits mos${mos1input}-exp-im-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos1-sp_selfscal_v2.fits DIV
farith rate-$LOWENE-$HIGHENE-mos1-sp_selfscal_v1.fits ${mos1prefix}-cheese.fits rate-$LOWENE-$HIGHENE-mos1-sp_selfscal_v1.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-mos1-sp_selfscal_v2.fits ${mos1prefix}-cheese.fits rate-$LOWENE-$HIGHENE-mos1-sp_selfscal_v2.fits.mask MUL

rm -f tmp1.fits tmp2.fits tmp3.fits rate-$LOWENE-$HIGHENE-mos2-sp_self_v1.fits rate-$LOWENE-$HIGHENE-mos2-sp_self_v2.fits rate-$LOWENE-$HIGHENE-mos2-sp_selfscal_v1.fits rate-$LOWENE-$HIGHENE-mos2-sp_selfscal_v2.fits rate-$LOWENE-$HIGHENE-mos2-sp_selfscal_v1.fits.mask rate-$LOWENE-$HIGHENE-mos2-sp_selfscal_v2.fits.mask
farith comb-obj-im-$LOWENE-$HIGHENE-mos2.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos2-sp.fits tmp1.fits SUB
farith tmp1.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos2-sp.fits tmp2.fits SUB
farith tmp2.fits comb-exp-im-$LOWENE-$HIGHENE-mos2.fits rate-$LOWENE-$HIGHENE-mos2-sp_self_v1.fits DIV
farith tmp2.fits mos${mos2input}-exp-im-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos2-sp_self_v2.fits DIV
fcarith tmp2.fits 2073600.0 tmp3.fits MUL
farith tmp3.fits comb-exp-im-$LOWENE-$HIGHENE-mos2.fits rate-$LOWENE-$HIGHENE-mos2-sp_selfscal_v1.fits DIV
farith tmp3.fits mos${mos2input}-exp-im-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos2-sp_selfscal_v2.fits DIV
farith rate-$LOWENE-$HIGHENE-mos2-sp_selfscal_v1.fits ${mos2prefix}-cheese.fits rate-$LOWENE-$HIGHENE-mos2-sp_selfscal_v1.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-mos2-sp_selfscal_v2.fits ${mos2prefix}-cheese.fits rate-$LOWENE-$HIGHENE-mos2-sp_selfscal_v2.fits.mask MUL

rm -f tmp1.fits tmp2.fits tmp3.fits rate-$LOWENE-$HIGHENE-pn-sp_self_v1.fits rate-$LOWENE-$HIGHENE-pn-sp_self_v2.fits rate-$LOWENE-$HIGHENE-pn-sp_selfscal_v1.fits rate-$LOWENE-$HIGHENE-pn-sp_selfscal_v2.fits rate-$LOWENE-$HIGHENE-pn-sp_selfscal_v1.fits.mask rate-$LOWENE-$HIGHENE-pn-sp_selfscal_v2.fits.mask
farith comb-obj-im-$LOWENE-$HIGHENE-pn.fits comb-back-im-sky-$LOWENE-$HIGHENE-pn-sp.fits tmp1.fits SUB
farith tmp1.fits comb-prot-im-sky-$LOWENE-$HIGHENE-pn-sp.fits tmp2.fits SUB
farith tmp2.fits comb-exp-im-$LOWENE-$HIGHENE-pn.fits rate-$LOWENE-$HIGHENE-pn-sp_self_v1.fits DIV
farith tmp2.fits pn${pninput}-exp-im-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-pn-sp_self_v2.fits DIV
fcarith tmp2.fits 2073600.0 tmp3.fits MUL
farith tmp3.fits comb-exp-im-$LOWENE-$HIGHENE-pn.fits rate-$LOWENE-$HIGHENE-pn-sp_selfscal_v1.fits DIV
farith tmp3.fits pn${pninput}-exp-im-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-pn-sp_selfscal_v2.fits DIV
farith rate-$LOWENE-$HIGHENE-pn-sp_selfscal_v1.fits ${pnprefix}-cheese.fits rate-$LOWENE-$HIGHENE-pn-sp_selfscal_v1.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-pn-sp_selfscal_v2.fits ${pnprefix}-cheese.fits rate-$LOWENE-$HIGHENE-pn-sp_selfscal_v2.fits.mask MUL

rm -f tmp-obj.fits tmp-back.fits tmp-prot.fits tmp-obj_sub.fits tmp-obj_sub_sub.fits tmp-obj_sub_sub_scal.fits rate-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits tmp-sigma.fits sigma-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits
farith comb-obj-im-$LOWENE-$HIGHENE-mos1.fits comb-obj-im-1400-1600-mos1.fits tmp-obj.fits SUB 
farith comb-back-im-sky-$LOWENE-$HIGHENE-mos1-sp.fits comb-back-im-sky-1400-1600-mos1.fits tmp-back.fits SUB 
farith comb-prot-im-sky-$LOWENE-$HIGHENE-mos1-sp.fits comb-prot-im-sky-1400-1600-mos1.fits tmp-prot.fits SUB 
farith tmp-obj.fits tmp-back.fits tmp-obj_sub.fits SUB 
farith tmp-obj_sub.fits tmp-prot.fits tmp-obj_sub_sub.fits SUB 
fcarith tmp-obj_sub_sub.fits 2073600.0 tmp-obj_sub_sub_scal.fits MUL
ftimgcalc tmp-sigma.fits "(DATA)**0.5*60*60*60*60/2.5/2.5" a="DATA=tmp-obj.fits"
farith tmp-obj_sub_sub_scal.fits mos${mos1input}-exp-im-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits DIV
farith tmp-sigma.fits mos${mos1input}-exp-im-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits DIV
#ftimgcalc sigma-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits "DATA!=-999?DATA/EXP:0.0" a="DATA=tmp-sigma.fits" b="EXP=mos${mos1input}-exp-im-$LOWENE-$HIGHENE.fits"

rm -f tmp-obj.fits tmp-back.fits tmp-prot.fits tmp-obj_sub.fits tmp-obj_sub_sub.fits tmp-obj_sub_sub_scal.fits rate-$LOWENE-$HIGHENE-mos2-sp_Alemit.fits tmp-sigma.fits sigma-$LOWENE-$HIGHENE-mos2-sp_Alemit.fits
farith comb-obj-im-$LOWENE-$HIGHENE-mos2.fits comb-obj-im-1400-1600-mos2.fits tmp-obj.fits SUB
farith comb-back-im-sky-$LOWENE-$HIGHENE-mos2-sp.fits comb-back-im-sky-1400-1600-mos2.fits tmp-back.fits SUB
farith comb-prot-im-sky-$LOWENE-$HIGHENE-mos2-sp.fits comb-prot-im-sky-1400-1600-mos2.fits tmp-prot.fits SUB
farith tmp-obj.fits tmp-back.fits tmp-obj_sub.fits SUB
farith tmp-obj_sub.fits tmp-prot.fits tmp-obj_sub_sub.fits SUB
fcarith tmp-obj_sub_sub.fits 2073600.0 tmp-obj_sub_sub_scal.fits MUL
ftimgcalc tmp-sigma.fits "DATA**0.5*60*60*60*60/2.5/2.5" a="DATA=tmp-obj.fits"
#ftimgcalc tmp-sigma.fits "(DATA-BACK)**0.5*60*60*60*60/2.5/2.5" a="DATA=tmp-obj.fits" b="BACK=tmp-back.fits"
farith tmp-obj_sub_sub_scal.fits mos${mos2input}-exp-im-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos2-sp_Alemit.fits DIV
farith tmp-sigma.fits mos${mos2input}-exp-im-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos2-sp_Alemit.fits DIV

rm -f tmp-obj.fits tmp-back.fits tmp-prot.fits tmp-obj_sub.fits tmp-obj_sub_sub.fits tmp-obj_sub_sub_scal.fits rate-$LOWENE-$HIGHENE-pn-sp_Alemit.fits tmp-sigma.fits sigma-$LOWENE-$HIGHENE-pn-sp_Alemit.fits
farith comb-obj-im-$LOWENE-$HIGHENE-pn.fits comb-obj-im-1400-1600-pn.fits tmp-obj.fits SUB
farith comb-back-im-sky-$LOWENE-$HIGHENE-pn-sp.fits comb-back-im-sky-1400-1600-pn.fits tmp-back.fits SUB
farith comb-prot-im-sky-$LOWENE-$HIGHENE-pn-sp.fits comb-prot-im-sky-1400-1600-pn.fits tmp-prot.fits SUB
farith tmp-obj.fits tmp-back.fits tmp-obj_sub.fits SUB
farith tmp-obj_sub.fits tmp-prot.fits tmp-obj_sub_sub.fits SUB
fcarith tmp-obj_sub_sub.fits 2073600.0 tmp-obj_sub_sub_scal.fits MUL
ftimgcalc tmp-sigma.fits "(DATA**0.5)*60*60*60*60/2.5/2.5" a="DATA=tmp-obj.fits"
#ftimgcalc tmp-sigma.fits "(DATA-BACK)**0.5*60*60*60*60/2.5/2.5" a="DATA=tmp-obj.fits" b="BACK=tmp-back.fits"
farith tmp-obj_sub_sub_scal.fits pn${pninput}-exp-im-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-pn-sp_Alemit.fits DIV
farith tmp-sigma.fits pn${pninput}-exp-im-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-pn-sp_Alemit.fits DIV

#Remask
foreach prefix ( "mos1" "mos2" "pn" )
if ( -f rate-${LOWENE}-${HIGHENE}-${prefix}-sp_Alemit.fits.mask) then
rm rate-${LOWENE}-${HIGHENE}-${prefix}-sp_Alemit.fits.mask
endif
if ( -f sigma-${LOWENE}-${HIGHENE}-${prefix}-sp_Alemit.fits.mask) then
rm sigma-${LOWENE}-${HIGHENE}-${prefix}-sp_Alemit.fits.mask
endif
end

rm rate-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits.mask rate-$LOWENE-$HIGHENE-mos2-sp_Alemit.fits.mask rate-$LOWENE-$HIGHENE-pn-sp_Alemit.fits.mask
rm sigma-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits.mask sigma-$LOWENE-$HIGHENE-mos2-sp_Alemit.fits.mask sigma-$LOWENE-$HIGHENE-pn-sp_Alemit.fits.mask
farith rate-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits ${mos1prefix}-cheese.fits rate-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits ${mos1prefix}-cheese.fits sigma-$LOWENE-$HIGHENE-mos1-sp_Alemit.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-mos2-sp_Alemit.fits ${mos2prefix}-cheese.fits rate-$LOWENE-$HIGHENE-mos2-sp_Alemit.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-mos2-sp_Alemit.fits ${mos2prefix}-cheese.fits sigma-$LOWENE-$HIGHENE-mos2-sp_Alemit.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-pn-sp_Alemit.fits ${pnprefix}-cheese.fits rate-$LOWENE-$HIGHENE-pn-sp_Alemit.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-pn-sp_Alemit.fits ${pnprefix}-cheese.fits sigma-$LOWENE-$HIGHENE-pn-sp_Alemit.fits.mask MUL


#---------------------------------------------------------------------------------------------
foreach prefix ( "mos1" "mos2" "pn" )

if (! -f rate-${LOWENE}-${HIGHENE}-${prefix}-sp_Alemit.fits.mask) then
echo "ERROR : no file of "rate-${LOWENE}-${HIGHENE}-${prefix}-sp_Alemit.fits.mask
exit 
endif
if (! -f sigma-${LOWENE}-${HIGHENE}-${prefix}-sp_Alemit.fits.mask) then
echo "ERROR : no file of "sigma-${LOWENE}-${HIGHENE}-${prefix}-sp_Alemit.fits.mask
exit 
endif
end
#---------------------------------------------------------------------------------


#copy file
cp rate-**.fits.mask sigma-**.fits.mask ../data/

