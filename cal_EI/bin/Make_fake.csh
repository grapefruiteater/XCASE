#!/bin/csh -f

if ( $#argv != 1 ) then
cat <<EOF

bin/Make_fake.csh redshift

   redsfhit - cluster redshift

Use:  bin/Make_fake.csh 0.1832


EOF

exit

endif

set z=$1

#-----------------------------------------------------------------------------------------
#setup prefix
#e.g. mos1prefix=mos1S001 mos1input=1S001 
set prefixlog="../spec_ana/log/prefix.log"
set mos1prefix=`head -1 $prefixlog | gawk '{print $1}'`
set mos2prefix=`head -2 $prefixlog | tail -1 |  gawk '{print $1}'`
set pnprefix=`tail -1 $prefixlog | gawk '{print $1}'`

set mos1input=`head -1 $prefixlog | gawk '{print $2}'`
set mos2input=`head -2 $prefixlog | tail -1 |  gawk '{print $2}'`
set pninput=`tail -1 $prefixlog | gawk '{print $2}'`
set mosinput=`echo "'"$mos1input" "$mos2input"'"`
#-----------------------------------------------------------------------------------------
# best fit file : e.g. save-all.xcm
set xspecoutput=`gawk '{if($1=="Best") print "../spec_ana/"$3}' ../spec_ana/log/sp_partial-inputfile.log`
# err fit file : e.g. save-err.xcm
set xspecoutputErr=`gawk '{if($1=="Err") print "../spec_ana/"$3}' ../spec_ana/log/sp_partial-inputfile.log`
#-----------------------------------------------------------------------------------------
set mos1qdp=fake-data_${mos1prefix}.qdp
set mos2qdp=fake-data_${mos2prefix}.qdp
set pnqdp=fake-data_${pnprefix}.qdp
#------------------------------------------------------------------------------------------
echo "Making Xspec scripts..."
bin/make_fake_xcm.py $mos1prefix $z $xspecoutput
bin/make_fake_xcm.py $mos2prefix $z $xspecoutput
bin/make_fake_xcm.py $pnprefix $z $xspecoutput

#---------------------------------------------------------------------------
echo "Making fake spectrum..."
foreach qdpfile ( $mos1qdp $mos2qdp $pnqdp )
if ( -f $qdpfile) then
rm $qdpfile
endif
end
rm *fake*.pi
csh fake_step1_${mos1prefix}.csh
csh fake_step2_${mos1prefix}.csh
csh fake_step1_${mos2prefix}.csh
csh fake_step2_${mos2prefix}.csh
csh fake_step1_${pnprefix}.csh
csh fake_step2_${pnprefix}.csh
#---------------------------------------------------------------------------
# check fake pi file
set regfile="../spec_ana/SetUpFile/Region.dat"
set Nreg=`more $regfile | wc -l | gawk '{print $1-2}'`
foreach prefix ( $mos1prefix $mos2prefix $pnprefix )
set i=1
while ( $i <= $Nreg )
set i2=`echo $i | gawk '{print $1+1}'`
set rinn=`head -$i $regfile | tail -1 | gawk '{print $1}'`
set rout=`head -$i2 $regfile | tail -1 | gawk '{print $1}'`
set region=${rinn}-${rout}
set file=${prefix}fake-${region}.pi
#echo $file
if ( -z $file ) then
   echo "ERROR : no file "$file
exit 
endif
  @ i ++
end
end


#---------------------------------------------------------------------------
#cautions!!
#Sometime fails interpolation
echo "Calculate a conversion factor from Sx to ne...."
foreach qdpfile ( $mos1qdp $mos2qdp $pnqdp )
if (! -f $qdpfile) then
echo "ERROR : no file of "$qdpfile
exit 
endif
echo "NO NO NO NO" >> $qdpfile
end

bin/lambda_cal.py $mos1prefix $z $mos1qdp 
bin/lambda_cal.py $mos2prefix $z $mos2qdp 
bin/lambda_cal.py $pnprefix $z $pnqdp 

#---------------------------------------------------------------------------
bin/temp_cal.py $z $xspecoutput $xspecoutputErr

#--------------------------------------------------------------------------
# Check whether output files exist or not
#-------------------------------------------------------------------------
foreach prefix ( $mos1prefix $mos2prefix $pnprefix )
set output=../result/Lambda_${prefix}.dat
if (! -f $output) then
echo "ERROR : no file "$output
exit
endif
if ( -z $output) then
echo "ERROR : "$output" exists but zero size"
exit
endif

end

echo "======================================================="
echo "Success  : all files are completed"
echo "Next run : "
echo "cd ../radial_fit"
echo "bin/fittting.csh cogname z rinn rout Nbin"
echo "======================================================="
