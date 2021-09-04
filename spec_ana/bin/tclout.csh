#!/bin/csh -f

rm xspec_fit_result.dat
cat <<EOF|xspec 
@save-all_cross.xcm
fit

log xspec_fit_result.dat
show free

exit
exit
EOF
