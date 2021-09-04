#!/bin/csh -f

set prefix=$1
set region_name=$2

rm ../data/${prefix}-obj-${region_name}-grp_noback.pi
grppha ../data/${prefix}-obj-${region_name}.pi ../data/${prefix}-obj-${region_name}-grp_noback.pi <<EOF
chkey RESPFILE ../data/${prefix}-${region_name}.rmf 
chkey ANCRFILE ../data/${prefix}-${region_name}.arf 
group min 25
show group
exit
EOF

