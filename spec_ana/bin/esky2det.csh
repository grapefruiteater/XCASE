#!/bin/csh -f

set RA  = $argv[1]
set DEC = $argv[2]
#set DET = $argv[3]
set IMAGEFILE = $argv[3]

#source ${SAS_DIR}/setsas.csh >& /dev/null
#setenv SAS_CCF ccf.cif
#setenv SAS_ODF `ls -1 *SUM.SAS`

esky2det datastyle=user ra=${RA} dec=${DEC} checkfov=no outunit=det withheader=no calinfostyle=set calinfoset=${IMAGEFILE} verbosity=0 | head -1
#instrument=${DET} 
