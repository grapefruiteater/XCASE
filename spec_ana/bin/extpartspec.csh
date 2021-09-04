#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`
#======================================================================================================
if ($#argv != 3) then
    echo "#----------------usage----------------#"
    echo "extpartspec.csh RA Dec r"
    echo ""
    echo "RA  : selected center [ deg ]"
    echo "Dec : selected center [ deg ]"
    echo "r : radius [ arcsec ]"
    echo ""
    echo "Example :"
    echo "bin/extpartspec.csh 120.2 -1.12 60  "
    echo "#-------------------------------------#"
exit
endif

set RA=$1
set Dec=$2
set secradius=$3
set degradius=`echo "scale=8; $3 / 60 / 60" |bc `

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

set MOS1clean = ${mos1prefix}-clean.fits
set MOS2clean = ${mos2prefix}-clean.fits
set PNclean = ${pnprefix}-clean.fits
set PNcleanOOT = ${pnprefix}-clean-oot.fits

#set DetX = `bin/esky2det2.csh $RA $Dec ${mos1prefix}-obj-image-sky.fits | gawk '{if(NF==2) print $1}'`
#set DetY = `bin/esky2det2.csh $RA $Dec ${mos1prefix}-obj-image-sky.fits | gawk '{if(NF==2) print $2}'`
#echo $DetX
#echo $DetY

evselect table=${MOS1clean} imagebinning=binSize imageset=MOS1image.fits withimageset=yes \
   xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80
evselect table=${MOS2clean} imagebinning=binSize imageset=MOS2image.fits withimageset=yes \
   xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80
evselect table=$PNclean imagebinning=binSize imageset=PNimage.fits withimageset=yes \
   xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80
evselect table=$PNcleanOOT imagebinning=binSize imageset=PNimage-oot.fits withimageset=yes \
   xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80

rm PNimage-oot_rescaled.fits PNclean_image.fits
set ootscale=`echo "scale=8; 0.063*1.0" |bc `
farith PNimage-oot.fits $ootscale PNimage-oot_rescaled.fits MUL
farith PNimage.fits PNimage-oot_rescaled.fits PNclean_image.fits SUB

#MOS1
ds9 -xpa localhost &
sleep 10
cat MOS1image.fits | xpaset ds9 fits
xpaset -p ds9 scale log
os.system("xpaset -p ds9 cmap b")
xpaset -p ds9 regions system physical
echo "fk5; circle $RA $Dec $degradius" | xpaset ds9 regions
rm tmp.reg
xpaset -p ds9 regions save ./tmp.reg -format ds9 -system wcs -physical
sleep 3
xpaset -p ds9 exit

set X=` tail -1 tmp.reg | gawk '{print $1}'`
set arr=( `echo $X | tr -s ',' ' '| tr -s '(' ' '| tr -s ')' ' '`)

set X=${arr[2]}
set Y=${arr[3]}
set R1=${arr[4]}
set R2=`echo "scale=8; $R1 + 1200" |bc `

rm MOS1source_spectrum.fits MOS1background_spectrum.fits MOS1.rmf MOS1.arf MOS1_spectrum_grp.fits
evselect table=$MOS1clean withspectrumset=yes spectrumset=MOS1source_spectrum.fits \
   energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999 \
   expression="#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN circle($X,$Y,$R1))"
evselect table=$MOS1clean withspectrumset=yes spectrumset=MOS1background_spectrum.fits \
   energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999 \
   expression="#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN annulus($X,$Y,$R1,$R2))"

#MOS2
ds9 -xpa localhost &
sleep 10
cat MOS2image.fits | xpaset ds9 fits
xpaset -p ds9 scale log
xpaset -p ds9 regions system physical
echo "fk5; circle $RA $Dec $degradius" | xpaset ds9 regions
rm tmp.reg
xpaset -p ds9 regions save ./tmp.reg -format ds9 -system wcs -physical
sleep 3
xpaset -p ds9 exit

set X=` tail -1 tmp.reg | gawk '{print $1}'`
set arr=( `echo $X | tr -s ',' ' '| tr -s '(' ' '| tr -s ')' ' '`)

set X=${arr[2]}
set Y=${arr[3]}
set R1=${arr[4]}
set R2=`echo "scale=8; $R1 + 1200" |bc `

rm MOS2source_spectrum.fits MOS2background_spectrum.fits MOS2.rmf MOS2.arf MOS2_spectrum_grp.fits
evselect table=$MOS2clean withspectrumset=yes spectrumset=MOS2source_spectrum.fits \
   energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999 \
   expression="#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN circle($X,$Y,$R1))"
evselect table=$MOS2clean withspectrumset=yes spectrumset=MOS2background_spectrum.fits \
   energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999 \
   expression="#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN annulus($X,$Y,$R1,$R2))"

#PN
ds9 -xpa localhost &
sleep 10
cat PNimage.fits | xpaset ds9 fits
xpaset -p ds9 scale log
xpaset -p ds9 regions system physical
echo "fk5; circle $RA $Dec $degradius" | xpaset ds9 regions
rm tmp.reg
xpaset -p ds9 regions save ./tmp.reg -format ds9 -system wcs -physical
sleep 3
xpaset -p ds9 exit

set X=` tail -1 tmp.reg | gawk '{print $1}'`
set arr=( `echo $X | tr -s ',' ' '| tr -s '(' ' '| tr -s ')' ' '`)

set X=${arr[2]}
set Y=${arr[3]}
set R1=${arr[4]}
set R2=`echo "scale=8; $R1 + 1200" |bc `

rm PNsource_spectrum.fits PNbackground_spectrum.fits PN.rmf PN.arf PN_spectrum_grp.fits
evselect table=$PNclean withspectrumset=yes spectrumset=PNsource_spectrum.fits \
   energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 \
   expression="(FLAG==0) && (PATTERN<=4) && ((X,Y) IN circle($X,$Y,$R1))"
rm PN_oot_spectrum.fits PNbackground_spectrum-oot.fits
evselect table=$PNcleanOOT withspectrumset=yes spectrumset=PN_oot_spectrum.fits \
   energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 \
   expression="(FLAG==0) && (PATTERN<=4) && ((X,Y) IN circle($X,$Y,$R1))"
fparkey value=CTS_OOT fitsfile=PN_oot_spectrum.fits keyword=TTYPE2
faddcol infile=PNsource_spectrum.fits colfile=PN_oot_spectrum.fits colname=CTS_OOT
fcalc clobber=yes infile=PNsource_spectrum.fits outfile=PNsource_spectrum.fits \
  clname=CTS_OOT expr=CTS_OOT*0.063
fcalc clobber=yes infile=PNsource_spectrum.fits outfile=PNsource_spectrum.fits \
  clname=COUNTS expr=COUNTS-CTS_OOT

evselect table=$PNclean withspectrumset=yes spectrumset=PNbackground_spectrum.fits \
   energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 \
   expression="(FLAG==0) &&(PATTERN<=4) && ((X,Y) IN annulus($X,$Y,$R1,$R2))"
evselect table=$PNcleanOOT withspectrumset=yes spectrumset=PNbackground_spectrum-oot.fits \
   energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 \
   expression="(FLAG==0) && (PATTERN<=4) && ((X,Y) IN annulus($X,$Y,$R1,$R2))"
fparkey value=CTS_OOT fitsfile=PNbackground_spectrum-oot.fits keyword=TTYPE2
faddcol infile=PNbackground_spectrum.fits colfile=PNbackground_spectrum-oot.fits colname=CTS_OOT
fcalc clobber=yes infile=PNbackground_spectrum.fits outfile=PNbackground_spectrum.fits \
  clname=CTS_OOT expr=CTS_OOT*0.063
fcalc clobber=yes infile=PNbackground_spectrum.fits outfile=PNbackground_spectrum.fits \
  clname=COUNTS expr=COUNTS-CTS_OOT

backscale spectrumset=MOS1source_spectrum.fits badpixlocation=$MOS1clean
backscale spectrumset=MOS1background_spectrum.fits badpixlocation=$MOS1clean
backscale spectrumset=MOS2source_spectrum.fits badpixlocation=$MOS2clean
backscale spectrumset=MOS2background_spectrum.fits badpixlocation=$MOS2clean
backscale spectrumset=PNsource_spectrum.fits badpixlocation=$PNclean
backscale spectrumset=PNbackground_spectrum.fits badpixlocation=$PNclean

rmfgen spectrumset=MOS1source_spectrum.fits rmfset=MOS1.rmf
rmfgen spectrumset=MOS2source_spectrum.fits rmfset=MOS2.rmf
rmfgen spectrumset=PNsource_spectrum.fits rmfset=PN.rmf

arfgen spectrumset=MOS1source_spectrum.fits arfset=MOS1.arf withrmfset=yes rmfset=MOS1.rmf \
   badpixlocation=$MOS1clean detmaptype=psf
arfgen spectrumset=MOS2source_spectrum.fits arfset=MOS2.arf withrmfset=yes rmfset=MOS2.rmf \
   badpixlocation=$MOS2clean detmaptype=psf
arfgen spectrumset=PNsource_spectrum.fits arfset=PN.arf withrmfset=yes rmfset=PN.rmf \
   badpixlocation=$PNclean detmaptype=psf

grppha MOS1source_spectrum.fits MOS1_spectrum_grp.fits <<EOF
chkey BACKFILE MOS1background_spectrum.fits
chkey RESPFILE MOS1.rmf
chkey ANCRFILE MOS1.arf
group min 35
exit
EOF

grppha MOS2source_spectrum.fits MOS2_spectrum_grp.fits <<EOF
chkey BACKFILE MOS2background_spectrum.fits
chkey RESPFILE MOS2.rmf
chkey ANCRFILE MOS2.arf
group min 35
exit
EOF

grppha PNsource_spectrum.fits PN_spectrum_grp.fits <<EOF
chkey BACKFILE PNbackground_spectrum.fits
chkey RESPFILE PN.rmf
chkey ANCRFILE PN.arf
group min 35
exit
EOF


rm pointspec_apec_mos1.ps pointspec_apec_mos1.qdp pointspec_apec_mos1.pco save-pointspec_apec_mos1.xcm
cat <<EOF|xspec 
log log/xspec.tbl
@pointspec_apec_mos1.xcm
log none
pl ldata del
iplot

wen pointspec_apec_mos1
h pointspec_apec_mos1.ps/cps
exit
exit
EOF

rm pointspec_powerlaw_mos1.ps pointspec_powerlaw_mos1.qdp pointspec_powerlaw_mos1.pco save-pointspec_powerlaw_mos1.xcm
cat <<EOF|xspec 
log log/xspec.tbl
@pointspec_powerlaw_mos1.xcm
log none
pl ldata del
iplot

wen pointspec_powerlaw_mos1
h pointspec_powerlaw_mos1.ps/cps
exit
exit
EOF

rm pointspec_apec_mos2.ps pointspec_apec_mos2.qdp pointspec_ape_mos2c.pco save-pointspec_apec_mos2.xcm
cat <<EOF|xspec 
log log/xspec.tbl
@pointspec_apec_mos2.xcm
log none
pl ldata del
iplot

wen pointspec_apec_mos2
h pointspec_apec_mos2.ps/cps
exit
exit
EOF
rm pointspec_powerlaw_mos2.ps pointspec_powerlaw_mos2.qdp pointspec_powerlaw_mos2.pco save-pointspec_powerlaw_mos2.xcm
cat <<EOF|xspec 
log log/xspec.tbl
@pointspec_powerlaw_mos2.xcm
log none
pl ldata del
iplot

wen pointspec_powerlaw_mos2
h pointspec_powerlaw_mos2.ps/cps
exit
exit
EOF

rm pointspec_apec_pn.ps pointspec_apec_pn.qdp pointspec_apec_pn.pco save-pointspec_apec_pn.xcm
cat <<EOF|xspec 
log log/xspec.tbl
@pointspec_apec_pn.xcm
log none
pl ldata del
iplot

wen pointspec_apec
h pointspec_apec_pn.ps/cps
exit
exit
EOF

rm pointspec_powerlaw_pn.ps pointspec_powerlaw_pn.qdp pointspec_powerlaw_pn.pco save-pointspec_powerlaw_pn.xcm
cat <<EOF|xspec 
log log/xspec.tbl
@pointspec_powerlaw_pn.xcm
log none
pl ldata del
iplot

wen pointspec_powerlaw_pn
h pointspec_powerlaw_pn.ps/cps
exit
exit
EOF
