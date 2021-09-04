#!/bin/csh -f

source ${SAS_DIR}/setsas.csh >& /dev/null
setenv SAS_CCF ccf.cif 
setenv SAS_ODF `ls -1 *SUM.SAS`
 
#protonコマンドの出力 ”Scaled Normalization” を入れる。
set mos1pnorm = 
set mos2pnorm = 
set pnpnorm   =
#スペクトルフィットから得られるベストフィットパラメータを入れる。
set mos1pindex = 
set mos2pindex = 
set pnpindex  = 

proton prefix=1S001 caldb=$SAS_ESAS_CALDB specname=mos1S001-obj-full-grp.pi ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 elow=$LOWENE ehigh=$HIGHENE spectrumcontrol=1 pindex=$mos1pindex pnorm=$mos1pnorm clobber=1 >& log/20sp_partial_mos1.log
proton prefix=2S002 caldb=$SAS_ESAS_CALDB specname=mos2S002-obj-full-grp.pi ccd1=1 ccd2=1 ccd3=1 ccd4=1 ccd5=1 ccd6=1 ccd7=1 elow=$LOWENE ehigh=$HIGHENE spectrumcontrol=1 pindex=$mos2pindex pnorm=$mos2pnorm clobber=1 >& log/20sp_partial_mos2.log
proton prefix=S003 caldb=$SAS_ESAS_CALDB specname=pnS003-obj-full-grp.pi  elow=$LOWENE ehigh=$HIGHENE spectrumcontrol=1 pindex=$pnpindex pnorm=$pnpnorm clobber=1 >& log/20sp_partial_pn.log

rot-im-det-sky prefix=1S001 mask=1 elow=$LOWENE ehigh=$HIGHENE mode=2
rot-im-det-sky prefix=2S002 mask=1 elow=$LOWENE ehigh=$HIGHENE mode=2
rot-im-det-sky prefix=S003 mask=1 elow=$LOWENE ehigh=$HIGHENE mode=2

###---make image---

#MOS1&MOS2&PN
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="1S001 2S002 S003" >& log/21comb.log
bin_image thresholdmasking=0.02 detector=0 prefix="1S001 2S002 S003" binning=1 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=no >& log/22bin_image.log
adapt smoothingcounts=50 thresholdmasking=0.02 detector=0 binning=2 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=0 >& log/23adapt.log
mv adapt-$LOWENE-$HIGHENE.fits adapt-$LOWENE-$HIGHENE-all-sp.fits
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-all-sp.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-all-sp.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-all-sp.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-all-sp.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-all-sp.fits
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-all-sp.fits
#MOS1&MOS2
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="1S001 2S002" >& log/21comb_mos.log
bin_image thresholdmasking=0.02 detector=0 prefix="1S001 2S002" binning=1 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=no >& log/23bin_image_mos.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos-sp.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos-sp.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos-sp.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos-sp.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos-sp.fits
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos-sp.fits
#MOS1
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="1S001" >& log/21comb_mos1.log
bin_image thresholdmasking=0.02 detector=1 prefix="1S001" binning=1 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=no >& log/23bin_image_mos1.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos1-sp.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos1-sp.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos1-sp.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos1-sp.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos1-sp.fits
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos1-sp.fits
#MOS2
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="2S002" >& log/21comb_mos2.log
bin_image thresholdmasking=0.02 detector=1 prefix="2S002" binning=1 elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=no >& log/23bin_image_mos2.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-mos2-sp.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-mos2-sp.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-mos2-sp.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-mos2-sp.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-mos2-sp.fits
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-mos2-sp.fits
#PN
comb caldb=$SAS_ESAS_CALDB withpartcontrol=1 withsoftcontrol=1 withswcxcontrol=0 elowlist=$LOWENE ehighlist=$HIGHENE mask=1 prefixlist="S003" >& log/21comb_pn.log
bin_image thresholdmasking=0.02 detector=2 binning=1 prefix="S003" elow=$LOWENE ehigh=$HIGHENE withpartcontrol=yes withsoftcontrol=yes withswcxcontrol=no >& log/23bin_image_pn.log
mv rate-$LOWENE-$HIGHENE.fits rate-$LOWENE-$HIGHENE-pn-sp.fits
mv sigma-$LOWENE-$HIGHENE.fits sigma-$LOWENE-$HIGHENE-pn-sp.fits
mv comb-obj-im-$LOWENE-$HIGHENE.fits comb-obj-im-$LOWENE-$HIGHENE-pn-sp.fits
mv comb-back-im-sky-$LOWENE-$HIGHENE.fits comb-back-im-sky-$LOWENE-$HIGHENE-pn-sp.fits
mv comb-exp-im-$LOWENE-$HIGHENE.fits comb-exp-im-$LOWENE-$HIGHENE-pn-sp.fits
mv comb-prot-im-sky-$LOWENE-$HIGHENE.fits comb-prot-im-sky-$LOWENE-$HIGHENE-pn-sp.fits

#Remask

farith rate-$LOWENE-$HIGHENE-mos1-sp.fits mos1S001-cheese.fits rate-$LOWENE-$HIGHENE-mos1-sp.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-mos1-sp.fits mos1S001-cheese.fits sigma-$LOWENE-$HIGHENE-mos1-sp.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-mos2-sp.fits mos2S002-cheese.fits rate-$LOWENE-$HIGHENE-mos2-sp.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-mos2-sp.fits mos2S002-cheese.fits sigma-$LOWENE-$HIGHENE-mos2-sp.fits.mask MUL
farith rate-$LOWENE-$HIGHENE-pn-sp.fits pnS003-cheese.fits rate-$LOWENE-$HIGHENE-pn-sp.fits.mask MUL
farith sigma-$LOWENE-$HIGHENE-pn-sp.fits pnS003-cheese.fits sigma-$LOWENE-$HIGHENE-pn-sp.fits.mask MUL


##---------------------------------------------------------------------------

#DET座標の銀河団中心位置を入れる。

set CentX_MOS1 = 
set CentY_MOS1 = 

set CentX_MOS2 = 
set CentY_MOS2 = 

set CentX_PN   =
set CentY_PN   =
 
rm mos1_psf.fits mos2_psf.fits pn_psf.fits
psfgen region="(DETX, DETY) IN box(CentX_MOS1,CentY_MOS1,2000,2000)" instrument=M1 output=mos1_psf.fits level=ELLBETA energy=1400 coortype=XY 
psfgen region="(DETX, DETY) IN box(CentX_MOS2,CentY_MOS2,2000,2000)" instrument=M2 output=mos2_psf.fits level=ELLBETA energy=1400 coortype=XY
psfgen region="(DETX, DETY) IN box(CentX_PN,CentY_PN,2000,2000)" instrument=PN output=pn_psf.fits level=ELLBETA energy=1400 coortype=XY

rm mos1_psf_2bin.fits mos2_psf_2bin.fits pn_psf_2bin.fits 
fimgbin mos1_psf.fits mos1_psf_2bin.fits 2
fimgbin mos2_psf.fits mos2_psf_2bin.fits 2
fimgbin pn_psf.fits pn_psf_2bin.fits 2


