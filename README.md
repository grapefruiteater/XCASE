### User Manual
#### Original Developer of XCASE: K. Miyaoka, N. Okabe and T. Kitaguchi
 
# 推奨環境
* HEASoft 6.19
* SAS xmmsas_20160201_1833
* SAS Caldb ccf 2017/8 ごろに取得したもの
* ESAS Caldb esas-caldb-sasv13.tar.gz
* ds9 Version 7.4
* python module 
 - pyfits 3.4
 - scipy 0.19.1
 - numpy 1.12.0
 - astropy 2.0.1
 - pywcs　1.12
 - matplotlib==2.0.0
* library for c program
 - cfitsio 5.3.39
 - gsl 
 - MINUIT2
 - Eigen

# 環境変数の設定
### csh  
* HEASoft  
 - setenv HEADAS /hoge/x86_64-unknown-linux-gnu-libc2.17  
 - source ${HEADAS}/headas-init.csh  
* SAS  
 - setenv SAS_DIR /hoge/xmmsas_20160201_1833  
 - setenv SAS_CCFPATH /hoge/ccf  
 - setenv SAS_ESAS_CALDB /hoge/esas

### bash  
* HEASoft  
 - export HEADAS=/hoge/x86_64-unknown-linux-gnu-libc2.17  
 - source ${HEADAS}/headas-init.sh  
* SAS 
 - export SAS_DIR=/hoge/xmmsas_20160201_1833  
 - export SAS_CCFPATH=/hoge/ccf  
 - export SAS_ESAS_CALDB=/hoge/esas  


# フィッティング手法等
* Galactic absorption : phabs, nH : LAB survey
* Abundance table : Anders E. & Grevesse N. Geochimica et Cosmochimica Acta 53, 197 (1989)
* XSPEC 12.9.0o
* AtomDB : 2.0.2
* redshift : MCXC catalog
* chi2 fitting, 1 sigma statistical error
* background spectrum : RASS(annulus of 1~2 deg)
* 1d radial fitting for surface brightness and projected temperature 

# XCASEの動作確認
* sample data (ObsID)
 - Abell 1795 (0097820101)
 - Abell 1689 (0093030101)
 - Abell 1835 (0551830201)  

URL : [Test data for XCASE](http://home.hiroshima-u.ac.jp/m161855/index.html)  
URL : [MIRROR @ okabe HP ](http://theo.phys.sci.hiroshima-u.ac.jp/~okabe/files/hsc-x/)


# ディレクトリ構造  


~/XCASE/spec_ana  
　　　　/image_ana  
　　　　/cal_EI   
　　　　/radial_fit  
　　　　/odf    
　　　　/result  
　　　　/data  

# スクリプト群


  - setup.sh: 環境変数の設定
  - filter.csh: イベントデータを作成、2シグマの範囲外のカウントレートを除去する（フレアカット）
  - cheese.csh: 点源除去のためのマスクリストを作成  
  - make_mask.csh: マスクイメージを作成
  - make_image.csh: 各 EPIC 検出器のイメージを作成 (検出器全面のスペクトルファイルも作成)
  - SearchCentroid.py: 500 kpcの半径内のカウントで重み付けした銀河団中心位置
       		     + CAMIRA cluster を中心にした 2 次モーメント、3 次モーメント、歪度
  - make_spec.csh: 指定した領域のスペクトル、レスポンスファイルの作成
  - makecross_arf_***_part*.csh: 指定した領域のcross arfの作成
  - make_fitcode.py: XSPECのスペクトルフィットコードの作成
  - sp_partial.csh: フレアの引き残りイメージ作成に必要なスケールファクターを計算  
  - proton.csh: フレアの引き残りイメージ作成
  - remake_image.csh: フレアの引き残りイメージを差し引いたレートイメージを作成
  - psfgen.csh: psfmap作成
  - make_fake_xcm.py: 正味の密度への変換係数の計算に用いるシュミレーションコードの作成
  - lambda_cal.py: 正味の密度への変換係数の計算
  - fitting.csh: 表面輝度分布、温度分布のフィッティングおよびガス質量や全質量の計算

# スクリプト実行手順

   0. odfの展開  
      `cd XCASE`  
      `$cp odffile.tar.gz odf/`  
      `cd odf`  
      `tar zxvf odffile.tar.gz`  
      `tar xvf ****_********.TAR `  
      `cd ../spec_ana`  

   1. データスクリーニング  
      `~/XCASE/spec_ana$bin/filter.csh`  
      イベントデータの作成、2シグマの範囲外のカウントレートの除去（フレアカット）を行う。
      
   2. 点源除去  
      `~/XCASE$echo "400\n2300" >& SetUpFile/eRange.dat`  
      spec_ana/SetUpFile/eRange.datに、イメージで用いるエネルギーの範囲を書き込む。  
      上記は400~2300 eVのイメージを作成する場合の例。  
      > 例)eRange.datの中身  
      >　 400  
      >　 2300
	    
      `~/XCASE/spec_ana$bin/cheese.csh scale rate dist`  			        
      点源除去のための点源を検出し、点源のリストを作成する。      
      * scale : 点源半径を調節するパラメータ。明るさがバックグラウンドのscale倍になる半径。   
      * rate  : 点源を検出するカウントレートの最小スレッショルド値。   
      * dist  : 点源間の距離の最小値(arcsec)。distより近い距離にある点源達は検出されない。
      
      > 例）~/XCASE/spec_ana$bin/cheese.csh 0.25 0.001 0
      
      `~/XCASE/spec_ana$ds9 mos1S001-obj-image-sky.fits &`   
      銀河団領域や検出器の溝など誤って点源と検出された点源をリージョンファイルから省く。  
	    ds9上でregion file "mos1S001-bkg_region-sky.reg"を開く。  
	    省きたい点源領域を削除し、fk5のunitでregion fileとして保存する。	     
 	   > 例）points.regという名前で保存した場合のその中身  
	       Region file format: DS9 version 4.1  
 	       global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1  
	       fk5  
	       circle(13:11:19.794,-1:21:49.73,20.7084")  
	       circle(13:11:20.488,-1:23:04.38,18.5953")  
	       
     `~/XCASE/spec_ana$bin/make_mask.csh points.reg`  
     点源リストを更新し、マスクイメージを作成する。  

   3. イメージ作成    
      `~/XCASE/spec_ana$bin/make_image.csh binning`  
      make_image.cshを実行し、各EPIC検出器のイメージを作成する。  
      * binnning : イメージのビンニング(1~4の整数)。  
      1は900×900、2は450×450、3は300×300、4は225×225 [pixel×pixel]のイメージを作成する。

   4. 中心決定  
      `~/XCASE/$emacs SetUpFile/Region.dat`     
      spec_ana/SetUpFile/Region.datに、銀河団中心から、スペクトルを抽出する同心円円環状の領域をarcsec unitで指定する。
      >例)Region.datの中身  
      >　 0  
      >　 20  
      >　 40  
      >　 60  
      >　 80  
      >　 100  
      >　 140  
      >　 180  
      >　 270  
      >　 360  
      >　 540  
      
      `~/XCASE/spec_ana$bin/SearchCentroid.py rate_image CentX CentY z pointsource_reg.file`      
      SearchCentroid.pyを実行し、半径500 kpc内でカウントレートで重み付けした銀河団中心を決める。  
      * rate_image : EPICの複合カウントレートイメージ。  
      * CentX : 中心探査の初期値。rate_imageのcoordinateでのX座標  
      * CentY : 中心探査の初期値。rate_imageのcoordinateでのY座標
      * z     : 銀河団の赤方偏移。
      * pointsource_reg.file : 点源のリージョンファイル。
	
	>例) ~/XCASE/spec_ana$bin/SearchCentroid.py rate-400-2300-all.fits 450 450 0.189 points.reg

   5. スペクトル抽出   
      `~/XCASE/spec_ana$bin/make_spec.csh Ngrp`		
      指定した同心円円環領域のスペクトルを抽出する。   
      spec_ana/SetUpFile/GrpBin.datに、Ngrpをアウトプット。	
      >例)~/XCASE/spec_ana$bin/make_spec.csh 50	
				      
      いくつかのスペクトルの抽出に失敗した場合は、該当の領域のスペクトルを再度抽出する。  		
      `~/XCASE/spec_ana$bin/remake_spec.csh det rinn rout`					
      `~/XCASE/spec_ana$bin/run_grppha.csh Ngrp`	    
      >例)~/XCASE/spec_ana$bin/remake_spec.csh pn 0 20  
        　 ~/XCASE/spec_ana$bin/run_grppha.csh 50	

      他円環からの漏れ込みをモデルとして組みたい場合、cross arfを作成する。  
      `~/XCASE/spec_ana$bin/makecross_arf_mos1_part1.csh`  
      `~/XCASE/spec_ana$bin/makecross_arf_mos1_part2.csh`  
      `~/XCASE/spec_ana$bin/makecross_arf_mos1_part3.csh`  
      `~/XCASE/spec_ana$bin/makecross_arf_mos2_part1.csh`  
      `~/XCASE/spec_ana$bin/makecross_arf_mos2_part2.csh`  
      `~/XCASE/spec_ana$bin/makecross_arf_mos2_part3.csh`  
      `~/XCASE/spec_ana$bin/makecross_arf_pn_part1.csh`  
      `~/XCASE/spec_ana$bin/makecross_arf_pn_part2.csh`  
      `~/XCASE/spec_ana$bin/makecross_arf_pn_part3.csh`      
      中心5領域で相互にそれより外側の円環で隣り合う円環どうしのみの漏れこみを考慮したarfを作成する。
      6.のフィットコード作成では、中心4領域の相互漏れこみそれより外側の円環で隣り合う円環どうしの漏れこみを考慮したフィットコードを作成する。

      `~/XCASE/spec_ana$bin/arf_plot.py`  
      cross arfが正しく作成されているか確認する。作成されていないものは"******** does not exits"と表示される。  
      `~/XCASE/spec_ana$bin/remake_crossarf.csh prefix imagebin rinn1 rout1 rinn2 rout2`  
      作成されなかったcross arfを再度、作成する。  
      * prefix      : 検出器の種類(mos1 or mos2 or pn)
      * imagebin    : arf作成に用いるカウントイメージのビンサイズ(50~500)
      * rinn1-rout1 : 漏れ出す側の円環
      * rinn2-rout2 : 漏れ込む側の円環        
     > 例）~/XCASE/spec_ana$bin/remake_crossarf.csh mos1 50 0-40-40-60

    グルーピングする。  
     `~/XCASE/spec_ana$bin/grp_corssarf.csh 50` 

     各円環領域のスペクトルのカウント数を確認しても良い。  
      `~/XCASE/spec_ana$bin/Calcount.py`					
      calc_******.datにカウント数を出力している。
      
   6. スペクトルフィットコード作成  
      `~/XCASE/spec_ana$bin/make_fitcode.py z nH FitcodeSampleFile/ample_num.xcm mode`	     
      XSPECのスペクトルフィッティングコードを作成する。  
      * z               : redshift
      * nH[10^22 cm^2]  : colum density　HEASRCのツールから21 cmの観測値を用いる。  
        URL : [HEASRC nH cluculation](https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl)   
      * sample_num.xcm  : 円環の数(num)に対応したサンプルフィットコード。
      * mode            : 0=usual mode, 1=cross arf mode  
     > 例）~/XCASE/spec_ana$bin/make_fitcode.py 0.182 0.0121 FitcodeSampleFile/sample_10_cross.xcm 0
     
      ROAST衛星から銀河団中心から1~2degの円環のスペクトルをバックグラウンドとして用いる。  
        URL : [HEASRC X-ray background Tool](http://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xraybg/xraybg.pl?Entry=13+49+00%2E50%2C+%2B26+35+07%2E0&CoordSys=J2000&NR=GRB%2FSIMBAD%2BSesame%2FNED&radius=2.00&region=annulus&inner_radius=1.0&spectrum=create&scaling=hist)   
      RASSスペクトルをgrpphaでグルーピングする。  
      `~/XCASE/spec_ana$grppha xrbg_orig.pi rass.pi`  
      `GRPPHA[]cheky RESPFILE pspcc.rsp`  
      
   7. スペクトルフィッティング  
      Xspecでフィッティングコードを読み込み、フィッティングする。  
      (フレアの影響が少ないデータはSoft proton backgroundのpowerlawのパラメータが定まらない場合があり、normを0にしてパラメータを固定してフィッティングを行う。）  
      `~/XCASE/spec_ana]$xspec`  
      `XSPEC12>@spectrum.xcm`  
      `XSPEC12>cpd /xs`  
      `XSPEC12>plot ldata ratio`  
      `XSPEC12>show free`   
      `XSPEC12>paral leven 8`    
      `XSPEC12>fit`  
      `XSPEC12>paral error 8`    
      `XSPEC12>tclout param 45`  
      `XSPEC12>echo "45 $xspec_tclout" >> observed_temp_value.dat`  
      `XSPEC12>error 1.0 45`  
      `XSPEC12>tclout error 45`  
      `XSPEC12>echo "45 $xspec_tclout" >> observed_temp_error.dat`  
      `XSPEC12>tclout param 192`  
      `XSPEC12>echo "192 $xspec_tclout" >> observed_temp_value.dat`  
      `XSPEC12>error 1.0 192`  
      `XSPEC12>tclout error 192`  
      `XSPEC12>echo "192 $xspec_tclout" >> observed_temp_error.dat`  
      `XSPEC12>save all save-all.xcm`  

   8. フレアの引き残り成分のイメージ作成に必要なスケールファクターの計算  
   Soft proton background image作成のため、 Caldb 内の Soft proton imageに対する Scaled Normalizationの計算結果を出力する。  
   引数は上で保存したsave file。  
   ※たとえSoft proton backgroundのpowerlawのnormを0にして固定し、モデルを無くした場合でも実行してください。  
   `~/XCASE/spec_ana$bin/sp_partial.csh save-all.xcm`

    output  
    log/sp_partial-mos1.log  
    log/sp_partial-mos2.log  
    log/sp_partial-pn.log  
    スクリプト内部では、自動でスペクトルフィットから得られるソフトプロトンバックグラウンドのpowerlawのnormalizationのベストフィットパラメータを読み込んである。 
    set mos1rnorm =  
    set mos2rnorm =   
    set pnrnorm =  


   9. フレアの引き残り成分のイメージ作成  
   `cd ../image_ana`  
   カレントディレクトリをimage_anaに変更する。

   proton.cshを実行し、ソフトプロトンの引きのこりのバックグラウンドイメージを作成する。    

   `~/XCASE/image_ana$bin/proton.csh`  

   スクリプト内部では以前のスクリプトのアウトプットを元にして以下の入力を自動してある。
    sp_partialコマンドの出力 log/sp_partial-mos1.log から、”Scaled Normalization” の値を読み込む。  
    set mos1pnorm =   
    set mos2pnorm =   
    set pnpnorm   =  
    スペクトルフィットから得られるソフトプロトンのバックグラウンドのpowerlawのphotoindexのベストフィットパラメータを読み込む。  
    set mos1pindex =  
    set mos2pindex =  
    set pnpindex  =  


   11 . イメージ作成  
   `~/XCASE/image_ana$bin/remake_image.csh`  
   remake_image.cshを実行し、ソフトプロトンの引きのこりのバックグラウンドを引いたイメージを作成する。  

   ※Al輝線1.5 keVを除きたい場合は、以下のタスクを実行する。
   `~/XCASE/image_ana$bin/proton_Alemit.csh`   
   `~/XCASE/image_ana$bin/remake_image_Alemit.csh` 
   

   12 . 銀河団中心のpsfを作成  
   `~/XCASE/image_ana$bin/psfgen.csh`  
   銀河団中心から半径 1.5 arcmin のPSF mapを作成する。  
   PSF mapは選択したエネルギーレンジの中心のエネルギー値で作成する。  
   > 例）400~2300 eVのエネルギーレンジの場合は1350 eV。
      
   13 . 正味の密度への単位変換に必要な係数の計算  
   xspecによるfakeitシミュレーションスクリプトを作成する  
   `cd ../cal_EI/`  
    カレントディレクトリをcal_EIに変更する。     
    `~/XCASE/cal_EI$bin/make_fake_xcm.py mos1preix z ../spec_ana/save-all.xcm`  
    `~/XCASE/cal_EI$bin/make_fake_xcm.py mos2preix z ../spec_ana/save-all.xcm`  
    `~/XCASE/cal_EI$bin/make_fake_xcm.py pnpreix z ../spec_ana/save-all.xcm`   
     make_fake_xcm.pyを実行し、xspecによるfakeitシミュレーションスクリプトを作成する。  
    * mos1prefix(mos2prefix,pnprefix) : type of EPIC  
    * z  : Redshift  
    * save-all.xcm : xspecのsave allコマンドで保存したスペクトルベストフィット       
    `~/XCASE/cal_EI$rm *fake*.pi`  
    `~/XCASE/cal_EI$csh fake_step1_mos1prefix.csh`  
    `~/XCASE/cal_EI$csh fake_step2_mos1prefix.csh`  
    `~/XCASE/cal_EI$csh fake_step1_mos2prefix.csh`  
    `~/XCASE/cal_EI$csh fake_step2_mos2prefix.csh`  
    `~/XCASE/cal_EI$csh fake_step1_pnprefix.csh`  
    `~/XCASE/cal_EI$csh fake_step2_pnprefix.csh`  
    fakeitシミュレーションスクリプトを実行し、シミュレーションスペクトルデータを出力する。      
    `~/XCASE/cal_EI$bin/lambda_cal.py mos1prefix z fake-spectrum ../spec_ana/save-all.xcm`  
    `~/XCASE/cal_EI$bin/lambda_cal.py mos2prefix z fake-spectrum ../spec_ana/save-all.xcm`  
    `~/XCASE/cal_EI$bin/lambda_cal.py pnprefix z fake-spectrum ../spec_ana/save-all.xcm`  
    lambda.pyを実行し、EIの変換係数を出力する。  
    * fake-spectrum : fake-data_mos1prefix.qdp or fake-data_mos2prefix.qdp or fake-data_pnprefix.qdp   
    `~/XCASE/cal_EI$bin/temp_cal.py z ../spec_dana/save-all.xcm ../spec_ana/observed_temp_error.dat`  
    temp_cal.pyを実行し、温度、アバンダンス、normalization、検出器間の系統差を示す定数を出力する。  
    * observed_temp_error.dat : xspecで計算した 1σ error を保存したテキストファイル  

   14 . 表面輝度分布およびプロジェクテッド温度のフィッティング  
   radial_fit/binディレクトリ内にあるfitting.cshに引数を使ってSx profileとT profileのfittingを行う。  
   Sx profileは銀河団中心から 2 acrmin の半径内でPSF convolutionしたモデルでfittingしている。  
   src 内のMakefileを編集し、適切なパスを通す。(デフォルトMakefile内の例では、HEASoftとGSLのライブラリを用いている。GSLがない場合はインストールが必要。）  
  
  カレントディレクトリをradial_fitに変更し、fitting.cshを実行する。  
  `~/XCASE$cd radial_fit`  
  `~/XCASE/radial_fit$bin/fittting.csh clustername z inner_radius outer_radius binning`  
  >例）~/XCASE/radial_fit$bin/fittting.csh ABELL1689 0.1832 0.1 10 22`  
    * name : 銀河団名（任意）  
    * rinn,rout,Nbin : フィット範囲の下限と上限とビン数  
    * z : Redshift  

    以下は自動で読み込むパラメーター
    * xc,yc : pixel単位の銀河団中心  
    * con_mos1,con_mos2,con_pn : result/constant.dat 内の各検出器の系統差を示す定数  
    
    outputはradial_fitのディレクトリに保存され、スクリプト内で../resultにコピーされる。


# アウトプットファイル



# トラブルシューティング
*filter.cshでEPIC/MOS1,MOS2,PNのprefix(e.g.mos1S001)は任意のprefixを取ってくるため、１つの観測データに複数のprefix(e.r.mos1U002)が含まれている場合、log/prefix.logの中身を解析で用いるprefixに書き換える必要がある。  
*cheese.cshで視野内に銀河団が複数存在する場合、解析に必要のないdifuseな銀河団放射を完全には除去できない。表面輝度フィッテイングコード内でdifuseな銀河団の領域のpixelのカウントレートを使わないようにするか、サブ銀河団としてモデルを加えるか。

# 参考文献
* [ESAS Cookbook](ftp://heasarc.gsfc.nasa.gov/xmm/software/xmm-esas/old/xmm-esas-v2/xmm-esas.pdf)
* 