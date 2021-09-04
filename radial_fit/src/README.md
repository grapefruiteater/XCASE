# Make

 `src]$make`
 
 radial_fit/binに実行ファイルができます。　
 コンパイルが上手くいかない場合はMakefileの中身を修正してください。
 
# Sx profile fittingに考慮しているもの。

 1. 温度の半径依存性による放射率の変化。温度測定がない半径は内側から外挿（２次の多項式）。
 2. mos1,mos2,pnのそれぞれのPSFで畳み込み
 3. mos1,mos2,pnの各結果と平均の値の差を系統誤差としてエラーに追加。
 4. Vikhilinin modelでパラメーターを固定しているものがある。スクリプト上部を参照。
 5. もし、パラメーターの変更などしたい場合はコード書き直す必要がある。
 6. beta modelとdobule beta modelのソースコードもあり。
 7. フィティングには長時間かかるので並列計算推奨。スクリプト内の
    `setenv OMP_NUM_THREADS 4`
    でスレッド数の変更が可能。

# T profile fittingに考慮しているもの。

 2次元温度プロファイルかspectroscopic-like ( or emission weight)で３次元温度プロファイルを復元
 
# 他スクリプト内のpythonスクリプトで静水圧平衡質量やガス質量を計算

  cosmology : H0=70km/s/Mpc Om0=0.3 OL=0.7
  
# 温度、密度の動径方向プロファイルのフィティングコード

 1.Sx profile fitting : 

 fit_sx_vikhlinin_mos1+mos2+pn.cxx : Vikhlinin model

 fit_sx_beta_mos1+mos2+pn.cxx : beta model

 fit_sx_2beta_mos1+mos2+pn.cxx : double beta model

 fit_sx_2clusters.cxx : simultaneous-fit for Sx profiles of 2 clusters with beta model
 
 fit_sx_2clusters_bkgfixed.cxx : simultaneous-fit for Sx profiles of 2 clusters with beta model (bkg fixed model)

 2.output :
 
 Sxresult/Lambda_profile.dat     - radial depednece of emissivity  
 Sxresult/ne_bestfit.dat         - best-fit n profile in arbital unit  
 Sxresult/sx_bestfit.dat         - best-fit sx profile (normalized)  
 Sxresult/sx_bestfit.log         - best-fit log  
 Sxresult/sx_bestfit_noncor.dat  - best-fit raw sx profile   
 Sxresult/sx_profile.dat         - obs data (normalized)  
 Sxresult/sx_profile_noncor.dat  - obs data  

 3.Tx profile fitting :

 fit_Tx_ew.cxx

fit_Tx_sl_1styear.cxx 
 
 output :
 Tresult/T3d_bestfit.dat    - 3D best-fit T profile  
 Tresult/Tx_bestfit.dat     - 2D best-fit T profile  
 Tresult/Tx_bestfit.log     - best-fit log file  
 Tresult/Tx_profile.dat     - Tobs profile  
 
 4.Calculation by cal_ne0.py
 
 nresult/Mgas.dat  　　　　 - gas mass profile  
 nresult/ne.dat             - electron density profile  
 nresult/rho.dat            - mass density profile  
 
 5.Calculation of M_HE : 
  
  cal_MHE.py 

   output :
   MHEresult/MHE_prof.dat    - M_HE  profile   
   MHEresult/Mgas_prof.dat   - M_gas profile   
   HEresult/fgas_prof.dat    - f_gas profile 
   MHEresult/MHE_Delta.dat   - M_HE  at r_Delta  
   MHEresult/Mgas_Delta.dat  - M_gas at r_Delta  
   MHEresult/fgas_Delta.dat  - f_gas at r_Delta  
  
  
  