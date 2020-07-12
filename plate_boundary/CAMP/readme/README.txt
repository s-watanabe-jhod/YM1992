########################################################################
############              README.txt  (2004/12/10)              ############
########################################################################

名称    ：CAMP Standard Model Version 1.0
参考文献：[1] Hashimoto, C., Fukui, K., and Matsu'ura, M., 3-D modelling of plate interfaces and numerical simulation of long-term crustal deformation in and around Japan, Pure and Applied Geophysics 161, 2053-2067, 2004. [2] Fukui, K., Hashimoto, C., Sato, T., Iwasaki, T., and Matsu'ura, M., 3-D standard model of plate interfaces in and around Japan, in preparation.
ファイル：camp_standard_model_1.0.tar.gz (README.txt, camp_pac.bdy, camp_phs.bdy, camp_pac.f, camp_phs.f)

本モデルの使用に際しては，以下の事項を確認してください．
１．本モデルは，Crustal Activity Modelling Program (CAMP)の成果として公開しているものです．
２．研究成果を発表する場合には，CAMP Standard Modelを使用したことを明記してください．
３．本モデルを営利目的に使用することを固く禁止します．
４．本モデルの使用により生じた問題に関しては，一切責任を負いません．
５．本モデルを第三者に再配布することを固く禁止します．

【モデルの仕様】
●モデル領域：東経125-155度，北緯20-50度，深さ0-100km
●座標系：直交直角座標系(x，y，z)
※球座標系(経度，緯度)から直交直角座標系(X，Y)への変換には，Lambert Conformal Conic projectionを用いています．その際のcentral pointは，東経140度・北緯35度とし，standard parallelsは，北緯30度・北緯40度としています．更に，鉛直下向きがz軸の正方向となるようにx=-Y，y=-Xとしています．単位はkmとしています．
●基準面(z=0)：固体地球の表面
※プレート境界面の深さzは，海に於いては海底面からの距離を，陸に於いては地表面からの距離を，それぞれ表しています．
●基底関数：bicubic B-spline関数
※プレート境界面の深さzを，基底関数の重ね合わせにより表現しています．
●解像度(bicubic B-spline関数の節点間隔)：x，y方向共に8km

【ファイルの内容】
●bicubic B-spline関数の重ね合わせ係数データ(camp_pac.bdy，camp_phs.bdy)
※以下の形式のデータとなっています．
1列目 k   ：x方向の基底関数(bicubic B-spline関数)の位置 [km]
2列目 l	  ：y方向の基底関数(bicubic B-spline関数)の位置 [km]
3列目 a   ：bicubic B-spline関数の重ね合わせ係数 [km]
4列目 flag：モデルの有効領域の指定
●プレート境界面の深さの計算コード(camp_pac.f，camp_phs.f)
※それぞれ，camp_pac.bdy，camp_phs.bdyのデータを読み込んで，任意の点(経度，緯度)に於けるプレート境界面の深さ[km]を計算します．例えば，camp_pac.fについて，
145 40
144 40
143 40
140.9 40
140.8 40
145 19
145 21
0 0 !END
のような形式の入力ファイル(lon_lat.inp)を用意し，
a.out<lon_lat.inp
とします．結果は，
 OUT OF THE MODEL REGION
 0.77397635289170652
 10.817832199967201
 99.328552590396129
 OUT OF THE MODEL REGION
 OUT OF THE MODEL REGION
 80.963499905803786
のように表示されます．
