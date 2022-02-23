# -*- coding: utf-8 -*-

import pandas as pd
import os
import pyproj
import math

import warnings
warnings.filterwarnings('ignore')

from config import Config
        

if __name__ == '__main__':
    
    ## 設定ファイル
    conf=Config()
    grs80 = pyproj.Geod(ellps='GRS80') 
    try:
        os.mkdir("../Intensity_data/")
    except:
        None
    
    ## 出力ファイルのデータフレームの作成
    cols = ['発生時刻','震央経度','震央緯度','震源深さ','気象庁マグニチュード','最大震度(階級)','最大震度(計測震度)','観測記録数','最大震度の震央距離']
    output = pd.DataFrame(index=[], columns=cols)    
    
    ## 年ごとに処理をする
    for target_year in conf.target_year:
        
        EQdata_path = '../Catalog/i{}.dat'.format(target_year)
        eq_key=""
        eq_key_list=[]
        
        ## 地震ごとに別々のtxtファイルに書き出す
        with open(EQdata_path) as f:
            for i,line in enumerate(f):
                if line[0]=="A" or line[0]=="B":
                    eq_key=line[1:16]
                    eq_key_list.append(eq_key)
                    try:
                        os.remove('../Intensity_data/{}.txt'.format(eq_key))
                    except:
                        None
                path_w = '../Intensity_data/{}.txt'.format(eq_key)
                with open(path_w, mode='a') as f:
                    f.write(line)
        
        
        ## 地震情報の読み込み    
        for eq_key in eq_key_list:
            try:
                
                input_path='../Intensity_data/{}.txt'.format(eq_key)
                
                ## 変数の定義
                epi_lon,epi_lat,depth,Mj,I_class_max,records_num=0,0,0,0,0,0
                code_list, lon_list, lat_list, I_list, amp_list, organ_list = [],[],[],[],[],[]
                input_Observation=None
                i=0
                epi_lat,epi_lon,depth,Mw = 0,0,0,0
                output_path ,output_flag = "",""
                
                Make_DF = True
                with open(input_path) as f:
                    for i,line in enumerate(f):
                        if i==0:
                            #震源情報をヘッダーから読み取る
                            
                            print(f"eq_key {eq_key}")
                            # print("line.split()",line.split())
                            #地震カタログ震度データのフォーマット
                            #https://www.data.jma.go.jp/svd/eqev/data/bulletin/data/shindo/format_j.pdf
                            
                            eqi_lat_int=line[22-1:24].replace(" ", "")
                            eqi_lat_dec1=line[25-1:26] if line[25-1:26] != "  " else 0
                            eqi_lat_dec2=line[27-1:28] if line[27-1:28] != "  " else 0
                            epi_lat=float(eqi_lat_int) + (float(eqi_lat_dec1)+float(eqi_lat_dec2)/100)/60
                            
                            eqi_lon_int=line[33-1:36]
                            eqi_lon_dec1=line[37-1:38] if line[37-1:38] != "  " else 0
                            eqi_lon_dec2=line[39-1:40] if line[39-1:40] != "  " else 0
                            epi_lon=float(eqi_lon_int) + (float(eqi_lon_dec1)+float(eqi_lon_dec2)/100)/60
                            
                            depth_int=line[45-1:47].replace(" ", "",3) if line[45-1:47] != "   " else 9999
                            depth_dic=line[48-1:49].replace(" ", "",3) if line[48-1:49] != "  " else 9999
                            if depth_int == 9999 or depth_dic == 9999:
                                depth = 9999
                            else:
                                depth=float(depth_int) + float(depth_dic) /100 
                            
                            Mj = float(line[53-1:54])/10 if line[53-1:54] != "  " else 9999
                            
                            I_class_max = line[62-1:62]
                            # 最大震度(階級)
                            if I_class_max == "1":
                                I_class_max = "震度1"
                            elif I_class_max == "2":
                                I_class_max = "震度2"
                            elif I_class_max == "3":
                                I_class_max = "震度3"
                            elif I_class_max == "4":
                                I_class_max = "震度4"
                            elif I_class_max == "5":
                                I_class_max = "震度5"
                            elif I_class_max == "6":
                                I_class_max = "震度6"
                            elif I_class_max == "7":
                                I_class_max = "震度7"
                            elif I_class_max == "A":
                                I_class_max = "震度5弱"
                            elif I_class_max == "B":
                                I_class_max = "震度5強"
                            elif I_class_max == "C":
                                I_class_max = "震度6弱"
                            elif I_class_max == "D":
                                I_class_max = "震度6強"
                                
                            records_num = line.split()[-1].replace("K","") #line[91-1:95]
                            
                            print("震源情報 ")
                            print(f"経度{round(epi_lon,3)} 緯度{round(epi_lat,3)} 深さ{depth}km 気象庁マグニチュード{Mj} 震度階級{I_class_max} 観測記録数{records_num}")
                            
                            continue
                        
                        ## ヘッダー以外の処理
                        ##　各観測点の記録の読み込み
                        code=line.split()[0]#.replace(" ","")
                        I=float(line.split()[3])/10
                        Target_index=conf.matchingtable[conf.matchingtable["観測点番号"]==int(code)].index
                        Target_line = conf.matchingtable.loc[Target_index]
                        lon=float(Target_line["経度"])
                        lat=float(Target_line["緯度"])
                        AVS30=float(Target_line["AVS"])
                        # amp=float(Target_line["ARV"])
                        if AVS30 == 0:
                            continue
                        amp = 10 ** (2.367-0.852 * math.log10(AVS30)) if AVS30 != 0 else 0
                        # print(Target_line["機関"].values)
                        # print(str(Target_line["機関"]))
                        # print(str(Target_line["機関"]).split(" "))
                        organ=Target_line["機関"].values[0]
                        
                        code_list.append(code);lon_list.append(lon);lat_list.append(lat);I_list.append(I);amp_list.append(amp);organ_list.append(organ)
                    
                    
                    ## 最大震度の観測点を選択
                    max_I_index = I_list.index(max(I_list))
                    
                    lon_target = lon_list[max_I_index]
                    lat_target = lat_list[max_I_index]
                    
                    azimuth, bkw_azimuth, hyp_distance = grs80.inv(epi_lon,epi_lat,lon_target,lat_target)
                    
                    print(f"計測震度の最大値 {I_list[max_I_index]}")
                    print(f"最大震度の震央距離 {round(hyp_distance/1000)}km")
                    print()
                    
                    
                    record = pd.Series([eq_key,
                                        epi_lon,epi_lat,depth,Mj,
                                        I_class_max,I_list[max_I_index],
                                        records_num,
                                        hyp_distance/1000
                                        ] ,  index=output.columns)
                    output = output.append(record,ignore_index=True)
                
                    
            except Exception as e:
                print(str(e))
    
    output.to_csv(f'../震度情報_{min(conf.target_year)}~{max(conf.target_year)}年.csv',encoding="shift_jis",index=False) 
    
    
                
            #     None    
            #     if Make_DF:
            #         code_df, lon_df, lat_df ,I_df, amp_df, organ_df = pd.Series(code_list),pd.Series(lon_list),pd.Series(lat_list),pd.Series(I_list),pd.Series(amp_list),pd.Series(organ_list)
            #         input_Observation = pd.concat([code_df, lon_df, lat_df, I_df, amp_df, organ_df], axis=1)
            #         input_Observation.columns = ['観測点番号','経度','緯度','計測震度',"地盤増幅度","機関"]
                    
                    
            #         if input_Observation["計測震度"].max() >= conf.max_I_lowerlimit:
            #             try:
            #                 os.mkdir(output_path)
            #             except:
            #                 None
                            
            #             print(eq_key)
                            
            #             ## 計測震度から基盤面最大速度をもとめて観測記録のdfに追加する
            #             PV_list=[]
            #             for i,A in enumerate(input_Observation.itertuples()):
            #                 # print(A[4])
            #                 I_input=A[4]
            #                 amp = A[5]
            #                 def f(PGV):
            #                     ## 計測震度と最大速度の関係式
            #                     return 3.383-0.165*Mw + 2.254* math.log10(PGV) -0.082* math.log10(PGV)**2
            #                 up,low=200,0.001
            #                 err=0.001
            #                 count=0
            #                 while True:
            #                     count=count+1
            #                     mid=(up+low)/2
            #                     y=f(mid)-I_input
            #                     if abs(y)<err: #解が発見された
            #                         break
            #                     elif (f(low)-I_input)*y<0: #解は下限と中間点の間にある
            #                         up=mid
            #                     else:                     #解は上限と中間点の間にある
            #                         low=mid
            #                 _PV = mid/amp if amp !=0 else 0
            #                 PV_list.append(_PV)
            #                 # print("数値解は",_PV)
            #             input_Observation["PGV"]=PV_list    
                        
            #             ## 地表最大速度を地盤増幅度で除して基盤面最大速度を求める。
            #             input_Observation=input_Observation[input_Observation['地盤増幅度']>0]
            #             input_Observation["PGV_b600"]=input_Observation["PGV"]/input_Observation["地盤増幅度"]
                        
            #             print("input_Observation")
            #             print(input_Observation,"\n")
            #             # input_Observation.to_csv("input_Observation.csv",encoding="shift_jis",index=False)
            #             input_Observation.to_csv('{}/obs_plot.txt'.format(output_path),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','計測震度'])
                        
            #             ## 推計
            #             est_rogic = EST_Intensity()
                        
            #             flag = 1    
            #             if flag == 1:
                            
            #                 # input_Observation_JMA_KNET=input_Observation[input_Observation['機関'] != "自治体"]
            #                 # input_Observation_JMA_KNET = input_Observation_JMA_KNET.reset_index()
            #                 # input_Observation_KNET=input_Observation[input_Observation['機関'] == "防災科研"]
            #                 # input_Observation_KNET = input_Observation_KNET.reset_index()
                            
            #                 ## 推計される観測点を絞る
            #                 input_Observation_A=input_Observation.copy()
            #                 _lon_list = [epi_lon]*len(input_Observation_A)
            #                 _lat_list = [epi_lat]*len(input_Observation_A)
            #                 azimuth, bkw_azimuth, distance_list = grs80.inv(input_Observation_A['経度'].values.tolist(), input_Observation_A['緯度'].values.tolist(), _lon_list, _lat_list) 
            #                 input_Observation_A["震央距離"]=distance_list
            #                 input_Observation_A=input_Observation_A[input_Observation_A["震央距離"]<conf.distance_range*1000 * 3 ]#少し広めに設定する
            #                 # input_Observation_A=input_Observation_A.reset_index()
            #                 print("input_Observation_A")
            #                 print(input_Observation_A,"\n")
                            
            #                 if len(input_Observation_A)>2:
            #                     #気象庁 防災科研 自治体
                        
            #                     ## 最短距離からの推計
            #                     est_by_shortest_distance = est_rogic.est_Intensity_B_to_A(B_Obs=input_Observation,A_Obs=input_Observation_A)
            #                     est_by_shortest_distance.to_csv('{}/est_by_shortest_distance_{}.csv'.format(output_path,output_flag),encoding="shift_jis",index=False)
            #                     est_by_shortest_distance.to_csv('{}/est_by_shortest_distance_{}.txt'.format(output_path,output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
                        
            #                     # est_by_shortest_distance_JMA_KNET = est_rogic.est_Intensity_B_to_A(B_Obs=input_Observation_JMA_KNET,A_Obs=input_Observation)
            #                     # est_by_shortest_distance_JMA_KNET.to_csv('{}/est_by_shortest_distance_JMA_KNET_{}.csv'.format(output_path,output_flag),encoding="shift_jis",index=False)
            #                     # est_by_shortest_distance_JMA_KNET.to_csv('{}/est_by_shortest_distance_JMA_KNET_{}.txt'.format(output_path,output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
                                
            #                     # est_by_shortest_distance_KNET = est_rogic.est_Intensity_B_to_A(B_Obs=input_Observation_KNET,A_Obs=input_Observation)
            #                     # est_by_shortest_distance_KNET.to_csv('{}/est_by_shortest_distance_KNET_{}.csv'.format(output_path,output_flag),encoding="shift_jis",index=False)
            #                     # est_by_shortest_distance_KNET.to_csv('{}/est_by_shortest_distance_KNET_{}.txt'.format(output_path,output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
                                
                                
            #                     ## 観測点間距離が短い順に3点からの推計
                                
            #                     ## ドロネー図を使用した推計
            #                     est_by_delaunay = est_rogic.est_Intensity_B_to_A_delaunay(B_Obs=input_Observation,A_Obs=input_Observation_A)
            #                     est_by_delaunay.to_csv('{}/est_by_delaunay_{}.csv'.format(output_path,output_flag),encoding="shift_jis",index=False)
            #                     est_by_delaunay.to_csv('{}/est_by_delaunay_{}.txt'.format(output_path,output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
                            
            #                     # est_by_delaunay_JMA_KNET = est_rogic.est_Intensity_B_to_A_delaunay(B_Obs=input_Observation_JMA_KNET,A_Obs=input_Observation)
            #                     # est_by_delaunay_JMA_KNET.to_csv('{}/est_by_delaunay_JMA_KNET_{}.csv'.format(output_path,output_flag),encoding="shift_jis",index=False)
            #                     # est_by_delaunay_JMA_KNET.to_csv('{}/est_by_delaunay_JMA_KNET_{}.txt'.format(output_path,output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
                            
            #                     # est_by_delaunay_KNET = est_rogic.est_Intensity_B_to_A_delaunay(B_Obs=input_Observation_KNET,A_Obs=input_Observation)
            #                     # est_by_delaunay_KNET.to_csv('{}/est_by_delaunay_KNET_{}.csv'.format(output_path,output_flag),encoding="shift_jis",index=False)
            #                     # est_by_delaunay_KNET.to_csv('{}/est_by_delaunay_KNET_{}.txt'.format(output_path,output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
                            
                            
            #                     ##　ボロノイ図を使用した推計
            #                     est_by_voronoi = est_rogic.est_Intensity_B_to_A_voronoi(B_Obs=input_Observation,A_Obs=input_Observation_A)
            #                     est_by_voronoi.to_csv('{}/est_by_voronoi_{}.csv'.format(output_path,output_flag),encoding="shift_jis",index=False)
            #                     est_by_voronoi.to_csv('{}/est_by_voronoi_{}.txt'.format(output_path,output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
                            
            #                     # est_by_voronoi_JMA_KNET = est_rogic.est_Intensity_B_to_A_voronoi(B_Obs=input_Observation_JMA_KNET,A_Obs=input_Observation)
            #                     # est_by_voronoi_JMA_KNET.to_csv('{}/est_by_voronoi_JMA_KNET_{}.csv'.format(output_path,output_flag),encoding="shift_jis",index=False)
            #                     # est_by_voronoi_JMA_KNET.to_csv('{}/est_by_voronoi_JMA_KNET_{}.txt'.format(output_path,output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
                            
            #                     # est_by_voronoi_KNET = est_rogic.est_Intensity_B_to_A_voronoi(B_Obs=input_Observation_KNET,A_Obs=input_Observation)
            #                     # est_by_voronoi_KNET.to_csv('{}/est_by_voronoi_KNET_{}.csv'.format(output_path,output_flag),encoding="shift_jis",index=False)
            #                     # est_by_voronoi_KNET.to_csv('{}/est_by_voronoi_KNET_{}.txt'.format(output_path,output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
            # except:
            #     None                
                            
                            
                            