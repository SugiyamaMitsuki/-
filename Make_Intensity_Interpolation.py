# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:14:57 2020

@author: sugiyama

クラスを使おうとしたけど逆にわかりにくくなったかも
"""

import matplotlib.pyplot as plt
import pandas as pd
import os
import folium
import glob
import numpy as np
import sys
import pyproj
from scipy.spatial import Delaunay, delaunay_plot_2d, Voronoi, voronoi_plot_2d
from sympy.geometry import Point, Polygon
import sympy as sp
import math
import time
import warnings
warnings.filterwarnings('ignore')


from config import Config
from estimate_Intensity_rogic import EST_Intensity
    

if __name__ == '__main__':
    
    ## 設定ファイル
    conf=Config()
    
    
    ## 震度データを読み込んでDFに変換
    code_list, lon_list, lat_list, I_list, amp_list, organ_list = [],[],[],[],[],[]
    with open(conf.input_path) as f:
        for i,line in enumerate(f):
            if i==0: 
                continue
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
    code_df, lon_df, lat_df ,I_df, amp_df, organ_df = pd.Series(code_list),pd.Series(lon_list),pd.Series(lat_list),pd.Series(I_list),pd.Series(amp_list),pd.Series(organ_list)
     
    input_Observation = pd.concat([code_df, lon_df, lat_df, I_df, amp_df, organ_df], axis=1)
    input_Observation.columns = ['観測点番号','経度','緯度','計測震度',"地盤増幅度","機関"]
    
    
    ## 計測震度から基盤面最大速度をもとめて観測記録のdfに追加する
    PV_list=[]
    for i,A in enumerate(input_Observation.itertuples()):
        # print(A[4])
        I_input=A[4]
        amp = A[5]
        def f(PGV):
            ## 計測震度と最大速度の関係式
            return 3.383-0.165*conf.Mw + 2.254* math.log10(PGV) -0.082* math.log10(PGV)**2
        up,low=200,0.001
        err=0.001
        count=0
        while True:
            count=count+1
            mid=(up+low)/2
            y=f(mid)-I_input
            if abs(y)<err: #解が発見された
                break
            elif (f(low)-I_input)*y<0: #解は下限と中間点の間にある
                up=mid
            else:                     #解は上限と中間点の間にある
                low=mid
        _PV = mid/amp if amp !=0 else 0
        PV_list.append(_PV)
        # print("数値解は",_PV)
    input_Observation["PGV"]=PV_list    
    
    ## 地表最大速度を地盤増幅度で除して基盤面最大速度を求める。
    input_Observation=input_Observation[input_Observation['地盤増幅度']>0]
    input_Observation["PGV_b600"]=input_Observation["PGV"]/input_Observation["地盤増幅度"]
    
    print("input_Observation")
    print(input_Observation,"\n")
    input_Observation.to_csv("input_Observation.csv",encoding="shift_jis",index=False)
    input_Observation.to_csv('{}/obs_plot.txt'.format(conf.output_path),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','計測震度'])
    
    ## 推計
    est_rogic = EST_Intensity()
    
    
    if conf.est:
        
        input_Observation_JMA_KNET=input_Observation[input_Observation['機関'] != "自治体"]
        input_Observation_JMA_KNET = input_Observation_JMA_KNET.reset_index()
        input_Observation_KNET=input_Observation[input_Observation['機関'] == "防災科研"]
        input_Observation_KNET = input_Observation_KNET.reset_index()
        
        #気象庁 防災科研 自治体

        ## 最短距離からの推計
        est_by_shortest_distance = est_rogic.est_Intensity_B_to_A(B_Obs=input_Observation,A_Obs=input_Observation)
        est_by_shortest_distance.to_csv('{}/est_by_shortest_distance_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        est_by_shortest_distance.to_csv('{}/est_by_shortest_distance_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])

        # est_by_shortest_distance_JMA_KNET = est_rogic.est_Intensity_B_to_A(B_Obs=input_Observation_JMA_KNET,A_Obs=input_Observation)
        # est_by_shortest_distance_JMA_KNET.to_csv('{}/est_by_shortest_distance_JMA_KNET_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        # est_by_shortest_distance_JMA_KNET.to_csv('{}/est_by_shortest_distance_JMA_KNET_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
        
        # est_by_shortest_distance_KNET = est_rogic.est_Intensity_B_to_A(B_Obs=input_Observation_KNET,A_Obs=input_Observation)
        # est_by_shortest_distance_KNET.to_csv('{}/est_by_shortest_distance_KNET_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        # est_by_shortest_distance_KNET.to_csv('{}/est_by_shortest_distance_KNET_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
        
        
        ## 観測点間距離が短い順に3点からの推計
        
        ## ドロネー図を使用した推計
        est_by_delaunay = est_rogic.est_Intensity_B_to_A_delaunay(B_Obs=input_Observation,A_Obs=input_Observation)
        est_by_delaunay.to_csv('{}/est_by_delaunay_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        est_by_delaunay.to_csv('{}/est_by_delaunay_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
    
        # est_by_delaunay_JMA_KNET = est_rogic.est_Intensity_B_to_A_delaunay(B_Obs=input_Observation_JMA_KNET,A_Obs=input_Observation)
        # est_by_delaunay_JMA_KNET.to_csv('{}/est_by_delaunay_JMA_KNET_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        # est_by_delaunay_JMA_KNET.to_csv('{}/est_by_delaunay_JMA_KNET_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
    
        # est_by_delaunay_KNET = est_rogic.est_Intensity_B_to_A_delaunay(B_Obs=input_Observation_KNET,A_Obs=input_Observation)
        # est_by_delaunay_KNET.to_csv('{}/est_by_delaunay_KNET_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        # est_by_delaunay_KNET.to_csv('{}/est_by_delaunay_KNET_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
    
    
        ##　ボロノイ図を使用した推計
        est_by_voronoi = est_rogic.est_Intensity_B_to_A_voronoi(B_Obs=input_Observation,A_Obs=input_Observation)
        est_by_voronoi.to_csv('{}/est_by_voronoi_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        est_by_voronoi.to_csv('{}/est_by_voronoi_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
    
        # est_by_voronoi_JMA_KNET = est_rogic.est_Intensity_B_to_A_voronoi(B_Obs=input_Observation_JMA_KNET,A_Obs=input_Observation)
        # est_by_voronoi_JMA_KNET.to_csv('{}/est_by_voronoi_JMA_KNET_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        # est_by_voronoi_JMA_KNET.to_csv('{}/est_by_voronoi_JMA_KNET_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
    
        # est_by_voronoi_KNET = est_rogic.est_Intensity_B_to_A_voronoi(B_Obs=input_Observation_KNET,A_Obs=input_Observation)
        # est_by_voronoi_KNET.to_csv('{}/est_by_voronoi_KNET_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        # est_by_voronoi_KNET.to_csv('{}/est_by_voronoi_KNET_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
    
    
    
    
    ## メッシュの推計震度を求める
    
    
    if conf.est_mesh:
        
        # conf.JSHIS_ARVをinput_Observationと同じDFにする
        # name,lon,lat,obs_I,AMP = A[1],A[2],A[3],A[4],A[5]
        
        ## 推計範囲のメッシュのみにする
        JSHIS_ARV = conf.JSHIS_ARV.copy()
        JSHIS_ARV.columns = ['経度','緯度','AVS30'] 
        JSHIS_ARV=JSHIS_ARV[JSHIS_ARV['経度']<conf.est_mesh_lon+conf.set_area]
        JSHIS_ARV=JSHIS_ARV[JSHIS_ARV['経度']>conf.est_mesh_lon-conf.set_area]
        JSHIS_ARV=JSHIS_ARV[JSHIS_ARV['緯度']<conf.est_mesh_lat+conf.set_area]
        JSHIS_ARV=JSHIS_ARV[JSHIS_ARV['緯度']>conf.est_mesh_lat-conf.set_area]
        JSHIS_ARV=JSHIS_ARV[JSHIS_ARV['AVS30']>10] #J-SHISの表層地盤情報が割り当てられていないものを排除
       
        ## AVS30をVs600m/sの地盤増幅度に変換する
        def AVS30_to_AMP(x):
            return 10 ** (2.367-0.852 * math.log10(x))
        JSHIS_ARV['地盤増幅度'] = JSHIS_ARV['AVS30'].apply(AVS30_to_AMP)
        JSHIS_ARV['観測点番号']=range(len(JSHIS_ARV))
        JSHIS_ARV['計測震度']=0
        JSHIS_ARV=JSHIS_ARV.reindex(columns=['観測点番号', '経度', '緯度', '計測震度', '地盤増幅度'])
        
    
        print("JSHIS_ARV","\n")
        print(JSHIS_ARV)
        
        ## 最短距離からの推計
        est_by_shortest_distance_mesh = est_rogic.est_Intensity_B_to_A(B_Obs=input_Observation,A_Obs=JSHIS_ARV)
        est_by_shortest_distance_mesh.to_csv('{}/est_by_shortest_distance_mesh_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        est_by_shortest_distance_mesh.to_csv('{}/est_by_shortest_distance_mesh_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
        
        ## 観測点間距離が短い順に3点からの推計
        
        ## ドロネー図を使用した推計
        est_by_delaunay_mesh = est_rogic.est_Intensity_B_to_A_delaunay(B_Obs=input_Observation,A_Obs=JSHIS_ARV)
        est_by_delaunay_mesh.to_csv('{}/est_by_delaunay_mesh_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        est_by_delaunay_mesh.to_csv('{}/est_by_delaunay_mesh_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
    
        ##　ボロノイ図を使用した推計
        est_by_voronoi_mesh = est_rogic.est_Intensity_B_to_A_voronoi(B_Obs=input_Observation,A_Obs=JSHIS_ARV)
        est_by_voronoi_mesh.to_csv('{}/est_by_voronoi_mesh_{}.csv'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False)
        est_by_voronoi_mesh.to_csv('{}/est_by_voronoi_mesh_{}.txt'.format(conf.output_path,conf.output_flag),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
    
    