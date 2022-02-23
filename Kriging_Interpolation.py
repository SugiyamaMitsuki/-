# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 17:44:06 2020

@author: mitsuki
"""
import pandas as pd
import os
import numpy as np
import sys
import pyproj
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from gstools import Gaussian
# from pykrige.ok import OrdinaryKriging
from matplotlib import pyplot as plt
import sympy
import math
import time
from estimate_Intensity_rogic import EST_Intensity



from config import Config


"""
距離減衰式に断層最短距離を使うか否かで結果が大きく変わる気がする。
断層最短距離の検討をするためには断層モデルが必要！！！

元データを変更する場合は
DF.itertuples()の箇所で列番号があっているか注意
"""

if __name__ == '__main__':
    ## 設定ファイル
    conf=Config()
    
    est_rogic = EST_Intensity()
    
    ## 設定
    ## クリギング法の参照範囲 [km]
    limit_distance = 15
    
    ## 気象庁は1kmメッシュだが250mメッシュを使用する
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
        
    print(JSHIS_ARV)
    
    
    
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
            organ=Target_line["機関"].values[0]
            
            code_list.append(code);lon_list.append(lon);lat_list.append(lat);I_list.append(I);amp_list.append(amp);organ_list.append(organ)
    code_df, lon_df, lat_df ,I_df, amp_df, organ_df = pd.Series(code_list),pd.Series(lon_list),pd.Series(lat_list),pd.Series(I_list),pd.Series(amp_list),pd.Series(organ_list)
     
    input_Observation = pd.concat([code_df, lon_df, lat_df, I_df, amp_df, organ_df], axis=1)
    input_Observation.columns = ['観測点番号','経度','緯度','計測震度',"地盤増幅度","機関"]
    
    
    
        

    ## 空の1kmメッシュを作成
    # lon_min,lon_max ,lat_min,lat_max = 136.25,137.25 ,34.75,35.75
    # lon_length = int((lon_max-lon_min)//(45/60/60))
    # lat_length = int((lat_max-lat_min)//(30/60/60))
    # delta_lon =45/60/60
    # delta_lat =30/60/60
    # cols = ['経度','緯度','地盤増幅度']
    # grid_1km = pd.DataFrame(index=[], columns=cols)
    # for x in range(lon_length):
    #     for y in range(lat_length):
    #         lon=lon_min+x*delta_lon
    #         lat=lat_min+y*delta_lat                
    #         distance = np.empty(len(ampdata))
    #         distance = pow(ampdata[:,0]-float(lon),2) + pow(ampdata[:,1]-float(lat),2)
    #         min_arg = np.argmin(distance)
    #         # min_distance = distance.min()
    #         amp = ampdata[min_arg][2]
    #         record = pd.Series([lon,lat,amp] ,  index=grid_1km.columns)
    #         grid_1km = grid_1km.append(record,ignore_index=True)
    #         print(lon,lat,amp)
    # print(grid_1km)
    # grid_1km.to_csv('../Estimation_list/grid_1km.csv',encoding="shift_jis",index=False)
    
    
    
    
    cols = ['経度','緯度','推計震度',"参照観測点数"]
    estimate_df = pd.DataFrame(index=[], columns=cols)  
    
    
    ## 計測震度から基盤面最大速度をもとめて観測記録のdfに追加する
    PV_list=[]
    fault_distance_list=[]
    hyp_distance_list=[]
    for i,A in enumerate(input_Observation.itertuples()):
        # print(A[4])
        I_input=A[4]
        amp = A[5]
        def f(PGV):
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
        
        #断層最短距離を計算
        X=0
        if conf.Distance_attenuation_type=="Shortest_fault_distance":
            X=est_rogic.calculate_Shortest_fault_distance(A[2], A[3])/1000
            fault_distance_list.append(X)
        else:fault_distance_list.append(X)
        grs80 = pyproj.Geod(ellps='GRS80')  # GRS80楕円体
        azimuth, bkw_azimuth, hyp_distance = grs80.inv(A[2], A[3], conf.epi_lon,conf.epi_lat)
        hyp_distance_list.append(hyp_distance/1000)
        
        
    input_Observation["PV"]=PV_list
    input_Observation["断層最短距離[km]"]=fault_distance_list
    input_Observation["震源距離[km]"]=hyp_distance_list
    
    input_Observation.to_csv("input_Observation.csv",encoding="shift_jis",index=False)
    
    
    PV_list_for_Ch=[]
    dx,distance_max=0.2,1000
    distance = np.arange(5, distance_max, dx)
    #等価震源距離
    for i in distance:
        Xeq=i
        PV=0
        # 等価震源距離
        if conf.Distance_attenuation_type=="Using_Equivalent_hypocentral_distance":
            PV = 10**(conf.a*conf.Mw + conf.h*conf.D + conf.d +conf.e -np.log10(Xeq) -0.002*Xeq)
        
        # 断層最短震源距離
        if conf.Distance_attenuation_type=="Shortest_fault_distance":
            PV = 10**(conf.a*conf.Mw + conf.h*conf.D + conf.d +conf.e -np.log10(Xeq+0.0028*10**(0.50*conf.Mw)) -0.002*Xeq)
        # print(PV)
        PV_list_for_Ch.append(PV)
    # plt.plot(distance,PV_list_for_Ch)
    
    #観測点間距離ごとに共分散を求める 配列番号が距離[km]になっている
    Ch=[]
    step_arr = np.arange(0,max(distance)-10, 1)            #numpy.arrange(a, b, d)メソッドでa以上b未満で間隔dの等差数列を作ることができます。   
    
    average=sum(PV_list_for_Ch)/len(PV_list_for_Ch)
    for step in step_arr:
        hh=int(step/dx)
        N=len(PV_list_for_Ch)
        data_list=[PV_list_for_Ch[i]+PV_list_for_Ch[i+hh] for i in range(N-hh)]
        average=sum(data_list)/(len(data_list)*2)
        # print("サンプル数",len(data_list),"平均",average)
        # print(sum([(PV_list_for_Ch[i]-average)*(PV_list_for_Ch[i+hh]-average) for i in range(N-hh)])/(N-hh))
        Ch.append(sum([(PV_list_for_Ch[i]-average)*(PV_list_for_Ch[i+hh]-average) for i in range(N-hh)])/(N-hh))
    
    
    ## ここまでは共通の操作
    
    for k,A in enumerate(JSHIS_ARV.itertuples()):
        # print(A)
        est_lon,est_lat,est_amp=A[2],A[3],A[5]
        
        ## 推定地点からの距離で推計に使用する観測点を絞る
        _lon_list = [est_lon]*len(input_Observation)
        _lat_list = [est_lat]*len(input_Observation)
        grs80 = pyproj.Geod(ellps='GRS80')  # GRS80楕円体
        azimuth, bkw_azimuth, distance_from_est = grs80.inv(input_Observation['経度'].values.tolist(), input_Observation['緯度'].values.tolist(), _lon_list, _lat_list) 
        input_Observation["distance_from_est"]=distance_from_est
        input_Observation_limited=input_Observation[input_Observation['distance_from_est'] < limit_distance*1000]
        
        ## i番目とすべての観測点間の距離
        grs80 = pyproj.Geod(ellps='GRS80')  # GRS80楕円体
        distance_list=[]
        for A in input_Observation_limited.itertuples():
            _lon_list = [est_lon]*len(input_Observation_limited)
            _lat_list = [est_lat]*len(input_Observation_limited)
            azimuth, bkw_azimuth, A_distance_list = grs80.inv(input_Observation_limited['経度'].values.tolist(), input_Observation_limited['緯度'].values.tolist(), _lon_list, _lat_list) 
            distance_list.append(A_distance_list)        
        distance_list = np.array(distance_list)/1000
        # print(distance_list)
        
        C= [[0 for j in range(len(input_Observation_limited)+2)] for i in range(len(input_Observation_limited)+2)]
        for i in range(len(input_Observation_limited)):
            for j in range(len(input_Observation_limited)):
                C[i][j]=Ch[int(distance_list[i][j])]
        
                
        for i,B in enumerate(input_Observation_limited.itertuples()): 
            
            ####
            ## 断層面を考える場合はここに工夫が必要
            ####
            # print(B)
            # sys.exit()
            azimuth, bkw_azimuth, hyp_distance = grs80.inv(B[2], B[3], conf.epi_lon,conf.epi_lat)
            epi_distance = pow(hyp_distance**2 + conf.depth**2,0.5)
            Xeq=epi_distance/1000
            
            PV=0
            # 等価震源距離
            if conf.Distance_attenuation_type=="Using_Equivalent_hypocentral_distance":
                PV = 10**(conf.a*conf.Mw + conf.h*conf.D + conf.d +conf.e -np.log10(Xeq) -0.002*Xeq)
            # 断層最短震源距離
            if conf.Distance_attenuation_type=="Shortest_fault_distance":
                X=B[8]
                PV = 10**(conf.a*conf.Mw + conf.h*conf.D + conf.d +conf.e -np.log10(X+0.0028*10**(0.50*conf.Mw)) -0.002*X)
            
            C[i][len(input_Observation_limited)]=PV
            C[len(input_Observation_limited)][i]=PV
            # print(A[1], A[2], epi_lon,epi_lat,Xeq,PV)
            C[i][len(input_Observation_limited)+1]=1
            C[len(input_Observation_limited)+1][i]=1
            
        # print(C)
    
        Cx= [0 for j in range(len(input_Observation_limited)+2)] 
    
        
        #推定地点とすべての観測点間の距離
        _lon_list = [est_lon]*len(input_Observation_limited)
        _lat_list = [est_lat]*len(input_Observation_limited)
        azimuth, bkw_azimuth, distance_list_Cx = grs80.inv(input_Observation_limited['経度'].values.tolist(), input_Observation_limited['緯度'].values.tolist(), _lon_list, _lat_list)     
        distance_list_Cx.append(0)
        distance_list_Cx.append(0)
        distance_list_Cx = np.array(distance_list_Cx)/1000
        
        for i in range(len(input_Observation_limited)):
            Cx[i]=Ch[int(distance_list_Cx[i])]
        
        ####
        ## 断層面を考える場合はここに工夫が必要
        ####
        
        azimuth, bkw_azimuth, hyp_distance = grs80.inv(est_lon,est_lat, conf.epi_lon,conf.epi_lat)
        epi_distance = pow(hyp_distance**2 + conf.depth**2,0.5)
        Xeq=epi_distance/1000
        
        PV=0
        # 等価震源距離
        if conf.Distance_attenuation_type=="Using_Equivalent_hypocentral_distance":
            PV = 10**(conf.a*conf.Mw + conf.h*conf.D + conf.d +conf.e -np.log10(Xeq) -0.002*Xeq)
        # 断層最短震源距離
        if conf.Distance_attenuation_type=="Shortest_fault_distance":
            X=est_rogic.calculate_Shortest_fault_distance(est_lon,est_lat)/1000
            PV = 10**(conf.a*conf.Mw + conf.h*conf.D + conf.d +conf.e -np.log10(X+0.0028*10**(0.50*conf.Mw)) -0.002*X)
        
        Cx[len(input_Observation_limited)]=PV
        Cx[len(input_Observation_limited)+1]=1
    
        # print(a,Mw,h,D,d,e)
        # print(0.58,5.8,0.0031,29,0.06,-1.25)
        # sys.exit()
        # print(Cx)
        # print(C)
        # C_inv = np.linalg.inv(C)
        C_inv = np.linalg.pinv(C)
        b=np.dot(Cx, C_inv)
        # print(b)
        # print("重みの和",sum(b[:-2]))
        
        # np.savetxt('C_savetxt_3d.txt', C)
        # np.savetxt('Cx_savetxt_3d.txt', Cx)
        # np.savetxt('b_savetxt_3d.txt',b)
        
        
        # print(input_Observation_limited)
        # print(b)
        Est_PV = np.dot(b[:-2], input_Observation_limited["PV"].values)
        # print(Est_PV)
        Est_PGV=Est_PV*est_amp
        Est_Intensity = 3.383-0.165*conf.Mw + 2.254* np.log10(Est_PGV) -0.082* np.log10(Est_PGV)**2
        print(k,"推計震度",round(Est_Intensity,2),round(est_lon,5),round(est_lat,5),len(input_Observation_limited),"重みの和",sum(b[:-2]))
        
    
        record = pd.Series([est_lon,est_lat,Est_Intensity,len(input_Observation_limited)],index=estimate_df.columns)        
        estimate_df = estimate_df.append(record,ignore_index=True)     
    
        # PV1 = 10**(conf.a*conf.Mw + conf.h*conf.D + conf.d +conf.e -np.log10(Xeq) -0.002*Xeq)
        # I1=3.383-0.165*conf.Mw + 2.254* np.log10(PV1) -0.082* np.log10(PV1)**2
        # X=est_rogic.calculate_Shortest_fault_distance(est_lon,est_lat)/1000
        # PV2 = 10**(conf.a*conf.Mw + conf.h*conf.D + conf.d +conf.e -np.log10(X+0.0028*10**(0.50*conf.Mw)) -0.002*X)
        # I2=3.383-0.165*conf.Mw + 2.254* np.log10(PV2) -0.082* np.log10(PV2)**2
        # print(PV1,PV2,round(PV2-PV1,5),I2-I1,"\n")
        
    print(estimate_df)
    # estimate_df.to_csv("kriging.csv",encoding="shift_jis",index=False)
    estimate_df.to_csv('{}/est_kriging_mesh_{}_{}_{}_{}.csv'.format(conf.output_path,conf.key,conf.mesh_size,conf.Distance_attenuation_type,limit_distance),encoding="shift_jis",index=False)
    estimate_df.to_csv('{}/est_kriging_mesh_{}_{}_{}_{}.txt'.format(conf.output_path,conf.key,conf.mesh_size,conf.Distance_attenuation_type,limit_distance),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])
    estimate_df_over4=estimate_df[estimate_df["推計震度"]>3.5]
    estimate_df_over4.to_csv('{}/est_kriging_mesh_{}_{}_{}_over4_{}.txt'.format(conf.output_path,conf.key,conf.mesh_size,conf.Distance_attenuation_type,limit_distance),encoding="shift_jis",index=False, sep='\t',header=None, columns=['経度','緯度','推計震度'])

    print(estimate_df["推計震度"].max)
    
    
    # fig = plt.figure()
    # ax = Axes3D(fig)
    # ax.set_xlabel("lon")
    # ax.set_ylabel("lat")
    # ax.set_zlabel("I")
    # ax.scatter(estimate_df['経度'],estimate_df['緯度'],estimate_df['推計震度'],marker="o",c="g",alpha=0.5)
    # # ax.scatter(np.array(lon_list),np.array(lat_list),np.array(depth_list)*(-1),marker="o",c="w")
    # ax.scatter(conf.epi_lon,conf.epi_lat,0,marker="*",c="r",s=8)
    # ax.text(conf.epi_lon,conf.epi_lat,0, "Epicenter")
    # plt.show()





    
    
    
    """
    
    # conditioning data
    data = np.array([[0.3, 1.2, 40],
                     [1.9, 0.6, 1],
                     [1.1, 3.2, 100]])
    # grid definition for output field
    gridx = np.arange(0.0, 5.5, 0.1)
    gridy = np.arange(0.0, 6.5, 0.1)
    # a GSTools based covariance model
    cov_model = Gaussian(dim=2, len_scale=1, anis=0.2, angles=-0.5, var=0.5, nugget=0.1)
    # ordinary kriging with pykrige
    OK1 = OrdinaryKriging(data[:, 0], data[:, 1], data[:, 2], cov_model)
    z1, ss1 = OK1.execute('grid', gridx, gridy)
    plt.imshow(z1, origin="lower")
    plt.show()
    
    
    print(cov_model)
    
    
    
    
    
    """
