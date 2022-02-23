# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:14:57 2020

@author: sugiy
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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


class EST_Intensity(Config):
    def __init__(self):
        super(EST_Intensity, self).__init__()
        
        self.grs80 = pyproj.Geod(ellps='GRS80')  # GRS80楕円体
    
        
    def inside_or_outside_judgment(self,Polygon_vertex_coordinates,Judgmentpoint):
        x_list ,y_list= [],[]
        for xy in Polygon_vertex_coordinates:
            x_list.append(xy[0]);y_list.append(xy[1])
        A,B=[Judgmentpoint[0],Judgmentpoint[1]], [min(x_list),Judgmentpoint[1]]
        
        crossing_count=int(0)
        # #AB×ACとAB×ADを計算する
        def f(p1, p2, p3):
            return (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])
        for i ,coordinates in  enumerate(Polygon_vertex_coordinates):
            C = Polygon_vertex_coordinates[i]
            D = Polygon_vertex_coordinates[0] if i+1 == len(Polygon_vertex_coordinates) else Polygon_vertex_coordinates[i+1]
            t1 = f(A, B, C)
            t2 = f(A, B, D)
            t3 = f(C, D, A)
            t4 = f(C, D, B)
            if t1 * t2 < 0.0 and t3 * t4 < 0.0:
                crossing_count+=1
        Judgment = False if crossing_count != 1 else True
        
        return Judgment
    
    def calculate_Shortest_fault_distance(self,lon,lat):
        fault_Panels_df=[]
        for m, fault_param in enumerate(self.Fault_model_list):
            #self.epi_lon,self.epi_lat,self.depth 震源 破壊開始点
            fault_lon,fault_lat,fault_depth,strike,dip,fault_length,fault_width,length_slide,width_slide=fault_param[0],fault_param[1],fault_param[2],fault_param[3],fault_param[4],fault_param[5],fault_param[6],fault_param[7],fault_param[8]
            # print(fault_lon,fault_lat,strike,dip,fault_length,fault_width,length_slide,width_slide)
            
            fault_depth=fault_depth/1000
            # 極半径Polar radius (km)
            Polar_radius = 6356.752314
            # 赤道半径Equatorial radius (km)
            Equatorial_radius = 6378.137
            lat_per_1km =360/( 2 * Polar_radius * np.pi )
            lon_per_1km =360/( 2 * Equatorial_radius * np.pi * np.cos(np.radians(fault_lat)) )
    
            
            ## 断層面の角の平面座標（緯度経度）Basis_pointを調べる                   
            w = (fault_width - width_slide)
            l = (fault_length - length_slide)
            w2 = w * np.cos(np.radians(dip)) #真上から見た幅
            depth_basis = fault_depth + w * np.sin(np.radians(dip))
            
            # diagonal = np.sqrt( (w*np.cos(np.radians(dip)))**2 + l**2 )
            lat_basis = fault_lat + \
                        (- w2 * np.sin(np.radians(strike)) + l * np.cos(np.radians(strike))) *lat_per_1km
            lon_basis = fault_lon +\
                        (w2 * np.cos(np.radians(strike)) + l * np.sin(np.radians(strike))) * lon_per_1km
            
            ## 断層面をself.fault_mesh_splitmに分割
            mesh_length = self.fault_mesh_split
            length_count = int(fault_length * 1000 / mesh_length)
            width_count = int(fault_width * 1000 / mesh_length)
            
            fault_Panels = np.zeros((length_count,width_count,3))
            lon_list,lat_list,depth_list = [],[],[]
            
            #i==0
            fault_Panels[0, 0, 0] = lat_basis;
            fault_Panels[0, 0, 1] = lon_basis;
            fault_Panels[0, 0, 2] = depth_basis;
            
            # print(fault_Panels)
            #jが幅方向,iが長さ方向
            for i in range(length_count):
                if i>0:
                    #列が入れ替わったときj=0について
                    fault_Panels[i, 0, 0] = fault_Panels[i - 1, 0, 0] - mesh_length * 0.001* np.cos(np.radians(strike)) * lat_per_1km;
                    fault_Panels[i, 0, 1] = fault_Panels[i - 1, 0, 1] - mesh_length * 0.001* np.sin(np.radians(strike)) * lon_per_1km;
                    fault_Panels[i, 0, 2] = fault_Panels[i - 1, 0, 2];
                    lat_temp,lon_temp,depth_temp=fault_Panels[i, 0, 0],fault_Panels[i, 0, 1],fault_Panels[i, 0, 2]
                    lon_list.append(lon_temp)
                    lat_list.append(lat_temp)
                    depth_list.append(depth_temp)
                
                for j in range(1,width_count):
                    fault_Panels[i, j, 0] = fault_Panels[i , j - 1, 0] + mesh_length * 0.001* np.cos(np.radians(dip)) * np.sin(np.radians(strike)) * lat_per_1km
                    fault_Panels[i, j, 1] = fault_Panels[i , j - 1, 1] - mesh_length * 0.001* np.cos(np.radians(dip)) * np.cos(np.radians(strike)) * lon_per_1km
                    fault_Panels[i, j, 2] = fault_Panels[i , j - 1, 2] - mesh_length * 0.001* np.sin(np.radians(dip))
    
                    lat_temp,lon_temp,depth_temp=fault_Panels[i, j, 0],fault_Panels[i, j, 1],fault_Panels[i, j, 2]
                    lon_list.append(lon_temp)
                    lat_list.append(lat_temp)
                    depth_list.append(depth_temp)
            
            df = pd.DataFrame({'緯度':lat_list,
                               '経度':lon_list,
                               '深さ':depth_list#km
                               })        
            fault_Panels_df.append(df)
        
        ## 断層最短距離を計算する
        X_candidate=[]
        for df in fault_Panels_df:
            _lon_list = [lon]*len(df)
            _lat_list = [lat]*len(df)
            azimuth, bkw_azimuth, hyp_distance_list = self.grs80.inv(df['経度'].values.tolist(), df['緯度'].values.tolist(), _lon_list, _lat_list)
            epi_distance_list=[]
            for dis,dep in zip(hyp_distance_list,df['深さ'].values.tolist()):
                epi_distance_list.append(pow(dis**2 + (dep*1000)**2,0.5))
            X_candidate.append(min(epi_distance_list))
        X=min(X_candidate)
        # print(epi_distance_list)
        
        # plt.figure(figsize=(5,6),dpi=300,facecolor='w',edgecolor='k')
        
        # fig = plt.figure(dpi=300)
        # ax = Axes3D(fig)
        # ax.set_xlabel("Lon.")
        # ax.set_ylabel("Lat.")
        # ax.set_zlabel("Depth [km]")
        # ax.scatter(fault_Panels_df[0]['経度'],fault_Panels_df[0]['緯度'],fault_Panels_df[0]['深さ']*(-1),marker="o",c="g",alpha=0.5)
        # ax.scatter(fault_Panels_df[1]['経度'],fault_Panels_df[1]['緯度'],fault_Panels_df[1]['深さ']*(-1),marker="o",c="b",alpha=0.5)
        # # ax.scatter(np.array(lon_list),np.array(lat_list),np.array(depth_list)*(-1),marker="o",c="w")
        # ax.scatter(self.epi_lon,self.epi_lat,self.depth*(-1)*0.001,marker="*",c="r",s=800)
        # ax.text(self.epi_lon,self.epi_lat,self.depth*(-1)*0.001, "Epicenter")
        
        # plt.show()
        # sys.exit()
                
                
        return X
        



        
        
    def get_ref_data(self,ref_B,epi_distance,AMP, lon, lat):
        """
        参照観測点のデータを用いて推計震度を計算する関数

        """
        grs80 = pyproj.Geod(ellps='GRS80')  # GRS80楕円体
        
        #reference
        ref_name,ref_lon,ref_lat,ref_amp=ref_B['観測点番号'],ref_B['経度'],ref_B['緯度'],ref_B['地盤増幅度']
        ref_sindo = ref_B["計測震度"]
        ref_PGVb600 = ref_B["PGV_b600"]
        
        azimuth, bkw_azimuth, ref_hyp_distance = grs80.inv(self.epi_lon, self.epi_lat, ref_lon, ref_lat)
        ref_epi_distance = pow(ref_hyp_distance**2 + self.depth**2,0.5)
        
        Est_Intensity_R_S, Est_Intensity_R, Est_Intensity_S = 0, 0, 0
        estPGVb600_R=0
        upR=0
        downR=0

        ## 距離減衰の補正
        ##　震源距離の逆数
        if self.Distance_attenuation_type ==  "Reciprocal_of_epicenter_distance":                
            estPGVb600_R = ref_PGVb600 * (ref_epi_distance/epi_distance)
            estPGV_R = estPGVb600_R * 1
            if self.Relationship_I_and_vel == "conversion":
                    Est_Intensity_R = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV_R) -0.082*math.log10(estPGV_R)**2
            if self.Relationship_I_and_vel == "simple":
                upR=ref_epi_distance
                downR=epi_distance
                Est_Intensity_R = ref_sindo + 2*math.log10(upR/downR) #mitsuki
        
            
        ##　震源距離の逆数 断層最短距離ver
        elif self.Distance_attenuation_type == "Reciprocal_of_Shortest_fault_distance":
            X_est=self.calculate_Shortest_fault_distance(lon,lat)/1000
            X_ref=self.calculate_Shortest_fault_distance(ref_lon, ref_lat)/1000
            estPGVb600_R = ref_PGVb600 * (X_ref/X_est)
            estPGV_R = estPGVb600_R * 1
            if self.Relationship_I_and_vel == "conversion":
                Est_Intensity_R = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV_R) -0.082*math.log10(estPGV_R)**2
            if self.Relationship_I_and_vel == "simple":
                upR=X_ref
                downR=X_est
                Est_Intensity_R = ref_sindo + 2*math.log10(upR/downR) #mitsuki
            
            
        ## 等価震源距離
        elif self.Distance_attenuation_type ==  "Using_Equivalent_hypocentral_distance":
            #推計地点
            Xeq = epi_distance/1000
            PV_est = 10**(self.a*self.Mw + self.h*self.D + self.d +self.e -np.log10(Xeq) -0.002*Xeq)
            
            #参照地点
            Xeq = ref_epi_distance/1000
            PV_ref = 10**(self.a*self.Mw + self.h*self.D + self.d +self.e -np.log10(Xeq) -0.002*Xeq)
            
            estPGVb600_R = ref_PGVb600 * (PV_est/PV_ref)
            estPGV_R = estPGVb600_R * 1
            
            if self.Relationship_I_and_vel == "conversion":
                Est_Intensity_R = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV_R) -0.082*math.log10(estPGV_R)**2
            if self.Relationship_I_and_vel == "simple":
                upR=PV_est
                downR=PV_ref
                Est_Intensity_R = ref_sindo + 2*math.log10(upR/downR) #mitsuki
            
            
        ## 断層最短距離(距離減衰式)
        elif self.Distance_attenuation_type ==  "Shortest_fault_distance":
            #推計地点
            X_est=self.calculate_Shortest_fault_distance(lon,lat)/1000
            PV_est = 10**(self.a*self.Mw + self.h*self.D + self.d +self.e -np.log10(X_est+0.0028*10**(0.50*self.Mw)) -0.002*X_est)
            
            #参照地点
            X_ref=self.calculate_Shortest_fault_distance(ref_lon, ref_lat)/1000
            PV_ref = 10**(self.a*self.Mw + self.h*self.D + self.d +self.e -np.log10(X_ref+0.0028*10**(0.50*self.Mw)) -0.002*X_ref)
            
            estPGVb600_R = ref_PGVb600 * (PV_est/PV_ref)
            estPGV_R = estPGVb600_R * 1
            
            if self.Relationship_I_and_vel == "conversion":
                Est_Intensity_R = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV_R) -0.082*math.log10(estPGV_R)**2
            if self.Relationship_I_and_vel == "simple":
                upR=PV_est
                downR=PV_ref
                Est_Intensity_R = ref_sindo + 2*math.log10(upR/downR) #mitsuki
                
                
            
           
        ## 表層増幅の補正
        estPGVb600_S = ref_PGVb600 * (1/1)
        estPGV_S = estPGVb600_S * AMP/ref_amp
        if self.Relationship_I_and_vel == "conversion":
            Est_Intensity_S = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV_S) -0.082*math.log10(estPGV_S)**2
        if self.Relationship_I_and_vel == "simple":
            Est_Intensity_S = ref_sindo + 2*math.log10(AMP/ref_amp) #mitsuki
            
        ## 距離減衰と表層増幅の補正
        estPGVb600 = estPGVb600_R #距離補正の方法に応じた基盤面最大速度
        estPGV = estPGVb600 * AMP/ref_amp
        if self.Relationship_I_and_vel == "conversion":
            Est_Intensity_R_S = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV) -0.082*math.log10(estPGV)**2
        if self.Relationship_I_and_vel == "simple":
            Est_Intensity_R_S = ref_sindo + 2*math.log10(AMP/ref_amp) + 2*math.log10(upR/downR)#mitsuki
        
        # ref_name,ref_lon,ref_lat,ref_amp,ref_sindo=0,0,0,0,0,0
        
        return ref_name,ref_lon,ref_lat,ref_amp,ref_sindo,float(Est_Intensity_R_S),float(Est_Intensity_R),float(Est_Intensity_S),upR/downR,AMP/ref_amp
    


    def est_Intensity_B_to_A(self,B_Obs,A_Obs):
        #B→A
        #最も近い観測点から補間
        
        cols = ['経度','緯度','観測点名','震央距離','地盤増幅度',
                '参考観測点名','参考観測点経度','参考観測点緯度','参考観測点地盤増幅度','観測点間距離','参考観測点計測震度',
                '観測計測震度','推計震度','推計震度距離のみ','推計震度地盤増幅のみ']
        estimate_A_from_B = pd.DataFrame(index=[], columns=cols)    
        
        for A in A_Obs.itertuples(): #推計される側
            B_Observation = B_Obs.copy() #推計する側(参照観測点)
            
            name,lon,lat,obs_I,AMP = A[1],A[2],A[3],A[4],A[5]
            # print(name,lon,lat,obs_I,AMP,PGV,PGVb600)
            
            B_Observation=B_Observation[B_Observation['観測点番号'] != name] #A=Bの場合もあるので推計地点を除いた最近傍を探す
    
            ## 最も近い観測点を探す
            _lon_list = [lon]*len(B_Observation)
            _lat_list = [lat]*len(B_Observation)
            azimuth, bkw_azimuth, distance_list = self.grs80.inv(B_Observation['経度'].values.tolist(), B_Observation['緯度'].values.tolist(), _lon_list, _lat_list) 
            min_arg = np.argmin(distance_list)
            min_distance = min(distance_list)
            ref_B = B_Observation.iloc[min_arg]
            
            ## 参照地点のデータ
            ref_name,ref_lon,ref_lat,ref_amp,distance=ref_B['観測点番号'],ref_B['経度'],ref_B['緯度'],ref_B['地盤増幅度'],min_distance
            ref_sindo = ref_B["計測震度"]
            ref_PGVb600 = ref_B["PGV_b600"]
            
            ## 震央距離
            azimuth, bkw_azimuth, hyp_distance = self.grs80.inv(self.epi_lon, self.epi_lat, lon, lat)
            azimuth, bkw_azimuth, ref_hyp_distance = self.grs80.inv(self.epi_lon, self.epi_lat, ref_lon, ref_lat)
            
            ## 震源距離
            epi_distance = pow(hyp_distance**2 + self.depth**2,0.5)
            ref_epi_distance = pow(ref_hyp_distance**2 + self.depth**2,0.5)
            
            
            Est_Intensity_R_S, Est_Intensity_R, Est_Intensity_S = 0, 0, 0
            estPGVb600_R=0
            
            upR=0
            downR=0
                
            ## 距離減衰の補正
            ##　震源距離の逆数
            if self.Distance_attenuation_type ==  "Reciprocal_of_epicenter_distance": 
                estPGVb600_R = ref_PGVb600 * (ref_epi_distance/epi_distance)
                estPGV_R = estPGVb600_R * 1
                if self.Relationship_I_and_vel == "conversion":
                    Est_Intensity_R = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV_R) -0.082*math.log10(estPGV_R)**2
                if self.Relationship_I_and_vel == "simple":
                    upR=ref_epi_distance
                    downR=epi_distance
                    Est_Intensity_R = ref_sindo + 2*math.log10(upR/downR) #mitsuki
            
                
            ##　震源距離の逆数 断層最短距離ver
            elif self.Distance_attenuation_type == "Reciprocal_of_Shortest_fault_distance":
                X_est=self.calculate_Shortest_fault_distance(lon,lat)/1000
                X_ref=self.calculate_Shortest_fault_distance(ref_lon, ref_lat)/1000
                estPGVb600_R = ref_PGVb600 * (X_ref/X_est)
                estPGV_R = estPGVb600_R * 1
                if self.Relationship_I_and_vel == "conversion":
                    Est_Intensity_R = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV_R) -0.082*math.log10(estPGV_R)**2
                if self.Relationship_I_and_vel == "simple":
                    upR=X_ref
                    downR=X_est
                    Est_Intensity_R = ref_sindo + 2*math.log10(upR/downR) #mitsuki
                
            ## 等価震源距離(距離減衰式)
            elif self.Distance_attenuation_type ==  "Using_Equivalent_hypocentral_distance":
                #推計地点
                Xeq = epi_distance/1000
                PV_est = 10**(self.a*self.Mw + self.h*self.D + self.d +self.e -np.log10(Xeq) -0.002*Xeq)
                
                #参照地点
                Xeq = ref_epi_distance/1000
                PV_ref = 10**(self.a*self.Mw + self.h*self.D + self.d +self.e -np.log10(Xeq) -0.002*Xeq)
                
                estPGVb600_R = ref_PGVb600 * (PV_est/PV_ref)
                estPGV_R = estPGVb600_R * 1
                
                if self.Relationship_I_and_vel == "conversion":
                    Est_Intensity_R = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV_R) -0.082*math.log10(estPGV_R)**2
                if self.Relationship_I_and_vel == "simple":
                    upR=PV_est
                    downR=PV_ref
                    Est_Intensity_R = ref_sindo + 2*math.log10(upR/downR) #mitsuki
                
                
            ## 断層最短距離(距離減衰式)
            elif self.Distance_attenuation_type ==  "Shortest_fault_distance":
                #推計地点
                X_est=self.calculate_Shortest_fault_distance(lon,lat)/1000
                PV_est = 10**(self.a*self.Mw + self.h*self.D + self.d +self.e -np.log10(X_est+0.0028*10**(0.50*self.Mw)) -0.002*X_est)
                
                #参照地点
                X_ref=self.calculate_Shortest_fault_distance(ref_lon, ref_lat)/1000
                PV_ref = 10**(self.a*self.Mw + self.h*self.D + self.d +self.e -np.log10(X_ref+0.0028*10**(0.50*self.Mw)) -0.002*X_ref)
                
                estPGVb600_R = ref_PGVb600 * (PV_est/PV_ref)
                estPGV_R = estPGVb600_R * 1
                
                if self.Relationship_I_and_vel == "conversion":
                    Est_Intensity_R = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV_R) -0.082*math.log10(estPGV_R)**2
                if self.Relationship_I_and_vel == "simple":
                    upR=PV_est
                    downR=PV_ref
                    Est_Intensity_R = ref_sindo + 2*math.log10(upR/downR) #mitsuki
                
                
            ## 表層増幅の補正
            estPGVb600_S = ref_PGVb600 * (1/1)
            estPGV_S = estPGVb600_S * AMP/ref_amp
            if self.Relationship_I_and_vel == "conversion":
                Est_Intensity_S = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV_S) -0.082*math.log10(estPGV_S)**2
            if self.Relationship_I_and_vel == "simple":
                Est_Intensity_S = ref_sindo + 2*math.log10(AMP/ref_amp) #mitsuki
            
            ## 距離減衰と表層増幅の補正
            estPGVb600 = estPGVb600_R #距離補正の方法に応じた基盤面最大速度
            estPGV = estPGVb600 * AMP/ref_amp
            if self.Relationship_I_and_vel == "conversion":
                Est_Intensity_R_S = 3.383-0.165*self.Mw + 2.254*math.log10(estPGV) -0.082*math.log10(estPGV)**2
            if self.Relationship_I_and_vel == "simple":
                Est_Intensity_R_S = ref_sindo + 2*math.log10(AMP/ref_amp) + 2*math.log10(upR/downR)#mitsuki
            
                
            print(name,"実際の計測震度",obs_I,"参考観測点震度",ref_sindo,"推計震度",Est_Intensity_R_S,"推計震度_距離のみ",Est_Intensity_R,"推計震度_地盤増幅のみ",Est_Intensity_S,"観測点間距離",distance)
            
            # cols = ['経度','緯度','観測点名','震央距離','地盤増幅度',
            #     '参考観測点名','参考観測点経度','参考観測点緯度','参考観測点地盤増幅度','観測点間距離','参考観測点計測震度',
            #     '観測計測震度','推計震度','推計震度距離のみ','推計震度地盤増幅のみ']
            record = pd.Series([lon,lat,name,hyp_distance,AMP,\
                                ref_name,ref_lon,ref_lat,ref_amp,distance,ref_sindo,\
                                obs_I,Est_Intensity_R_S,Est_Intensity_R,Est_Intensity_S] ,  index=estimate_A_from_B.columns)
            estimate_A_from_B = estimate_A_from_B.append(record,ignore_index=True)
            
        return estimate_A_from_B
    
                  
    def est_Intensity_B_to_A_delaunay(self,B_Obs,A_Obs):
        """
        ドロネー図の頂点から推計震度を計算する
        """
        cols = ['経度','緯度','観測点名','震央距離','地盤増幅度',
                '参考観測点名1','参考観測点経度1','参考観測点緯度1','参考観測点地盤増幅度1','観測点間距離1','参考観測点計測震度1',
                '参考観測点名2','参考観測点経度2','参考観測点緯度2','参考観測点地盤増幅度2','観測点間距離2','参考観測点計測震度2',
                '参考観測点名3','参考観測点経度3','参考観測点緯度3','参考観測点地盤増幅度3','観測点間距離3','参考観測点計測震度3',
                '観測計測震度','推計震度','推計震度距離のみ','推計震度地盤増幅のみ','面積割る長さの二乗','重心からの偏在具合z',
                '距離の補正係数1','地盤増幅の補正係数1',
                '距離の補正係数2','地盤増幅の補正係数2',
                '距離の補正係数3','地盤増幅の補正係数3',
                '距離地盤補正なしの推計震度']
        
        estimate_A_from_B = pd.DataFrame(index=[], columns=cols)    
        grs80 = pyproj.Geod(ellps='GRS80')  # GRS80楕円体
        
        for A in A_Obs.itertuples():
            B_Observation = B_Obs.copy() #推計する側(参照観測点)
            # print(A)
            name,lon,lat,obs_I,AMP = A[1],A[2],A[3],A[4],A[5]
            # print(name,lon,lat,obs_I,AMP,PGV,PGVb600)
            
            azimuth, bkw_azimuth, hyp_distance = grs80.inv(self.epi_lon, self.epi_lat, lon, lat)
            epi_distance = pow(hyp_distance**2 + self.depth**2,0.5)
            
            B_Observation=B_Observation[B_Observation['観測点番号'] != name] #推計地点を除く
    
            #推計に使う参照観測点を探す
            
            #緯度経度を平面直角座標（XY座標）に変換する
            lon_list=B_Observation['経度'].values.tolist()
            lat_list=B_Observation['緯度'].values.tolist()
            EPSG4612 = pyproj.Proj("+init=EPSG:4612")
            EPSG2451 = pyproj.Proj("+init=EPSG:2451")   
            y,x = pyproj.transform(EPSG4612, EPSG2451,lon_list,lat_list)
            x_1d,y_1d = np.array(x),np.array(y)
            # pts = B_Observation[['緯度', '経度']].values
            pts=np.stack([x_1d, y_1d],1)
    
            #検索する観測点の緯度経度を平面直角座標（XY座標）に変換する
            lon_y,lat_x = pyproj.transform(EPSG4612, EPSG2451,lon,lat)
            # print(lon_y,lat_x)
            
            #ドロネー図作成
            tri = Delaunay(pts)
            # delaunay_plot_2d(tri)
            
            inside = -1
            # p1 = Point(lat_x,lon_y)#判断したいポイント  
            for i,vertices in enumerate(pts[tri.simplices]):#三角形の頂点
                result = self.inside_or_outside_judgment(vertices,[lat_x,lon_y])
                if result == True:
                    inside = i
                    break
            
            ## ここまでが参照観測点選び
            
            ## 震度の推計
            num1 = tri.simplices[inside][0]
            # print(num)
            ref_B1 = B_Observation.iloc[num1]
            ref1_lon,ref1_lat=ref_B1['経度'],ref_B1['緯度']
            azimuth, bkw_azimuth, ref1_distance = grs80.inv(ref1_lon, ref1_lat, lon, lat) 
            ref1_name,ref1_lon,ref1_lat,ref1_amp,ref1_sindo,Est1_Intensity_R_S,Est1_Intensity_R,Est1_Intensity_S\
                ,Distance_correction_coefficient1,AMP_correction_coefficient1\
                    =self.get_ref_data(ref_B1,epi_distance,AMP, lon, lat)
            # print(name,"実際の計測震度",obs_I,"参考観測点震度",ref1_sindo,"推計震度",Est1_Intensity_R_S,"推計震度_距離のみ",Est1_Intensity_R,"推計震度_地盤増幅のみ",Est1_Intensity_S,"観測点間距離",ref1_distance)
            
            num2 = tri.simplices[inside][1]
            # print(num)
            ref_B2 = B_Observation.iloc[num2]
            ref2_lon,ref2_lat=ref_B2['経度'],ref_B2['緯度']
            azimuth, bkw_azimuth, ref2_distance = grs80.inv(ref2_lon, ref2_lat, lon, lat) 
            ref2_name,ref2_lon,ref2_lat,ref2_amp,ref2_sindo,Est2_Intensity_R_S,Est2_Intensity_R,Est2_Intensity_S\
                ,Distance_correction_coefficient2,AMP_correction_coefficient2\
                    =self.get_ref_data(ref_B2,epi_distance,AMP, lon, lat)
            # print(name,"実際の計測震度",obs_I,"参考観測点震度",ref2_sindo,"推計震度",Est2_Intensity_R_S,"推計震度_距離のみ",Est2_Intensity_R,"推計震度_地盤増幅のみ",Est2_Intensity_S,"観測点間距離",ref2_distance)
            
            num3 = tri.simplices[inside][2]
            # print(num)
            ref_B3 = B_Observation.iloc[num3]
            ref3_lon,ref3_lat=ref_B3['経度'],ref_B3['緯度']
            azimuth, bkw_azimuth, ref3_distance = grs80.inv(ref3_lon, ref3_lat, lon, lat) 
            ref3_name,ref3_lon,ref3_lat,ref3_amp,ref3_sindo,Est3_Intensity_R_S,Est3_Intensity_R,Est3_Intensity_S\
                ,Distance_correction_coefficient3,AMP_correction_coefficient3\
                    =self.get_ref_data(ref_B3,epi_distance,AMP, lon, lat)
            # print(name,"実際の計測震度",obs_I,"参考観測点震度",ref3_sindo,"推計震度",Est3_Intensity_R_S,"推計震度_距離のみ",Est3_Intensity_R,"推計震度_地盤増幅のみ",Est3_Intensity_S,"観測点間距離",ref3_distance)
            
            
            
            ## 観測点間距離で重みづけ平均
            if ref1_distance!=0 and ref2_distance!=0 and ref3_distance!=0:
                
                Est_Intensity_R_S\
                    =(Est1_Intensity_R_S*(1/ref1_distance)+Est2_Intensity_R_S*(1/ref2_distance)+Est3_Intensity_R_S*(1/ref3_distance))\
                        /(1/ref1_distance + 1/ref2_distance + 1/ref3_distance)
                
                Est_Intensity_R\
                    =(Est1_Intensity_R*(1/ref1_distance)+Est2_Intensity_R*(1/ref2_distance)+Est3_Intensity_R*(1/ref3_distance))\
                        /(1/ref1_distance + 1/ref2_distance + 1/ref3_distance)
                
                Est_Intensity_S\
                    =(Est1_Intensity_S*(1/ref1_distance)+Est2_Intensity_S*(1/ref2_distance)+Est3_Intensity_S*(1/ref3_distance))\
                        /(1/ref1_distance + 1/ref2_distance + 1/ref3_distance)
                #補正をしない場合
                Intensity_without_correction\
                    =(ref1_sindo*(1/ref1_distance)+ref2_sindo*(1/ref2_distance)+ref3_sindo*(1/ref3_distance))\
                        /(1/ref1_distance + 1/ref2_distance + 1/ref3_distance)
            
            else:
                Est_Intensity_R_S,Est_Intensity_R,Est_Intensity_S=0,0,0
                Intensity_without_correction=0
                print("推計失敗")
            
            
            
            ## ドロネー図の形状と推計精度の関係
            
            ## ドロネー図の面積と長さの二乗の比
            x1,y1=pts[num1][0],pts[num1][1]
            x2,y2=pts[num2][0],pts[num2][1]
            x3,y3=pts[num3][0],pts[num3][1]
            S=abs(x1*y2+x2*y3+x3*y1-y1*x2-y2*x3-y3*x1)/2#面積
            L=pow((x1-x2)**2+(y1-y2)**2,0.5)+pow((x2-x3)**2+(y2-y3)**2,0.5)+pow((x3-x1)**2+(y3-y1)**2,0.5) #長さの二乗
            L2=L*L
            SperL2=S/L2
            # print("面積割る長さの二乗",SperL2)
            
            ## 推計地点の重心からの偏在具合
            Qx,Qy=lat_x,lon_y
            α1 = sp.Symbol('α1')
            α2 = sp.Symbol('α2')
            α3 = sp.Symbol('α3')
            equation1 = α1*x1 + α2*x2 + α3*x3 - Qx
            equation2 = α1*y1 + α2*y2 + α3*y3 - Qy
            equation3 = α1 + α2 + α3 -1
            a=sp.solve([equation1, equation2, equation3])
            a1=a[α1]; a2=a[α2]; a3=a[α3]
            z = a1 + a2 * np.exp(2*np.pi*1j/3) + a3 * np.exp(4*np.pi*1j/3)
            
            
            
            print(name,"推計",Est_Intensity_R_S,"実際",obs_I,"偏在具合",abs(z),"\n")
    
            record = pd.Series([lon,lat,name,hyp_distance,AMP,\
                                ref1_name,ref1_lon,ref1_lat,ref1_amp,ref1_distance,ref1_sindo,\
                                ref2_name,ref2_lon,ref2_lat,ref2_amp,ref2_distance,ref2_sindo,\
                                ref3_name,ref3_lon,ref3_lat,ref3_amp,ref3_distance,ref3_sindo,\
                                obs_I,Est_Intensity_R_S,Est_Intensity_R,Est_Intensity_S,SperL2,abs(z),\
                                Distance_correction_coefficient1,AMP_correction_coefficient1,\
                                Distance_correction_coefficient2,AMP_correction_coefficient2,\
                                Distance_correction_coefficient3,AMP_correction_coefficient3,\
                                Intensity_without_correction],index=estimate_A_from_B.columns)        
            estimate_A_from_B = estimate_A_from_B.append(record,ignore_index=True)     
            
            
        return estimate_A_from_B
    
    def est_Intensity_B_to_A_voronoi(self,B_Obs,A_Obs):
        """
        ボロノイ図から推計震度を計算する
        """
        cols = ['経度','緯度','観測点名','震央距離','地盤増幅度',
                '観測計測震度','推計震度','推計震度距離のみ','推計震度地盤増幅のみ','参照観測点数','平均距離','最大距離','最小距離','ボロノイ面積',
                '参照地点','距離地盤補正なしの推計震度',
                '距離の補正係数平均','地盤増幅の補正係数平均']
    
        estimate_A_from_B = pd.DataFrame(index=[], columns=cols)    
        grs80 = pyproj.Geod(ellps='GRS80')  # GRS80楕円体
        for A in A_Obs.itertuples():
            B_Observation = B_Obs.copy() #推計する側(参照観測点)
            # print(A)
            name,lon,lat,obs_I,AMP = A[1],A[2],A[3],A[4],A[5]
            # print(name,lon,lat,obs_I,AMP,PGV,PGVb600)

            
            azimuth, bkw_azimuth, hyp_distance = grs80.inv(self.epi_lon, self.epi_lat, lon, lat)
            epi_distance = pow(hyp_distance**2 + self.depth**2,0.5)
            
            B_Observation=B_Observation[B_Observation['観測点番号'] != name] #推計地点を除く
    
            #推計に使う観測点を探す
            
            #緯度経度を平面直角座標（XY座標）に変換する
            lon_list=B_Observation['経度'].values.tolist()
            lat_list=B_Observation['緯度'].values.tolist()
            EPSG4612 = pyproj.Proj("+init=EPSG:4612")
            EPSG2451 = pyproj.Proj("+init=EPSG:2451")   
            y,x = pyproj.transform(EPSG4612, EPSG2451,lon_list,lat_list)
            x_1d,y_1d = np.array(x),np.array(y)
            # pts = B_Observation[['緯度', '経度']].values
            pts=np.stack([x_1d, y_1d],1)
    
            #検索する観測点の緯度経度を平面直角座標（XY座標）に変換する
            lon_y,lat_x = pyproj.transform(EPSG4612, EPSG2451,lon,lat)
            # print(lon_y,lat_x)
            
            #ボロノイ図作成
            vor = Voronoi(pts)
            inside = -1
            voronoi_area=0
            
            
            # p1 = Point(lat_x,lon_y)#判断したいポイントがどの領域に含まれるか 
            for i,region in enumerate(vor.regions):#[r for r in vor.regions if -1 not in r and r]):#ボロノイ領域ごとに内外判定
                if -1 not in region and region:
                    result = self.inside_or_outside_judgment(vor.vertices[region],[lat_x,lon_y])
                    if result == True:
                        poly = Polygon(*vor.vertices[region])
                        inside = i
                        voronoi_area=float(poly.area)
                        break
            
            
            if inside != -1:
                #辺を共有する領域の構成点から推定する vor.regions[inside]とregion
                share_points_list = []
                for i,region in enumerate(vor.regions):#[r for r in vor.regions if -1 not in r and r]):#ボロノイ領域ごとに内外判定
                    if -1 not in region and region:
                        share_point_num=0
                        for vertex in vor.regions[inside]:#推計地点が含まれる領域の頂点番号
                            share_point_num += 1 if vertex in region else 0
                        if share_point_num == 2:
                            share_points_list.append(i) 
                        
                share_points_list.append(inside)
                
                ref_name_list,ref_lon_list,ref_lat_list,ref_amp_list,ref_distance_list,ref_sindo_list=[],[],[],[],[],[]
                Est_Intensity_R_S_list,Est_Intensity_R_list,Est_Intensity_S_list=[],[],[]
                
                code_list=[]
                Distance_correction_coefficient_list,AMP_correction_coefficient_list=[],[]
                for num in share_points_list:
                    if len(vor.vertices[vor.regions[num]]) >= 2:
                        B_num=-1
                        for n,C in enumerate(B_Observation.itertuples()):#vor と B_Observationの順番が変わるから再び内外判定をする
                            result = self.inside_or_outside_judgment(vor.vertices[vor.regions[num]],[x[n],y[n]])
                            if result == True:
                                B_num = n
                                break
                        if (B_num == -1) or (B_num >= len(B_Observation)):
                            print("B_Observationの範囲外")
                            with open("error.txt", mode='w') as f:
                                f.write("{},{}".format(B_num,num))
                        else:
                            ref_B = B_Observation.iloc[B_num]
                            ref_lon,ref_lat=ref_B['経度'],ref_B['緯度']
                            azimuth, bkw_azimuth, ref_distance = grs80.inv(ref_lon, ref_lat, lon, lat) 
                            ref_name,ref_lon,ref_lat,ref_amp,ref_sindo,Est_Intensity_R_S,Est_Intensity_R,Est_Intensity_S\
                                 ,Distance_correction_coefficient,AMP_correction_coefficient\
                                     =self.get_ref_data(ref_B,epi_distance,AMP, lon, lat)
                            # print(name,"実際の計測震度",obs_I,"参考観測点震度",ref_sindo,"推計震度",Est_Intensity_R_S,"推計震度_距離のみ",Est_Intensity_R,"推計震度_地盤増幅のみ",Est_Intensity_S,"観測点間距離",ref_distance)
                            ref_name_list.append(ref_name)
                            ref_lon_list.append(ref_lon)
                            ref_lat_list.append(ref_lat)
                            ref_amp_list.append(ref_amp)
                            ref_distance_list.append(ref_distance)
                            ref_sindo_list.append(float(ref_sindo))
                            Est_Intensity_R_S_list.append(Est_Intensity_R_S)
                            Est_Intensity_R_list.append(Est_Intensity_R)
                            Est_Intensity_S_list.append(Est_Intensity_S)
                            
                            code_list.append(ref_B['観測点番号'])
                            Distance_correction_coefficient_list.append(Distance_correction_coefficient)
                            AMP_correction_coefficient_list.append(AMP_correction_coefficient)
                            # sys.exit()
    
                
                ref_distance_list = np.array(ref_distance_list)
                ref_per_distance_list = np.array(1/ref_distance_list)
                
                Est_Intensity_R_S_list = np.array(Est_Intensity_R_S_list)
                Est_Intensity_S_list = np.array(Est_Intensity_S_list)
                Est_Intensity_R_list = np.array(Est_Intensity_R_list)
                
                numerator=Est_Intensity_R_S_list*ref_per_distance_list
                numerator_R=Est_Intensity_R_list*ref_per_distance_list
                numerator_S=Est_Intensity_S_list*ref_per_distance_list
                
                # #観測点間距離で重みづけ平均
                Est_Intensity_R_S = numerator.sum()/ref_per_distance_list.sum()
                Est_Intensity_R = numerator_R.sum()/ref_per_distance_list.sum()
                Est_Intensity_S = numerator_S.sum()/ref_per_distance_list.sum()
                
                #補正をしない場合
                ref_sindo_list= np.array(ref_sindo_list)
                numerator_without_correction = ref_sindo_list*ref_per_distance_list
                Intensity_without_correction = numerator_without_correction.sum()/ref_per_distance_list.sum()
                
                #補正係数
                Distance_correction_coefficient_ave = sum(Distance_correction_coefficient_list) / len(Distance_correction_coefficient_list)
                AMP_correction_coefficient_ave = sum(AMP_correction_coefficient_list) / len(AMP_correction_coefficient_list)
                
                
                print(name,"実際の計測震度",obs_I,"推計震度",Est_Intensity_R_S,"推計震度_距離のみ",Est_Intensity_R,"推計震度_地盤増幅のみ",Est_Intensity_S,"参照観測点数",len(share_points_list)-1,"平均観測点間距離",ref_distance_list.mean())
                print("距離地盤補正なしの推計震度",Intensity_without_correction)
                print("\n")
                record = pd.Series([lon,lat,name,hyp_distance,AMP,\
                                    obs_I,Est_Intensity_R_S,Est_Intensity_R,Est_Intensity_S,len(share_points_list)-1,ref_distance_list.mean(),ref_distance_list.max(),ref_distance_list.min(),voronoi_area,\
                                    code_list,Intensity_without_correction,
                                    Distance_correction_coefficient_ave,AMP_correction_coefficient_ave],index=estimate_A_from_B.columns)        
                estimate_A_from_B = estimate_A_from_B.append(record,ignore_index=True)
                estimate_A_from_B.to_csv('test.csv',encoding="shift_jis",index=False)
                
        return estimate_A_from_B
    
            
