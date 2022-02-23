# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 22:15:52 2021

@author: sugiy

https://www.data.jma.go.jp/svd/eqev/data/bulletin/shindo.html
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial import Delaunay, delaunay_plot_2d, Voronoi, voronoi_plot_2d
from sympy.geometry import Point, Polygon
import sympy as sp
import math
import time
import folium
import glob
import pyproj
import sys
import os
import warnings
warnings.filterwarnings('ignore')

from config import Config
    

def grid250m(lon, lat):
        # https://www.stat.go.jp/data/mesh/pdf/gaiyo1.pdf
        p=int(lat*60/40)
        a=(lat*60)%40
        
        q=int(a/5)
        b=a%5
        
        r=int(b*60/30)
        c=(b*60)%30
        
        s=int(c/15)
        d=c%15
        
        t=int(d/7.5)
        f=d%7.5
        
        # あ=int(f/3.75)
        
        u=int(lon-100)
        f=lon-100-u
        
        v=int(f*60/7.5)
        g=f*60%7.5
        
        w=int(g*60/45)
        h=g*60%45
        
        x=int(h/22.5)
        i=h%22.5
        
        y=int(i/11.25)
        # j=i%11.25
        # い=int(j/5.625)
        
        m=s*2+x+1
        n=t*2+y+1
        # o=あ*2+い+1
        return( str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(m)+str(n))
        

if __name__ == '__main__':
    
    ## 設定ファイル
    conf=Config()
    grs80 = pyproj.Geod(ellps='GRS80') 
    
    MakeData = True
    MakeData = False
    
    if MakeData:
        ## 年ごとに処理をする
        for target_year in conf.target_year:
            
            cols = ['発生年','発生年月日時分','Mj','深さ[km]','最大震度(階級)','最大震度(計測震度)','最短震央距離[km]','最短震源距離[km]','観測記録数','陸海']
            Sindo_data = pd.DataFrame(index=[], columns=cols)    
        
            print("target_year",target_year)
            EQdata_path = '../Catalog/i{}.dat'.format(target_year)
            
            eq_key=""
            eq_key_list=[] #発生日時で地震を管理する
            
            ## 年単位のファイルの読み込み
            with open(EQdata_path) as f:
                
                ## 地震ごとに個別のファイルに書き出す
                for i,line in enumerate(f):
                    if line[0]=="A" or line[0]=="B":
                        eq_key=line[1:16]
                        eq_key_list.append(eq_key)
                        try: #同一ファイルがある場合(2回目の実行時など)に既存ファイルを削除する
                            os.remove('../Intensity_data_multipleEQ/{}.txt'.format(eq_key))
                        except:
                            None
                    path_w = '../Intensity_data_multipleEQ/{}.txt'.format(eq_key)
                    with open(path_w, mode='a') as f:
                        f.write(line)
                
            for eq_key in eq_key_list:
                ## 地震ごとの処理
                try:
                    input_path='../Intensity_data_multipleEQ/{}.txt'.format(eq_key)
    
                    ## 震度データをDFに変換
                    code_list, lon_list, lat_list, I_list,I_class_list, amp_list, organ_list = [],[],[],[],[],[],[]
                    input_Observation=None
                    i=0
                    epi_lat,epi_lon,depth,Mj = 0,0,0,0
                    I_class_max = 0
                    output_path ,output_flag = "",""
                    land_or_sea = 1
                    
                    Make_DF = True
                    print()
                    print(f"eq_key:{eq_key}")
                    ## 震度データ読み込み(単一地震ファイル)
                    with open(input_path) as f:
                        for i,line in enumerate(f):
                            if i==0:
                                
                                #震源情報をヘッダーから読み取る
                                try:
                                    print("line",line)
                                    print("line.split()",line.split())
                                    #https://www.data.jma.go.jp/svd/eqev/data/bulletin/data/shindo/format_j.pdf
                                    
                                    # eqi_lat_int=line.split()[2][:2]
                                    # eqi_lat_dec1=line.split()[2][2:4]
                                    # eqi_lat_dec2=line.split()[2][4:] if line.split()[2][4:] != "" else 0
                                    # # print(eqi_lat_int,eqi_lat_dec1,eqi_lat_dec2)
                                    # epi_lat=float(eqi_lat_int) + (float(eqi_lat_dec1)+float(eqi_lat_dec2)/100)/60
                                    
                                    eqi_lat_int=line[22-1:24].replace(" ", "")
                                    eqi_lat_dec1=line[25-1:26] if line[25-1:26] != "  " else 0
                                    eqi_lat_dec2=line[27-1:28] if line[27-1:28] != "  " else 0
                                    epi_lat=float(eqi_lat_int) + (float(eqi_lat_dec1)+float(eqi_lat_dec2)/100)/60
                                    
                                    # eqi_lon_int=line.split()[4][:3]
                                    # eqi_lon_dec1=line.split()[4][3:5]
                                    # eqi_lon_dec2=line.split()[4][5:] if line.split()[4][5:] != "" else 0
                                    # # print(eqi_lon_int,eqi_lon_dec1,eqi_lon_dec2)
                                    # epi_lon=float(eqi_lon_int) + (float(eqi_lon_dec1)+float(eqi_lon_dec2)/100)/60
                                    
                                    eqi_lon_int=line[33-1:36]
                                    eqi_lon_dec1=line[37-1:38] if line[37-1:38] != "  " else 0
                                    eqi_lon_dec2=line[39-1:40] if line[39-1:40] != "  " else 0
                                    epi_lon=float(eqi_lon_int) + (float(eqi_lon_dec1)+float(eqi_lon_dec2)/100)/60
                                    
                                    
                                    # depth = float(line.split()[6][:5])/1000 if line.split()[6][:5] != "" else 9999
                                    # Mj = float(line.split()[6][10:12])/10 if line.split()[6][10:12] != "" else 9999
                                    
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
                                    
                                    
                                    print("震源情報",epi_lon,epi_lat,depth,Mj,I_class_max)
                                    
                                    
                                    
                                #ヘッダーを読み込めない場合    
                                except Exception as e:
                                    print("震源情報が読み込めませんでした\n")
                                    print("message:{0}".format(str(e)))
                                    exc_type, exc_obj, exc_tb = sys.exc_info()
                                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                                    print(exc_type, fname, exc_tb.tb_lineno)
                                    print()
                                    Make_DF = False
                                    
                                    time.sleep(10)
                                    
                                    path_w = '../output/震源読み込めず_{}.txt'.format(eq_key)
                                    with open(path_w, mode='a') as f:
                                        f.write(line)
                
                                    break
                                
                                #陸域の地震か判定
                                mesh = grid250m(epi_lon,epi_lat)
                                index=conf.JSHIS_ARV[conf.JSHIS_ARV[0]==int(mesh)].index
                                grep_data = conf.JSHIS_ARV.loc[index]
                                # print(epi_lon,epi_lat)
                                # print(mesh)
                                # print(index)
                                AVS=0
                                if index > 1:
                                    AVS=float(grep_data[2])
                                if AVS == 0:
                                    land_or_sea = 0
                                
                                continue
                            
                            # print(line.split())
                            
                            code=line.split()[0]#.replace(" ","")
                            I_class=line.split()[2]
                            I=float(line.split()[3])/10 if line.split()[3] != "//" else 9999
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
                            
                            
                            
                                
                            ## 各観測点の記録をlistに追加
                            code_list.append(code);
                            lon_list.append(lon);
                            lat_list.append(lat);
                            I_list.append(I);
                            I_class_list.append(I_class);
                            amp_list.append(amp);
                            organ_list.append(organ)
                            
                            
                            
                            
                    ## 地震ごとの集計
                    code_df, lon_df, lat_df, I_df, I_class_df, amp_df, organ_df = pd.Series(code_list),pd.Series(lon_list),pd.Series(lat_list),pd.Series(I_list),pd.Series(I_class_list),pd.Series(amp_list),pd.Series(organ_list)
                    Obs = pd.concat([code_df, lon_df, lat_df, I_class_df, I_df, amp_df, organ_df], axis=1)
                    Obs.columns = ['観測点番号','経度','緯度','震度階級','計測震度',"地盤増幅度","機関"]
                    
                    print(Obs)
                    
                    
                    
                    # 最大震度(計測震度)
                    I_max=Obs["計測震度"].max()
                    
                    # 最大震度を記録した観測点の震央距離
                    # 最大震度を記録した観測点の震源距離　epi_lon,epi_lat,depth,Mj
                    
                    
                    # 最短震央距離
                    _lon_list = [epi_lon]*len(Obs)
                    _lat_list = [epi_lat]*len(Obs)
                    azimuth, bkw_azimuth, distance_list = grs80.inv(Obs['経度'].values.tolist(), Obs['緯度'].values.tolist(), _lon_list, _lat_list) 
                    min_distance = round(min(distance_list) /1000,2)
                    
                    # 最短震源距離
                    min_epi_distance = round(pow(min_distance**2 + depth**2,0.5),2)
                    
                    min_arg = np.argmin(distance_list)
                    max_Obs = Obs.iloc[min_arg]
                    
                    # 観測記録数
                    Obs_len = len(Obs)
                    
                    print(f"Mj{Mj} 深さ{depth}km 最大震度(階級){I_class_max} 最大震度(計測震度){I_max} 最短震央距離{min_distance}km 最短震源距離{min_epi_distance}km 観測記録数{Obs_len} 陸海{land_or_sea}")
                    
                    #cols = ['発生年','発生年月日時分','Mj','深さ[km]','最大震度(階級)','最大震度(計測震度)','最短震央距離[km]','最短震源距離[km]','観測記録数','陸海']
                    record = pd.Series([target_year,
                                        eq_key,
                                        Mj,
                                        depth,
                                        I_class_max,
                                        I_max,
                                        min_distance,
                                        min_epi_distance,
                                        Obs_len,
                                        land_or_sea] ,  
                                       index=Sindo_data.columns)
                    Sindo_data = Sindo_data.append(record,ignore_index=True)
    
                except Exception as e:
                    print("message:{0}".format(str(e)))
       
            Sindo_data.to_csv(f"../output/震度データ_{target_year}.csv",encoding="shift_jis",index=False)
        
        
    
    MakeGraph = True
    # MakeGraph = False
    if MakeGraph:
        
        # 一つのdfに結合
        df_sindodata = pd.DataFrame()
        for target_year in conf.target_year:
            df_temp = pd.read_csv(f"../output/震度データ_{target_year}.csv",encoding="shift_jis")
            print(df_temp)
            # http://sinhrks.hatenablog.com/entry/2015/01/28/073327
            df_sindodata = pd.concat([df_sindodata, df_temp])
        
        # 震度5を震度5弱へ 震度6を震度6弱へ
        df_sindodata=df_sindodata.replace('震度5', '震度5弱')
        df_sindodata=df_sindodata.replace('震度6', '震度6弱')
        
        
        # 集計dfにデータを集める
        df_aggregate = pd.DataFrame()
        sindo_list = ["震度1","震度2","震度3","震度4","震度5弱","震度5強","震度6弱","震度6強","震度7"]
        
        step=10
        age=range(1980,2020,step)#[1980,1990,2000,2010]
        
        for start_year in age:
            end_year=start_year+step-1

            print(f"{start_year}年~{end_year}年")
            df = df_sindodata.copy()
            df = df[df["陸海"]==1]
            df = df[df["発生年"]>=start_year]
            df = df[df["発生年"]<=end_year]
            sindo_count_list=[]
            for sindo in sindo_list :
                df_1 = df[df["最大震度(階級)"]==sindo]
                count=len(df_1)   #-1 if len(df_1)-1 != 0 else 0
                # print(f"{sindo}：{count}回")
                sindo_count_list.append(count)
            
            s = pd.DataFrame(sindo_count_list, index=sindo_list, columns=[f"{start_year}年~{end_year}年"])
            print(s)
            ##震度観測記録のヒストグラム
            fig = plt.figure(figsize=(10,4),dpi=300,facecolor='w',edgecolor='k')
            ax = fig.add_subplot(1,1,1)
            plt.rcParams["font.family"] = "MS Gothic"
            plt.rcParams['font.size'] = 12
            label=sindo_list
            data = sindo_count_list
            ax.bar(height=data,x=label)
            ax.set_ylim(0,30000)
            # ax.set_title('first histogram $\mu=100,\ \sigma=15$')
            # ax.set_xlabel('x')
            ax.set_ylabel('頻度')
            fig.show()
            df_aggregate = pd.concat([df_aggregate, s], axis=1)
        
        print()
        print("震度観測記録数")
        print(df_aggregate)
        
        
        ## 観測記録数ヒストグラム年代重ね書き
        # sindo_list = ["震度1","震度2","震度3","震度4","震度5弱","震度5強","震度6弱","震度6強","震度7"]
        sindo_list = ["震度4","震度5弱","震度5強","震度6弱","震度6強","震度7"]
        width=0.1
        fig = plt.figure(figsize=(10,4),dpi=300,facecolor='w',edgecolor='k')
        ax = fig.add_subplot(1,1,1)
        plt.rcParams["font.family"] = "MS Gothic"
        plt.rcParams['font.size'] = 12
        x_label=np.arange(len(df_aggregate.columns)) #,x=df_aggregate.columns
        ax.set_xticks(x_label+width*(len(sindo_list)/2)) #1980年~1989年  1990年~1999年  2000年~2009年  2010年~2019年の位置指定
        ax.set_xticklabels(labels=df_aggregate.columns)
        ax.set_ylabel('震度別　頻度')
        ax.grid(color='black', linestyle='dotted', linewidth=1, axis = "y")
        for i,target_sindo in enumerate(sindo_list):
            A=df_aggregate.loc[target_sindo]
            ax.bar(x=x_label+width*i,height=A.values.tolist(),label=target_sindo, align='edge', width=width, linewidth=0.5, edgecolor="black")
        ax.legend(sindo_list)#df_aggregate.index columns
        # ax.set_ylim(0,20000)
        ax.set_title('震度階級別観測記録数のヒストグラム (震央が陸域の地震)',y = -0.2)
        # ax.set_xlabel('x')
        fig.savefig("../output/fig/震度階級別観測記録数のヒストグラム (震央が陸域の地震).png")
        fig.show()
        
        
        
        
        ## 震度観測点数の推移
        # 000000000000 は記録がないが観測終了
        start_year = 1980
        end_year = 2018
        
        year_list=[]
        station_count_list=[]
        station_count_list_JMA=[]
        station_count_list_KNET=[]
        station_count_list_LG=[]
        
        for year in range(start_year,end_year+1,1):
            df = conf.matchingtable.copy()
            df = df[df["観測終了年月日時分"]!=0]
            #年末時点での数
            year_str = int(str(year+1)+"00000000")
            df = df[df["観測開始年月日時分"]<=year_str]
            year_end = int(str(year+1)+"00000000")
            df = df[df["観測終了年月日時分"]>=year_end]
            
            year_list.append(str(year)+"年")
            
            station_count=len(df)
            station_count_list.append(station_count)
            
            station_count_JMA=len(df[df["機関"]=="気象庁"])
            station_count_list_JMA.append(station_count_JMA)
            station_count_KNET=len(df[df["機関"]=="防災科研"])
            station_count_list_KNET.append(station_count_KNET)
            station_count_LG=len(df[df["機関"]=="自治体"])
            station_count_list_LG.append(station_count_LG)
            
            print(year,station_count,station_count_JMA,station_count_KNET,station_count_LG)
            
        # plt.style.use('fivethirtyeight')
        # plt.style.use('ggplot')
        # plt.style.use('seaborn')
        # plt.style.use('grayscale')
        
        organ_list=["気象庁","自治体","防災科研"]
        
        fig = plt.figure(figsize=(10,4),dpi=300,facecolor='w',edgecolor='k')
        ax = fig.add_subplot(1,1,1)
        plt.rcParams["font.family"] = "MS Gothic"
        plt.rcParams['font.size'] = 12
        # x_label=np.arange(len(year_list))-1 #,x=df_aggregate.columns
        # ax.set_xticks(x_label+(len(year_list)/2)) #1980年~1989年  1990年~1999年  2000年~2009年  2010年~2019年の位置指定
        fig.autofmt_xdate(rotation=80)
        # ax.set_xticklabels(labels=year_list)
        ax.set_ylabel('震度観測点数')
        ax.grid(color='black', linestyle='dotted', linewidth=1, axis = "y")
        
        ax.bar(x=year_list,height=station_count_list_JMA,label="気象庁", linewidth=0.5, edgecolor="black", align='center')
        ax.bar(x=year_list,height=station_count_list_LG,label="自治体",  linewidth=0.5, edgecolor="black", bottom=station_count_list_JMA, align='center')
        bottom_temp=[]
        for a,b in zip(station_count_list_JMA,station_count_list_LG):
            bottom_temp.append(a+b)
        ax.bar(x=year_list,height=station_count_list_KNET,label="防災科研",  linewidth=0.5, edgecolor="black", bottom=bottom_temp, align='center')
    
        ax.legend(organ_list)#df_aggregate.index columns
        fig.savefig("../output/fig/震度観測点数の推移域.png")
        fig.show()
        
        fig = plt.figure(figsize=(10,4),dpi=300,facecolor='w',edgecolor='k')
        ax = fig.add_subplot(1,1,1)
        ax.bar(x=year_list,height=station_count_list,label="観測点",  linewidth=0.5, edgecolor="black", align='center')
        fig.show()
        
        
        # df = df[df["最大震度(計測震度)"]<10]
        # fig = plt.figure()
        # ax = fig.add_subplot(1,1,1)
        # x=df["最大震度(計測震度)"]
        # ax.hist(x, bins=50)
        # ax.set_title('first histogram $\mu=100,\ \sigma=15$')
        # ax.set_xlabel('x')
        # ax.set_ylabel('freq')
        # fig.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
                            