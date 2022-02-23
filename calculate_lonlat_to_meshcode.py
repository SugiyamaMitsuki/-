# -*- coding: utf-8 -*-
import math as mt
import pandas as pd
import sys

# https://qiita.com/yuusei/items/549402a80efd7e7192ef


class LatLon2Grid():
    def __init__(self):
        super(LatLon2Grid, self).__init__()
        
        
        
        
    def grid1st(self,lon, lat):  # 1次メッシュ(4桁) 分割なし
        return int(mt.floor(lat*1.5)) * 100 + int(mt.floor(lon-100))

    def grid2nd(self,lon, lat):  # 2次メッシュ(6桁) 8分割
       return (int(mt.floor(lat*12       / 8))   * 10000 + int(mt.floor((lon-100)*8         / 8))  * 100   +   
               int(mt.floor(lat*12 %  8     ))   * 10    + int(mt.floor((lon-100)*8))  %  8               )  

    def grid3rd(self,lon, lat):  # 3次メッシュ(8桁) 8分割x10分割=80分割
        return (int(mt.floor(lat*120      / 80)) * 1000000 + int(mt.floor((lon-100))             ) * 10000 +  # 1次メッシュ
                int(mt.floor(lat*120 % 80 / 10)) * 1000    + int(mt.floor((lon-100)*80 % 80 / 10)) * 100 +    # 2次メッシュ
                int(mt.floor(lat*120 % 10))      * 10      + int(mt.floor((lon-100)*80)) % 10               ) 

    def grid4th(self,lon, lat):  # 4次メッシュ(9桁) 8分割x10分割x2分割=160分割
        return (int(mt.floor(lat*240       / 160)) * 10000000 + int(mt.floor((lon-100)*160       / 160)) * 100000 +    # 1次メッシュ
                int(mt.floor(lat*240 % 160 / 20))  * 10000    + int(mt.floor((lon-100)*160 % 160 / 20))  * 1000   +    # 2次メッシュ
                int(mt.floor(lat*240 % 20  / 2))   * 100      + int(mt.floor((lon-100)*160 % 20  / 2))   * 10     +    # 3次メッシュ
                int(mt.floor(lat*240)) % 2         * 2        + int(mt.floor((lon-100)*160)) % 2                  + 1) # 4次メッシュ
    
    def grid250m(self,lon, lat):
        # ４分の１地域メッシュ・コード
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
    
    def grid_11_digits(self,lon, lat):
        # https://www.stat.go.jp/data/mesh/pdf/gaiyo1.pdf
        
        ## lat
        p=int(lat*60/40)
        a=(lat*60)%40
        
        q=int(a/5)
        b=a%5
        
        r=int(b*60/30)
        c=(b*60)%30
        
        s=int(c/15)
        d=c%15
        
        t=int(d/7.5)
        tt=d%7.5
        
        a1=int(tt/3.75)
        aa1=tt%3.75
        
        b1=int(aa1/1.875)
        bb1=aa1%1.875
        
        c1=int(bb1/0.9375)
        # cc1=bb1%0.9375
        
        
        
        ## lon
        u=int(lon-100)
        f=lon-100-u
        
        v=int(f*60/7.5)
        g=f*60%7.5
        
        w=int(g*60/45)
        h=g*60%45
        
        x=int(h/22.5)
        i=h%22.5
        
        y=int(i/11.25)
        yy=i%11.25
        
        a2=int(yy/5.625)
        aa2=yy%5.625
        
        b2=int(aa2/2.8125)
        bb2=aa2%2.8125
        
        c2=int(bb2/1.40625)
        # cc2=bb2%1.40625
        
        
        m=s*2+x+1
        n=t*2+y+1
        
        o=a1*2+a2+1
        oo=b1*2+b2+1
        ooo=c1*2+c2+1
        # print(str(p),str(u),str(q),str(v),str(r),str(w),str(m),str(n),str(o),str(oo),str(ooo))
        # print(o,oo,ooo)
        
        return( str(p)+str(u)+str(q)+str(v)+str(r)+str(w)+str(m)+str(n)+str(o))
    
    
    
    
    
        

if __name__ == "__main__":
    latlon2grid=LatLon2Grid()
    lon = 141.3157
    lat = 43.1714
    print('1次メッシュ： ', latlon2grid.grid1st(lon, lat))  
    print('2次メッシュ： ', latlon2grid.grid2nd(lon, lat))  
    print('基準地域メッシュ： ', latlon2grid.grid3rd(lon, lat))  
    print('4次メッシュ： ', latlon2grid.grid4th(lon, lat))  
    print('11桁 ', latlon2grid.grid_11_digits(lon, lat))
    
    
    
    # 観測点temp.csv に エクセルからlon lat name　をコピー
    # Z-V3-JAPAN-AMP-VS400_M250.csvはJ-SHISからダウンロードした
    # Q-GISのpointssamplingtoosと同じ操作であることを確認済
    
    # stations = pd.read_csv('KiK-net観測点.csv',header = 0,encoding="shift_jis")
    # Mesh_data_250m = pd.read_csv('Z-V3-JAPAN-AMP-VS400_M250.csv',header=None,skiprows=7)
    
    # print(stations)
    # print(Mesh_data_250m)
    
    # cols = ['経度','緯度','メッシュコード','JCODE','AVS','ARV','観測点名']
    # output = pd.DataFrame(index=[], columns=cols) 
    
    # for i,A in enumerate(stations.itertuples()):
    #     # print(A)
    #     name=A[1]
    #     lon = A[4]
    #     lat = A[3]
    #     JCODE,AVS,ARV=0,0,0
    #     # print(lon,lat)
    #     if str(lon).replace('.', '').isdecimal() and lon > 0 and lat > 0:
    #         mesh = latlon2grid.grid250m(lon, lat)
    #         index=Mesh_data_250m[Mesh_data_250m[0]==int(mesh)].index
    #         # print(index)
    #         if index > 1:
    #             grep_data = Mesh_data_250m.loc[index]
    #             # print(grep_data)
    #             # print(grep_data.columns)
    #             # print(grep_data.index)
    #             JCODE=int(grep_data[1])
    #             AVS=float(grep_data[2])
    #             ARV=float(grep_data[3])
    #     else:
    #         lon,lat=0,0
    #         mesh = 0
            
    #     print(i,name,lon,lat,mesh,"JCODE",JCODE,"AVS",AVS,"ARV",ARV)
        
    #     record = pd.Series([lon,lat,mesh,JCODE,AVS,ARV,name] ,  index=output.columns)
    #     output = output.append(record,ignore_index=True)
        
    
    # print(output)
    # output.to_csv("KiK-net250mメッシュの表層地盤情報割り当て後.csv",encoding="shift_jis")
    
    
    
    
    
    
    