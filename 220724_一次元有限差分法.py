import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import ArtistAnimation
import matplotlib.animation as animation

import random
import sys 
from scipy import interpolate
import itertools



#パラメータの設定
dz = 0.1
dt = 0.0001

Number_of_timestep = int(5/dt)


#各層の物性値を設定
# read 地盤データ 
# 層厚[m]  減衰定数[-]  密度[g/cm3]  S波速度[m/s]
input_path = "地盤データサンプル.csv"
soil = pd.read_csv(input_path,header=0)
print(soil)
# 地盤条件の整理
rho = soil["密度[g/cm3]"].values
rho=np.round(rho,2)
Velocity = soil["S波速度[m/s]"].values
height = soil["層厚[m]"].values



_rho = []
_Velocity = []
#各層の物性値をdzで分割する
for i,h in enumerate(height):
    length = int(h/dz)
    _rho.append([rho[i]] * length)
    _Velocity.append([Velocity[i]]*length)

rho = np.array(list(itertools.chain.from_iterable(_rho)))
Velocity = np.array(list(itertools.chain.from_iterable(_Velocity)))

G = rho * Velocity **2 #(kN/m2)
depth = sum(height)
Number_of_layers = len(G) #int(depth/dz)



#入力地震動
t_data =np.linspace(0.0, dt*Number_of_timestep,Number_of_timestep)
# input_wave  = [1 * np.cos(2*np.pi*k*dt) for k in range(Number_of_timestep)]  
input_wave  = [0 for k in range(Number_of_timestep)]
input_wave[100] = 1

# 配列の用意
#速度 v[x,t]
v =  np.zeros([Number_of_layers,Number_of_timestep], dtype = 'float64') #配列の箱としては一つ余分に用意
#応力 T[x,t]
T =  np.zeros([Number_of_layers,Number_of_timestep], dtype = 'float64')


#初期条件
#境界条件 

fig = plt.figure()
ims = []


gamma2=9/8
gamma4=-1/24

alfa = np.exp(-1/(2*200)*dt) #減衰　Q＝200
# alfa = 1

#時刻ループ
for k in range(Number_of_timestep-1):

    print(f"step={k}/{Number_of_timestep}")
    
    # 応力 K+1

    #地表境界条件 Z=0
    z = 0
    T[z,k+1] =  0.0      #地表の応力は0  -T[z-1,k]

    #深さループ
    for z in range(1,Number_of_layers-1): # 1 ~ (Number_of_layers-1)-1
        # print(f"k={k} z={z} 深さ{z*dz}")
        # 4次差分　 
        T[z,k+1] = T[z,k]\
                    +G[z] * ( gamma2 * ( v[z,k]   - v[z-1,k] ) / dz \
                          +   gamma4 * ( v[z+1,k] - v[z-2,k] ) / dz ) * dt * alfa
        
    #最下層
    z=Number_of_layers - 1
    T[z,k+1] = T[z,k] + G[z] * ( 1 / ( 2 * dz ) ) * ( v[z,k] - v[z-1,k] ) * dt


    # 速度 K+1+1/2

    #地表境界条件
    z = 0
    v[z,k+1] = v[z,k] + (1/rho[z])*(1/(2*dz))*(T[z+1,k+1]-T[z,k+1])*dt 

    #深さループ
    for z in range(1,Number_of_layers-2): # 1 ~ (Number_of_layers-2)-1
        # print(f"k={k} z={z} 深さ{z*dz}")
        # 4次差分
        v[z,k+1] = v[z,k]\
                    +(1/rho[z]) * ( gamma2 * ( T[z+1,k+1] - T[z  ,k+1]) / dz \
                                +   gamma4 * ( T[z+2,k+1] - T[z-1,k+1]) / dz ) * dt
        
    #最下層の一つ上 2次差分
    z = Number_of_layers - 2 
    v[z,k+1] = v[z,k] + ( 1 / rho[z] ) * ( 1 / ( 2 * dz ) ) * ( T[z+1,k+1] - T[z,k+1] ) * dt

    #最下層
    z = Number_of_layers -1
    v[z,k+1] = input_wave[k]    #←入力地震動


    if k%10 == 0:
        x=[i*dz for i in range(Number_of_layers)]
        im = plt.plot(x,v[:,k+1], linewidth=1, color="blue")
        # im = plt.scatter( [i*dz for i in range(Number_of_layers)] , v[:,k+1], color="blue" )
        ims.append(im) 






# fig = plt.figure()
# ims = []
# im = plt.plot(x[1:N+1], u[1:N+1], linewidth=2, color="blue")
# ims.append(im) 
plt.title("")
plt.xlabel("depth")
plt.ylabel("v")
plt.ylim(-5,5)
ani = animation.ArtistAnimation(fig, ims, interval=100, repeat=False)
plt.show()
# ani.save("anim.gif", writer="pillow", fps=10)








plt.plot(input_wave)
plt.show()

for i in [-1,0]:
    fig_xy = plt.figure()
    plt.title(f"v[{i},:]")
    y=v[i,:]
    plt.plot(t_data,y)
    plt.xlabel('time')
    plt.show()
