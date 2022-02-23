"""
1D　Wave Equation
FTCS法
https://qiita.com/sci_Haru/items/8535f435ffa1febcd445
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import ArtistAnimation # アニメーション作成のためのメソッドをインポート
import time
import random
from scipy import interpolate
import math

def Draw_graph(title,x_label,y_label,x,y,output_path="temp"):
    plt.figure(figsize=(10.0, 4.0))
    title = "{}".format(title)
    x_label = x_label
    y_label = y_label
    x = np.array(x)
    y = np.array(y)
    plt.plot(x, y, label=title)
    plt.ylim(-max(abs(y)),max(abs(y)))
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    #plt.grid()
    leg = plt.legend(loc=1, fontsize=15)
    leg.get_frame().set_alpha(1)
    plt.savefig("./output/{0}.png".format(output_path))
    plt.show()
    
def save_gnu_data(x,y,columns,filename):
    df = pd.DataFrame()
    df[0]=x
    df[1]=y
    df.columns=columns
    filename = filename
    df.to_csv("{}".format(filename), header=False, index=False, sep="\t")
    with open("{}".format(filename)) as f:
        l = f.readlines()
    l.insert(0, '#.gnu_txt\n')
    with open("{}".format(filename), mode='w') as f:
        f.writelines(l)

t_max = 30
dx = 5

#各層の物性値を設定
N_layers = 2
#各層の物性値を設定
rho_input =  [1.5,5] #(g/cm^3)
V_input   =  [100,400] #(m/s)
H   =  [20,60] #(m)

G_depth = 5*400*400 
rho_depth=10000 #分母大きくしてに近くする

depth = sum(H)
dt = dx/max(V_input)/300
print("dt=",dt)
Number_of_timestep = int(t_max//dt)
 
print("Number_of_timestep=",Number_of_timestep)
num_lattice = depth/dx
print("num_lattice",num_lattice)
if num_lattice.is_integer():
    num_lattice = int(depth/dx)
else:
    print("stop")
    l
    time.sleep(100)
    
rho =  np.zeros([num_lattice])
V =  np.zeros([num_lattice])

#各層の物性値をdxで分割する
for i in range(len(H)):
    if i == 0:
        end = int(H[i]/depth*num_lattice)
        rho[:end] = rho_input[i]
        V[:end]=V_input[i]
    else:        
        start = int(sum(H[:i])/depth*num_lattice)
        end = int(sum(H[:i+1])/depth*num_lattice)
        rho[start:end] = rho_input[i]
        V[start:end] = V_input[i]
G = rho * (V **2)  

#time = np.arange(0, Number_of_timestep*dt, dt)  # 時間軸
t_data =np.linspace(0, dt*Number_of_timestep,Number_of_timestep) #numpy.linspace(a, b, n)メソッドでa以上b以下で個数nの等差数列（配列）を作成ができます     
print("tmax=",max(t_data))
freq = np.linspace(0, 1.0/dt, Number_of_timestep) 

#入力地震動
filename = "HYG025.vx.vel"
input_path = "{0}".format(filename)
df = pd.read_table(input_path,skiprows=[0])
x_observed = df.iloc[:,0]
y_observed = df.iloc[:,1]
dt_observed=df.iloc[1,0]-df.iloc[0,0]
x_ratio = dt_observed/dt
#時間刻みが違うから線形補完する
#ip2 = ["線形補間", interpolate.interp1d]
#ip7 = ["3次スプライン補間", lambda x, y: interpolate.interp1d(x, y, kind="cubic")]
x_latent = np.linspace(min(x_observed), max(x_observed),int(len(x_observed)*x_ratio)) 
fitted_curve = interpolate.interp1d(x_observed, y_observed)

input_wave = fitted_curve(x_latent)
print("len_input_wave=",len(input_wave))
plt.scatter(x_observed, y_observed, label="observed")
plt.plot(x_latent, input_wave, c="red", label="fitted")
plt.grid()
plt.legend()
plt.show()


#速度 v[x,t]
v =  np.zeros([num_lattice+1,Number_of_timestep]) #一層多く用意
#応力 T[x,t]
T =  np.zeros([num_lattice+1,Number_of_timestep])


#初期条件#初期条件
v[:,0]=0
T[:,0]=0

#境界条件 地表の応力０
T[0,:] = 0


gamma2=9/8
gamma4=-1/24
alfa = np.exp(-1/(2*400)*dt)
print("alfa",alfa)
alfa = 1
for k in range(Number_of_timestep-1):
    k_v=k-1/2
    #print("step=",k)
    
    #応力の最上層   
    T[0,k+1] = 0
    
    #速度の格子を下にずらして考える
    i=1
    #二次差分        
    #応力の更新式
    T[i,k+1] = T[i,k] + G[i-1]*(1/(2*dx))*(v[i,k]-v[i-1,k])*dt*alfa
    #速度の更新式
    i_v=i-1
    v[i_v,k+1] = v[i_v,k] + (1/((rho[i]+rho[i-1])/2))*(1/(2*dx))*(T[i,k]-T[i-1,k])*dt*alfa
    
    for i in range(2,num_lattice): #四次差分で計算できる範囲  2 <= i < num_lattice-1
        
        i_v=i-1/2-1
        #四次差分
        #応力の更新式
        T[i,k+1] = T[i,k]\
                    +G[i-1]*(gamma2*(v[int(i_v    +1/2),int(k_v+1/2)]-v[int(i_v-1+1/2),int(k_v+1/2)])/dx\
                        +gamma4*(v[int(i_v+2+ 1/2),int(k_v+1/2)]-v[int(i_v-2 +1/2),int(k_v+1/2)])/dx)*dt*alfa
        
        #速度の更新式
        v[int(i_v-1/2),int(k_v+1+1/2)] = v[int(i_v-1/2),int(k_v+1/2)]\
                    +(1/((rho[i]+rho[i-1])/2))*(gamma2*(T[i,k]-T[i-1,k])/dx\
                              +gamma4*(T[i+1,k]-T[i-2,k])/dx)*dt*alfa
        
        
    i=num_lattice #n層
    #二次差分        
    #応力の更新式
    T[i,k+1] = T[i,k] + G_depth*(1/(2*dx))*(v[i,k]-v[i-1,k])*dt*alfa
    #速度の更新式
    i_v=i-1
    v[i_v,k+1] = v[i_v,k] + (1/((rho_depth+rho[i-1])/2))*(1/(2*dx))*(T[i,k]-T[i-1,k])*dt*alfa
    
    #基盤　
    #速度の最下層
    #n+1層
    v[num_lattice,k+1] = input_wave[k+1]        #←入力地震動
        
#np.cos(2*np.pi*k)+np.sin(np.pi*k)#

#上をせいにする
v[:num_lattice-1,:]=v[:num_lattice-1,:]*10
Draw_graph(title="surface"  ,x_label="time [s]",y_label="vel [cm/s]",x=t_data,y=v[0,:],output_path="wave_1D_surface.png")
# Draw_graph(title="v[num_lattice//4,:]"  ,x_label="time [s]",y_label="vel [cm/s]",x=t_data,y=v[num_lattice//4,:],output_path="wave_1D")
# Draw_graph(title="v[num_lattice//2,:]"  ,x_label="time [s]",y_label="vel [cm/s]",x=t_data,y=v[num_lattice//2,:],output_path="wave_1D")
# Draw_graph(title="v[num_lattice-1,:]"  ,x_label="time [s]",y_label="vel [cm/s]",x=t_data,y=v[num_lattice-1,:],output_path="wave_1D")
Draw_graph(title="v[num_lattice,:]input"  ,x_label="time [s]",y_label="vel [cm/s]",x=t_data,y=-v[num_lattice,:],output_path="wave_1D")

save_gnu_data(x=t_data,y=v[0,:],columns=["time", "dis[cm]"] ,filename="./output/FDM1d_surface.dis")
save_gnu_data(x=t_data,y=v[num_lattice,:],columns=["time", "dis[cm]"] ,filename="./output/FDM1d_input_depth.dis")

#T=4*depth/Vs
T=4*depth/V[0]
print(T,1/T)
f_range = int(10//(1.0/dt))
data_length = int(1024*4)
freq = np.linspace(0, 1.0/dt, data_length)
FFT = np.fft.fft(v[0,:data_length])
Draw_graph(title="wave_surface_out"  ,x_label="f [Hz]",y_label="wave1D_fft  [velsec]",x=freq[:f_range],y=np.abs(FFT[:f_range])  ,output_path="wave1D_fft")


hk


# FFT = np.fft.fft(v[0,:])/(Number_of_timestep/2)
# plt.plot(freq[:Number_of_timestep//1000],np.abs(FFT[:Number_of_timestep//1000]))
# plt.show()
# a
x=list(range(num_lattice))
y=list(range(Number_of_timestep))
X, Y = np.meshgrid(x,y)

def functz(u):
    z=u[X,Y]
    return z

Z = functz(v)
fig = plt.figure()
ax = Axes3D(fig)
ax.view_init(elev=210., azim=60.) # アングル設定
ax.plot_wireframe(X,Y,Z, color='black')
ax.set_xlabel('depth')
ax.set_ylabel('time_shift')
ax.set_zlabel('Vel.')
plt.savefig("{0}.png".format("wave_1D_timeshift"))
plt.show()



fig = plt.figure()

anim = [] #アニメーション用に描くパラパラ図のデータを格納するためのリスト

for i in range(Number_of_timestep):
    U=list(v[:,i])
    x=list(range(num_lattice+1))
    if i % int(Number_of_timestep*0.1) ==0: 
        im=plt.plot(x,U, '-', color='red',markersize=10, linewidth = 2, aa=True)
        anim.append(im)

anim = ArtistAnimation(fig, anim) # アニメーション作成    
plt.xlabel('x')
plt.ylabel('v')
#fig.show() 
anim.save("t.gif", writer='imagemagick')   #アニメーションをt.gifという名前で保存し，gifアニメーションファイルを作成する。
