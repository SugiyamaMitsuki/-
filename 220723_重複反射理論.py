import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmath

class Multiple_Reflection_theory():
    def __init__(self):

        self.temp = None

        None
        
    def multiple_reflection_theory_rogic( self, N_layers, rho, Vs, H, h, data_length, dt ):
        """
        重複反射理論
        地表の複素振幅の値を1として各層の複素振幅を周波数領域で求める
        """
        
        #計算に使用する周波数軸
        freq = np.linspace(0, 1.0/dt, data_length)  #numpy.linspace(a, b, n)メソッドでa以上b以下で個数nの等差数列（配列）を作成ができます 
        #freq = np.arange(0, 1.0/dt, df)            #numpy.arrange(a, b, d)メソッドでa以上b未満で間隔dの等差数列を作ることができます。   
        
        #複素振幅の配列を準備  
        E =  np.empty([N_layers,data_length],dtype=np.complex128)   
        F =  np.empty([N_layers,data_length],dtype=np.complex128)
        
        # 各層の値
        Gstar =  np.empty([N_layers],dtype=np.complex128)
        R     =  np.empty([N_layers],dtype=np.complex128)
        
        # 地表の複素振幅E,Fを1とする
        j = 0
        E[j, :] = 1
        F[j, :] = 1
        
        # 複素振幅の計算に使用する各層固有の値
        for j in range(N_layers):
            #複素剛性率　(実部が剛性、虚部が非弾性減衰の大きさ)            
            Gj = rho[j] * (Vs[j] **2)            
            Gstar[j]=(1 + 2*h[j]*1j)*Gj

        # 複素振幅の計算に使用する隣接する層と層での関係値
        for j in range(N_layers-1):
            #インピーダンス比 第1層と第2層の関係　から　第n-1層と第n層の関係
            R[j] = cmath.sqrt((rho[j]*Gstar[j]) / (rho[j+1]*Gstar[j+1]))

        # 地表から下に向かって順次複素振幅を求める
        for j in range(N_layers-1):
            print(f"N_layers_step = {j}")
            
            for f_num,f in enumerate(freq):
                # print(f"{f_num} / {data_length}  f = {f}")

                # 伝播定数(減衰の効果を含めた複素数の波数)　周波数の関数
                pj= (2 * np.pi * f) * cmath.sqrt(rho[j]/Gstar[j])
                
                #複素振幅の漸化式 E上昇波　F下降波
                E[j+1,f_num] = 1/2 * ( 1 + R[j] ) * cmath.exp(  1j * pj * H[j] ) * E[j,f_num] \
                            +  1/2 * ( 1 - R[j] ) * cmath.exp( -1j * pj * H[j] ) * F[j,f_num]

                F[j+1,f_num] = 1/2 * ( 1 - R[j] ) * cmath.exp(  1j * pj * H[j] ) * E[j,f_num] \
                            +  1/2 * ( 1 + R[j] ) * cmath.exp( -1j * pj * H[j] ) * F[j,f_num]
                 
        return E,F,freq
    
    
        

if __name__ == '__main__':
    
    ## 重複反射理論用のクラス
    rogic = Multiple_Reflection_theory()
    
    # read 地盤データ 
    # 層厚[m]  減衰定数[-]  密度[g/cm3]  S波速度[m/s]
    input_path = "地盤データサンプル.csv"
    soil = pd.read_csv(input_path,header=0)

    print(soil)

    # soil = soil[0:2]

    # 地盤条件の整理
    N_layers = len(soil)
    rho = soil["密度[g/cm3]"].values
    Vs = soil["S波速度[m/s]"].values
    H = soil["層厚[m]"].values
    h = soil["減衰定数[-]"].values
    dt=0.02
    data_length=1024

    ## 重複反射理論で複素振幅を計算
    E,F,freq = rogic.multiple_reflection_theory_rogic(
                                                    N_layers = N_layers,
                                                    rho = rho,
                                                    Vs = Vs,
                                                    H = H,
                                                    h = h,
                                                    data_length = data_length,
                                                    dt = dt
                                                    )

    # 増幅スペクトル
    r=0
    s=N_layers-1

    #露頭基盤 2E波
    Z_out =  (2*E[r]) / (2*E[s])
    #内部基盤
    Z_in  =  (2*E[r]) / (E[s]+F[s])

    phaseSpectrum_out = [np.arctan2(c.imag, c.real) for c in Z_out]    # 位相スペクトル
    phaseSpectrum_in = [np.arctan2(c.imag, c.real) for c in Z_in]    # 位相スペクトル



    # 二層地盤でのプログラム確かめ用
    freq0 = Vs[0]/(4*H[0])
    x=freq/freq0
    x=freq
    
    plt.rcParams['font.family'] = 'Hiragino Sans' #日本語使用可

    fig = plt.figure() # Figureを作成
    ax1 = fig.add_subplot(3,1,1) # Axesを作成
    ax1.plot(x,abs(abs(Z_in)), label="内部基盤面（E+F）") 
    ax1.plot(x,abs(abs(Z_out)),label="解放基盤面（2E）")
    ax1.set_xlabel('f')
    ax1.set_ylabel('Zr/s(f)')
    # ax1.set_yscale('log')
    ax1.set_xlim(0,6)
    ax1.minorticks_on()
    ax1.legend()
    
    ax2 = fig.add_subplot(3,1,2) # Axesを作成
    ax2.plot(x,phaseSpectrum_out,label="解放基盤面（2E）位相スペクトル")
    ax2.set_xlabel('f')
    ax2.set_ylabel('位相差')
    ax2.minorticks_on()
    ax2.set_xlim(0,6)
    ax2.legend()

    ax3 = fig.add_subplot(3,1,3) # Axesを作成
    ax3.plot(x,phaseSpectrum_in,label="内部基盤面（E+F）位相スペクトル")
    ax3.set_xlabel('f')
    ax3.set_ylabel('位相差')
    ax3.minorticks_on()
    ax3.set_xlim(0,6)
    ax3.legend()

    plt.savefig("増幅スペクトル.png")
    plt.show()