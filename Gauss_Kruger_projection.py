import math
import numpy as np
# import sympy


def Positive_calculation(Ra: '参考椭球体长轴', Rb: '参考椭球体短轴', B: '待求点纬度，单位应为度', L: '待求点经度，单位应为度', L0: '中央子午线' = 0, Hp: '投影高，默认为0' = 0, zone: '投影带宽度' = 3, delt_y: '平面坐标平移量' = 500000):
    # 高斯投影正算方法

    # 扁率
    f = (Ra-Rb)/Ra

    # 椭球参数调整，抵偿坐标实现
    if Hp != 0:
        Ra = Hp + Ra
        Rb = Ra*(1-f)

    # 待求点经纬度
    B = math.radians(B)
    L = math.radians(L)

    # 第一偏心率
    e = math.sqrt(2.0*f - f**2)

    # 中间参数
    A1 = 1 + 3/4*e**2 + 45/64*e**4 + 175/256*e**6 + 11025 / \
        16384*e**8 + 43659/65536*e**10 + 693693/1048576 * e**12
    B1 = 3/8*e**2 + 15/32*e**4 + 525/1024*e**6 + 2205 / \
        4096*e**8 + 72765/131072*e**10 + 297297/524288*e**12
    C1 = 15/256*e**4 + 105/1024*e**6 + 2205/16384 * \
        e**8 + 10395/65536*e**10 + 1486485/8388608*e**12
    D1 = 35/3072*e**6 + 105/4096*e**8 + 10395/262144*e**10 + 55055/1048576*e**12
    E1 = 315/131072*e**8 + 3465/524288*e**10 + 99099/8388608*e**12
    F1 = 693/1310720*e**10 + 9009/5242880*e**12
    G1 = 1001/8388608*e**12

    # 子午线弧长
    X = Ra*(1-e**2)*(A1*B-B1*math.sin(2*B)+C1*math.sin(4*B)-D1 *math.sin(6*B)+E1*math.sin(8*B)-F1*math.sin(10*B)+G1*math.sin(12*B))

    # 中间参数
    W = math.sqrt(1-(e*math.sin(B))**2)

    # 卯酉圈半径
    N = Ra / W

    # 中间参数
    t = math.tan(B)

    # 第二偏心率
    e_1 = math.sqrt(Ra**2 - Rb**2)/Rb

    # 中间参数
    eta =  e_1 * math.cos(B)

    # 经度增量L-L0
    dd = 0
    if L0 == 0:
        if zone == 3:
            dd = round(math.degrees(L)/3, 0)
            L0 = dd * 3
        elif zone == 6:
            dd = int(math.degrees(L)/6)+1
            L0 = dd*6-3

    l_D = L - math.radians(L0)

    # 投影坐标x
    x_Gauss = X + N*math.sin(B)*math.cos(B)*l_D**2/2+N*math.sin(B)*math.cos(B)**3*(
        5-t**2+9*eta**2+4*eta**4)*l_D**4/24+N*math.sin(B)*math.cos(B)**5*(61-58*t**2+t**4)*l_D**6/720

    # 投影坐标y
    y_Gauss = N*math.cos(B)*l_D + N*(1-t**2+eta**2)*l_D**3*math.cos(B)**3 / 6 + N*math.cos(B)**5*(5-18*t**2+t**4+14*eta **2-58*t**2*eta**2)*l_D**5/120


    return {'x': x_Gauss, 'y': y_Gauss+delt_y}


def Back_calculation(Ra: '参考椭球体长轴', Rb: '参考椭球体短轴', x_Gauss: '待求点x自然坐标', y_Gauss: '待求点y自然坐标', h:'平面坐标的投影高',L0: '中央子午线',y_delt:'假东'):
    # 高斯投影反算方法
    
        # 扁率
    f = (Ra-Rb)/Ra
    #    print("扁率：" + str(f))
    if h != 0:
        Ra = Ra+h
        Rb = Ra*(1-f)

    y_Gauss = y_Gauss - y_delt
    
    print("待求x坐标：" + str(x_Gauss))
    print("待求y坐标：" + str(y_Gauss))

    # # 中央子午线
    # L0 = 111

    # 扁率
    f = (Ra-Rb)/Ra
    print("扁率：" + str(f))

    # 第一偏心率
    e = math.sqrt(2.0*f - f**2)
    print("第一偏心率：" + str(e))

    # 第二偏心率
    e_1 = math.sqrt(Ra**2 - Rb**2)/Rb
    print("第二偏心率：" + str(e_1))

    # 中间参数
    A1 = 1 + 3/4*e**2 + 45/64*e**4 + 175/256*e**6 + 11025 / \
        16384*e**8 + 43659/65536*e**10 + 693693/1048576.0 * e**12
    B1 = 3/8.0*e**2 + 15/32*e**4 + 525/1024*e**6 + 2205 / \
        4096*e**8 + 72765/131072*e**10 + 297297/524288.0*e**12
    C1 = 15/256*e**4 + 105.0/1024*e**6 + 2205/16384 * \
        e**8 + 10395/65536*e**10 + 1486485/8388608*e**12
    D1 = 35/3072*e**6 + 105/4096*e**8 + 10395/262144*e**10 + 55055/1048576.0*e**12
    E1 = 315/131072*e**8 + 3465/524288*e**10 + 99099/8388608.0*e**12
    F1 = 693/1310720.0*e**10 + 9009/5242880*e**12
    G1 = 1001.0/8388608*e**12

    # 子午线弧长初始化
    X = x_Gauss
    print("子午线弧长：" + str(X))

    # 中间参数
    a0 = (Ra*(1-e**2)*A1)

    # 待求点底点纬度初始化
    Bf0 = X/a0

    # Bf增量迭代算法

    def FB(B):
        FB = Ra*(1-e**2)*(A1*B-B1*math.sin(2*B)+C1*math.sin(4*B)-D1 *
                          math.sin(6*B)+E1*math.sin(8*B)-F1*math.sin(10*B)+G1*math.sin(12*B))
        return FB

    # #FB函数的导数
    def FB_1(B):
        FB_1 = Ra*(1-e**2)*(A1-2*B1*math.cos(2*B)+4*C1*math.cos(4*B)-6*D1 *
                            math.cos(6*B)+8*E1*math.cos(8*B)-10*F1*math.cos(10*B)+12*G1*math.cos(12*B))
        return FB_1

    def calculate_Bf(Bi):
        # 牛顿迭代法，求底点纬度收敛

        BB = (X-FB(Bi))/FB_1(Bi)
        delt = abs(BB)
        # if delt < 8e-17:
        #     return (BB + Bi)
        # else:
        #     return calculate_Bf(BB + Bi)
        while (delt > 8e-17):
            Bi = Bi + BB
            BB = (X-FB(Bi))/FB_1(Bi)
            delt = abs(BB)
        return Bi


    # d迭代计算Bf的值
    Bf = calculate_Bf(Bf0)

    # 卯酉圈半径初始化
    Nf = Ra/math.sqrt(1-e**2*(math.sin(Bf)**2))

    # 子午圈半径初始化
    Mf = Ra*(1-e**2)/(1-e**2*math.sin(Bf)**2)**(-1.5)

    # 中间参数
    tf = math.tan(Bf)
    eta_f = e_1*math.cos(Bf)

    # 待求点纬度
    B = math.degrees(Bf)-math.degrees(tf*y_Gauss**2/2/Mf/Nf)+math.degrees(tf*(5+3*tf**2+eta_f**2-9*eta_f **
                                                                              2*tf**2)*y_Gauss**4/24/Mf/Nf**3)-math.degrees(tf*(61+90*tf**2+45*tf**4)*y_Gauss**6/720/Mf/Nf**5)

    # 待求点经度
    L = math.degrees(y_Gauss/(math.cos(Bf)*Nf)*(1-(1+2*tf**2+eta_f**2)*y_Gauss/Nf**2/6) +
                     (5+28*tf**2+24*tf**4+6*eta_f**2+8*eta_f**2*tf**2)*(y_Gauss/Nf)**4/120)
    L = L+L0


    return {'B': B, 'L': L}
    
    def Positive_calculation_batch(Ra: '参考椭球体长轴', Rb: '参考椭球体短轴', B: '待求点纬度，单位应为度', L: '待求点经度，单位应为度', L0: '中央子午线' = 0, Hp: '投影高，默认为0' = 0, zone: '投影带宽度' = 3, delt_y: '平面坐标平移量' = 500000):
    # 高斯投影正算方法
    
    B.reshape([B.shape[0],1])
    L.reshape([L.shape[0],1])

    # 扁率   标量参数
    f = (Ra-Rb)/Ra

    # 椭球参数调整，抵偿坐标实现
    if Hp != 0:
        Ra = Hp + Ra      # 标量参数
        Rb = Ra*(1-f)     # 标量参数

    # 待求点经纬度   向量参数
    B = np.radians(B)
    L = np.radians(L)

    # 第一偏心率   标量参数
    e = math.sqrt(2.0*f - f**2)


    # 中间参数   标量参数
    A1 = 1 + 3/4*e**2 + 45/64*e**4 + 175/256*e**6 + 11025 / \
        16384*e**8 + 43659/65536*e**10 + 693693/1048576 * e**12
    B1 = 3/8*e**2 + 15/32*e**4 + 525/1024*e**6 + 2205 / \
        4096*e**8 + 72765/131072*e**10 + 297297/524288*e**12
    C1 = 15/256*e**4 + 105/1024*e**6 + 2205/16384 * \
        e**8 + 10395/65536*e**10 + 1486485/8388608*e**12
    D1 = 35/3072*e**6 + 105/4096*e**8 + 10395/262144*e**10 + 55055/1048576*e**12
    E1 = 315/131072*e**8 + 3465/524288*e**10 + 99099/8388608*e**12
    F1 = 693/1310720*e**10 + 9009/5242880*e**12
    G1 = 1001/8388608*e**12

    # 子午线弧长  向量参数
    X = Ra*(1-e**2)*(A1*B-B1*np.sin(2*B)+C1*np.sin(4*B)-D1 *np.sin(6*B)+E1*np.sin(8*B)-F1*np.sin(10*B)+G1*np.sin(12*B))
    
    # 中间参数  向量参数
    W = np.sqrt(1 - np.power(e * np.sin(B),2))
    
    # 卯酉圈半径  向量参数
    N = Ra / W
    
    # 中间参数  向量参数
    t = np.tan(B)

    # 第二偏心率   标量参数
    e_1 = math.sqrt(Ra**2 - Rb**2)/Rb
    
    # 中间参数  向量参数
    eta =  e_1 * np.cos(B)
    
    print("e: ",e)
    print("eta:",eta)
    print("e_1:",e_1)
    print("t:",t)
    print("N:",N)
    print("W:",W)
    print("X:",X)
    
    # 经度增量L-L0   向量参数
    dd = np.array([])
    if L0 == 0:
        if zone == 3:
            dd = np.round(np.degrees(L)/3, 0)
            L0 = dd * 3
        elif zone == 6:
            dd = np.trunc(np.degrees(L)/6)+1
            L0 = dd*6-3    
    l_D = L - math.radians(L0)

    # 投影坐标x
    x_Gauss = X + N * np.sin(B) * np.cos(B) * np.power(l_D,2)/2 + N * np.sin(B) * np.power(np.cos(B),3) * (5 - np.power(t,2) + 9 * np.power(eta,2) + 4 * np.power(eta,4)) * np.power(l_D,4) /24 + N * np.sin(B) * np.power(np.cos(B),5) * (61 - 58 * np.power(t,2) + np.power(t,4)) * np.power(l_D,6) / 720

    # 投影坐标y
    y_Gauss = N * np.cos(B) * l_D + N * (1 - np.power(t,2) + np.power(eta,2)) * np.power(l_D,3) * np.power(np.cos(B),3) / 6 + N * np.power(np.cos(B),5) * (5 - 18 * np.power(t,2) + np.power(t,4) + 14 * np.power(eta ,2) - 58 * np.power(t,2) * np.power(eta,2)) * np.power(l_D,5)/120
    
    x_Gauss = np.round(x_Gauss,4)
    y_Gauss = np.round(y_Gauss,4)

    return {'x': x_Gauss, 'y': y_Gauss+delt_y}
    
def Back_calculation_batch(Ra: '参考椭球体长轴', Rb: '参考椭球体短轴', x_Gauss: '待求点x自然坐标', y_Gauss: '待求点y自然坐标', h:'平面坐标的投影高',L0: '中央子午线',y_delt:'假东'):
    # 高斯投影反算方法

    # 扁率
    f = (Ra-Rb)/Ra
    
    # 修改椭球参数
    if h != 0:
        Ra = Ra+h
        Rb = Ra*(1-f)
    
    # 格式化矩阵维度
    x_Gauss.reshape([x_Gauss.shape[0],1])
    y_Gauss.reshape([y_Gauss.shape[0],1])
    
    # 去除假东
    y_Gauss = y_Gauss - y_delt

    # 第一偏心率
    e = math.sqrt(2.0*f - f**2)

    # 第二偏心率
    e_1 = math.sqrt(Ra**2 - Rb**2)/Rb

    # 中间参数
    A1 = 1 + 3/4*e**2 + 45/64*e**4 + 175/256*e**6 + 11025 / \
        16384*e**8 + 43659/65536*e**10 + 693693/1048576.0 * e**12
    B1 = 3/8.0*e**2 + 15/32*e**4 + 525/1024*e**6 + 2205 / \
        4096*e**8 + 72765/131072*e**10 + 297297/524288.0*e**12
    C1 = 15/256*e**4 + 105.0/1024*e**6 + 2205/16384 * \
        e**8 + 10395/65536*e**10 + 1486485/8388608*e**12
    D1 = 35/3072*e**6 + 105/4096*e**8 + 10395/262144*e**10 + 55055/1048576.0*e**12
    E1 = 315/131072*e**8 + 3465/524288*e**10 + 99099/8388608.0*e**12
    F1 = 693/1310720.0*e**10 + 9009/5242880*e**12
    G1 = 1001.0/8388608*e**12

    # 子午线弧长初始化
    X = x_Gauss
    print(X[0])

    # 中间参数
    a0 = (Ra * ( 1 - e**2 ) * A1)

    # 待求点底点纬度初始化
    Bf0 = X / a0
    print(Bf0[0])

    # Bf增量迭代算法

    def FB(B):
        FB = Ra * (1 - e**2) * (A1 * B - B1 * np.sin(2*B)+C1*np.sin(4*B)-D1 * np.sin(6*B) + E1 * np.sin(8*B) - F1 * np.sin(10*B) + G1 * np.sin(12*B))
        return FB

    # #FB函数的导数
    def FB_1(B):
        FB_1 = Ra * (1 - e**2) * (A1 - 2 * B1 * np.cos(2*B) + 4 * C1 * np.cos(4*B) - 6 * D1 * np.cos(6*B) + 8 * E1 * np.cos(8*B) - 10 * F1 * np.cos(10*B) + 12 * G1 * np.cos(12*B))
        return FB_1

    def calculate_Bf(X,Bi):
        # 迭代求底点纬度收敛
        BB = (X-FB(Bi))/FB_1(Bi)
        delt = abs(BB)
        while (delt > 8e-17):
            Bi = Bi + BB
            BB = (X-FB(Bi))/FB_1(Bi)
            delt = abs(BB)
        return Bi
    
    Bf = []
    for ii in range(0,Bf0.shape[0]):
        # d迭代计算Bf的值
        Bf.append(calculate_Bf(X[ii],Bf0[ii]))
    Bf = np.array(Bf)
    Bf = Bf.reshape([Bf.shape[0],1])
        
        
    # 卯酉圈半径初始化
    Nf = Ra / np.sqrt(1 - e**2 * np.power(np.sin(Bf),2))

    # 子午圈半径初始化
    Mf = Ra*(1 - e**2) / np.power((1 - e**2 * np.power(np.sin(Bf),2)),(-1.5))

    # 中间参数
    tf = np.tan(Bf)
    eta_f = e_1 * np.cos(Bf)

    # 待求点纬度
    B = np.degrees(Bf) - np.degrees(tf * np.power(y_Gauss,2) /2/Mf/Nf) + np.degrees(tf * (5 + 3 * np.power(tf,2) + np.power(eta_f,2) - 9 * np.power(eta_f,2) * np.power(tf,2)) * np.power(y_Gauss,4) /24/Mf/np.power(Nf,3)) - np.degrees(tf * (61+90* np.power(tf,2) +45* np.power(tf,4) ) * np.power(y_Gauss,6) /720/Mf/ np.power(Nf,5))

    # 待求点经度
    L = np.degrees(y_Gauss/(np.cos(Bf)*Nf)*(1-(1+2*tf**2+eta_f**2)*y_Gauss/Nf**2/6) + (5+28*tf**2+24*tf**4+6*eta_f**2+8*eta_f**2*tf**2)*(y_Gauss/Nf)**4/120)
    L = L+L0
    
    B = np.round(B,4)
    L = np.round(L,4)

    return {'B': B, 'L': L}


