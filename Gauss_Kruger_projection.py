import math
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


def Back_calculation(Ra: '参考椭球体长轴', Rb: '参考椭球体短轴', x_Gauss: '待求点x自然坐标', y_Gauss: '待求点y自然坐标', L0: '中央子午线'):
    # 高斯投影反算方法

    y_Gauss = y_Gauss-500000
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
        # 迭代求底点纬度收敛

        BB = (X-FB(Bi))/FB_1(Bi)
        delt = abs(Bi)
        if delt < 0.00000000000001:
            return (BB+Bi)
        else:
            return calculate_Bf(BB+Bi)


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
