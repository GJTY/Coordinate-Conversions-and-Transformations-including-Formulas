from math import radians, degrees, sqrt
from math import acos, asin, atan2
from math import cos, sin, tan


def BLH2XYZ(Ra: dict(type=float, help='参考椭球体的长轴'), Rb: dict(type=float, help='参考椭球体的短轴'), B: dict(type=float, help='待求点的纬度，度为单位'), L: dict(type=float, help='待求点的经度，度为单位'), H: dict(type=float, help='待求点的大地高高程')) -> dict:
    '''
    计算经纬度转换到到空间直角坐标的方法
    '''

    # 扁率
    f = (Ra-Rb)/Ra
    # print("扁率：" + str(f))

    # 第一偏心率
    e = sqrt(2*f - f**2)
    # print("第一偏心率：" + str(e))

    # 中间参数
    W = sqrt(1-(e*sin(B))**2)

    # 卯酉圈半径
    N = Ra / W
    # print("卯酉圈半径："+str(N))

    X = (N+H)*cos(B)*cos(L)
    Y = (N+H)*cos(B)*sin(L)
    Z = (N*(1-e**2)+𝐻)*sin(B)

    # print(str(X)+','+str(Y)+','+str(Z))
    # input()
    return {'X': X, 'Y': Y, 'Z': Z}


def XYZ2BLH(Ra: dict(type=float, help='参考椭球体的长轴'), Rb: dict(type=float, help='参考椭球体的短轴'), X: dict(type=float, help='待求点空间直角坐标系X'), Y: dict(type=float, help='待求点空间直角坐标系Y'), Z: dict(type=float, help='待求点空间直角坐标系Z')) -> dict:
    '''
    计算空间直角坐标转换到经纬度的方法
    '''

    # 扁率
    f = (Ra-Rb)/Ra
    print("扁率：" + str(f))

    # 第一偏心率
    e = sqrt(2*f - f**2)
    # e = round(e,7)
    print("第一偏心率：" + str(e))

    # 经度
    L = atan2(Y, X)

    # 接下来迭代计算B

    # 中间参数
    AA = sqrt(X**2+Y**2)
    # 初始值
    B0 = atan2(Z, AA)

    def caculat_B(Bi):
        # 递归方法，迭代计算B

        # 中间参数
        W = sqrt(1-(e*sin(Bi))**2)
        # 卯酉圈半径
        N = Ra / W
        B = atan2(Z+N*e**2*sin(Bi), AA)
        delt = abs(B - Bi)
        if delt < 1e-17:
            return B
        else:
            return caculat_B(B)

    B = caculat_B(B0)

    # 重新计算中间参数W
    W = sqrt(1-(e*sin(B))**2)
    # 重新计算卯酉圈半径N
    N = Ra / W

    H = Z/sin(B)-N*(1-e**2)
    # H = AA/cos(B) - N

    return {'B':B,'L':L,'H':H}
