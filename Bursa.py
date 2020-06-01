import numpy as np
from numpy import sin,cos
from numpy import arcsin as asin
from numpy import arccos as acos

class Bursa():
    '''
    布尔莎七参数模型，采用φ，Ω，κ转角系统；
    模型采用全微分最小二乘迭代拟合七参数，七参数改正数拟合精度上限为1.0e-8；
    此模型可适应大角度旋转下的参数计算，因此无应用方向上的限制
    
    Function：
    
        def bursa(T,phi,omiga,kappa,n,coords_S):
            计算七参数转换正算的函数。利用七参数和转换前的坐标计算出转换后的坐标。
        
        def adjustment(T,phi,omiga,kappa,n,coords_S,coords_T):
            平差计算七参数改正数的函数；该函数在计算七参数模型全微分的基础之上，在给定七参数的初始值T0，phi0，omiga0，kappa0，n0处
            执行泰勒展开，之后舍去高次小项，仅保留一阶项；在此基础之上对原线性化方程组进行改化并计算七参数在初始值下的改正数。函数同
            时返回当前参数下的精度评价指标。
            
        def bursa_Paramenter(coords_S,coords_T):
            利用七参数的改正数，在初始化（T= [0.0,0.0,0.0],phi = 0.0, omiga = 0.0, kappa=0.0, n = 1.0）的初始条件下进行迭代计算，
            直至七参数改正数全部小于1.0e-8，迭代结束，并返回最终的七参数优化结果。
        
        def accuracy_Evaluation(coords_S,coords_T,X_T):
            精度评价的函数，输入不同的点作为输入数据可进行內符精度和外符精度的校验。
    
    '''
    
    @staticmethod
    def bursa(T:np.array,phi:float,omiga:float,kappa:float,n:float,coords_S:np.array):
        '''
        布尔莎大角度七参数转换模型正算，采用phi，omiga，kappa转角系统，

        Parameter:
            T  - 平移参数向量，3 x 1 T = [Tx,Ty,Tz]
            phi  - 绕Z轴的旋转角φ
            omiga  - 绕Y轴的旋转角Ω
            kappa  - 绕X轴的旋转角κ
            n  - 尺度缩放，n = m + 1
            coords_S  - 转换前的点矩阵，(m * 3)  coords_S = [Xs  Ys  Zs].T
            
        Return:
            X_T  - 转换后的点位坐标，m * 3
            
        '''
        # 中间参数
        Xs = coords_S[:,0]
        Ys = coords_S[:,1]
        Zs = coords_S[:,2]

        num_p = coords_S.shape[0]


        # 中间参数
        R1 = cos(phi)*cos(omiga)
        R2 = cos(phi)*sin(omiga)*sin(kappa)+sin(phi)*cos(kappa)
        R3 = sin(phi)*sin(kappa)-cos(phi)*sin(omiga)*cos(kappa)
        R4 = -sin(phi)*cos(omiga)
        R5 = cos(phi)*cos(kappa)-sin(phi)*sin(omiga)*sin(kappa)
        R6 = sin(phi)*sin(omiga)*cos(kappa) + cos(phi)*sin(kappa)
        R7 = sin(omiga)
        R8 = -cos(omiga)*sin(kappa)
        R9 = cos(omiga)*cos(kappa)

        # 旋转矩阵R
        R = np.array([
            [R1,R2,R3],
            [R4,R5,R6],
            [R7,R8,R9]
        ])

        X_T = T.reshape([3,1]) + n*R.dot(coords_S.T)

        return X_T.T
    
    @staticmethod
    def adjustment(T:np.array,phi:float,omiga:float,kappa:float,n:float,coords_S:np.array,coords_T:np.array):
        '''
        布尔莎大角度七参数转换模型参数平差，采用phi，omiga，kappa转角系统，
        该函数采用全微分方程最小二乘迭代拟合法求取参数
        
        Parameter:
            T  - 平移参数向量初始值，3 x 1
            phi  - 绕Z轴的旋转角φ初始值
            omiga  - 绕Y轴的旋转角Ω初始值
            kappa  - 绕X轴的旋转角κ初始值
            n  - 尺度缩放参数初始值，n = m + 1
            coords_S  - 转换前的点矩阵，(m * 3)  coords_S = [Xs  Ys  Zs].T
            coords_T  - 转换后的点矩阵，(m * 3)  coords_S = [Xs  Ys  Zs].T
            
        Return:
            p  - 参数向量的改正数
            sg  - 精度评价指标，模型整体均方根误差，但并不除以m或者m-1，而是除以模型自由度3m - 7(以m为公共点个数)
            
        '''
        # 中间参数
        Xs = coords_S[:,0]
        Ys = coords_S[:,1]
        Zs = coords_S[:,2]

        num_p = coords_S.shape[0]


        # 中间参数
        M1 = (cos(phi)*(cos(kappa)*Ys + sin(kappa)*Zs) - Xs*sin(phi)*cos(omiga) + sin(phi)*sin(omiga)*(Zs*cos(kappa) - Ys*sin(kappa)))
        M2 = (cos(phi)*cos(omiga)*(Ys*sin(kappa) - Zs*cos(kappa)) - Xs*cos(phi)*sin(omiga))
        M3 = (cos(phi)*sin(omiga)*(Ys*cos(kappa) + Zs*sin(kappa))+sin(phi)*(Zs*cos(kappa) - Ys*sin(kappa)))
        M4 = (-sin(phi)*(Zs*sin(kappa) + Ys*cos(kappa))-Xs*(cos(phi)*cos(omiga))+cos(phi)*sin(omiga)*(Zs*cos(kappa)-Ys*sin(kappa)))
        M5 = (sin(phi)*cos(omiga)*(Zs*cos(kappa) - Ys*sin(kappa)) + Xs*(sin(phi)*sin(omiga)))
        M6 = (-sin(phi)*sin(omiga)*(Ys*cos(kappa) + Zs*sin(kappa)) + cos(phi)*(Zs * cos(kappa) - Ys*sin(kappa)))
        M7 = np.array([0]*num_p).reshape(num_p)
        M8 = (sin(omiga)*(Ys*sin(kappa) - Zs*cos(kappa)) + Xs*cos(omiga))
        M9 = (-cos(omiga)*(Ys*cos(kappa) + Zs*sin(kappa)))

        # M
        M = np.array([
            [M1,M2,M3],
            [M4,M5,M6],
            [M7,M8,M9]
            ])

        M = M.transpose(0,2,1).reshape(3*num_p,3)


        # 中间参数
        R1 = cos(phi)*cos(omiga)
        R2 = cos(phi)*sin(omiga)*sin(kappa)+sin(phi)*cos(kappa)
        R3 = sin(phi)*sin(kappa)-cos(phi)*sin(omiga)*cos(kappa)
        R4 = -sin(phi)*cos(omiga)
        R5 = cos(phi)*cos(kappa)-sin(phi)*sin(omiga)*sin(kappa)
        R6 = sin(phi)*sin(omiga)*cos(kappa) + cos(phi)*sin(kappa)
        R7 = sin(omiga)
        R8 = -cos(omiga)*sin(kappa)
        R9 = cos(omiga)*cos(kappa)

        # 旋转矩阵R
        R = np.array([
            [R1,R2,R3],
            [R4,R5,R6],
            [R7,R8,R9]
        ])


        # N
        N = R.dot(coords_S.T)

        N = N.reshape([num_p*3,1])


        # 单位矩阵I
        I= np.zeros([15,3])
        for ii in range(0,3):
            I[ii*5:ii*5+5,ii] = 1

        # R'
        R_ = np.hstack((I,n*M,N))

        
        # L 
        T = T.reshape([3,1])
        l = -T - (n*R.dot(coords_S.T))
        l = l.reshape([num_p*3,1])

        X_T = coords_T.T.reshape([num_p*3,1])

        L = l + X_T

        # V = R'x - (l + X-T)
        # x = (R'.T*R')(R'.T*L)
        p = np.linalg.inv(R_.T.dot(R_)).dot(R_.T.dot(L))

        # 精度评价
        V = R_.dot(p) - (L)
        # 自由度f = 3n - 7,n为使用的公共点个数
        f = 3*num_p - 7
        sg = np.sqrt(V.T.dot(V) / f)

        return p,sg
    
    
    @classmethod
    def bursa_Paramenter(cls,coords_S:np.array,coords_T:np.array):
        '''
        计算布尔莎模型参数的函数
        
        Parameter:
            coords_S  - 原坐标系控制点
            coords_T  - 对应目标坐标系控制点
            
        Return:
            T  - 平移参数向量，3 x 1
            phi  - 绕Z轴的旋转角φ
            omiga  - 绕Y轴的旋转角Ω
            kappa  - 绕X轴的旋转角κ
            n  - 尺度缩放，n = m + 1
            sg  - 模型拟合最终的精度评价指标，来源于Bursa.adjustment()函数
            
        '''
        T = np.array([0.0,0.0,0.0])
        phi = 0
        omiga = 0
        kappa = 0
        n = 1
        
        x,sg = cls.adjustment(T,phi,omiga,kappa,n,coords_S,coords_T)
        while (np.max(abs(x)) >= 1.0e-8):
            T += x[0:3].reshape([3])
            phi += float(x[3])
            omiga += float(x[4])
            kappa += float(x[5])
            n += float(x[6])
            x,sg= cls.adjustment(T,phi,omiga,kappa,n,coords_S,coords_T)
            
        T += x[0:3].reshape([3])
        phi += float(x[3])
        omiga += float(x[4])
        kappa += float(x[5])
        n += float(x[6])
        
        sg = round(float(sg),4)
        
        return T,phi,omiga,kappa,n,sg
    
    @staticmethod
    def accuracy_Evaluation(coords_S,coords_T,X_T):
        '''
        利用转换结果进行精度评价的函数。若利用求七参数的点进行评价则为內符精度，利用求七参数之外的点进行评价则为外符精度。
        Paramenter：
            coords_S  - 公共点转换前的点矩阵，(m * 3)  coords_S = [Xs  Ys  Zs].T
            X_T  - 经程序转换后的点位坐标，m * 3
            coords_T  - 公共点对应目标坐标系控制点坐标矩阵
        
        Return:
            v  - 各个点的X，Y，Z坐标残差
            sg  - 点位误差中误差绝对值
            sig  - X，Y，Z三个坐标的残差中误差的绝对值
            sig_b  - 各个点位的位置残差
        '''
        # 总点数
        n = coords_T.shape[0]
        
        # 残差
        v = X_T - coords_T
        
        # X，Y，Z坐标的残差中误差绝对值
        sig = np.sqrt(np.sum(v**2,axis = 0) / (n-1))
        
        # 点位中误差绝对值
        sg = np.sqrt(np.sum(sig**2))
        
        # 各个点位的位置残差
        sig_b  = np.sqrt(np.sum(np.power(v,2),axis = 1))
        
        return [sg,v,sig,sig_b]
    
    