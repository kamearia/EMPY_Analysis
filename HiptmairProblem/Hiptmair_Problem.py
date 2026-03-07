from netgen.occ import *
from ngsolve import *
import math
import sys
sys.path.append(r'../include')
from MatrixSolver import MatrixSolver as solver 
cpp_solver="EMPY"
#cpp_solver="JP_MARs"

class Hiptmair_Problem():
    import math
    c = 299792458.
    mu = 4*math.pi*1e-7
    eps = 1/(c*c*mu)
    
    def __init__(self, **kwargs):
        default_values = {"a":0.01,
                          "b":0.1,
                          "h":0.05,
                          "sigma":1.e6,
                          "jomega":True,
                          "SI":False,
                          "Darwin":False,
                          "Augumented":True
                         }
        default_values.update(kwargs)
        self.a=default_values["a"]
        self.b=default_values["b"]
        self.h=default_values["h"]
        self.sig=default_values["sigma"]
        self.jomega=default_values["jomega"]
        self.SI=default_values["SI"]
        self.Darwin=default_values["Darwin"]
        self.Augumented=default_values["Augumented"]

    def PrepareFESpace(self, feOrder=2):
        """FESpaceを一度だけ作成して保持する"""
        SI=self.SI
        Darwin=self.Darwin
        Augumented=self.Augumented
        mesh = self.mesh
        jomega = self.jomega
        self.SI=SI
        
        self.fesA = HCurl(mesh, order=feOrder, nograds=True, dirichlet="out|in|upper|lower|outer", complex=jomega)
        self.fesPhi = H1(mesh, order=feOrder, dirichlet="out|in|upper|lower", complex=jomega)
        
        if not Darwin or not Augumented:
            self.fesAPhi = self.fesA * self.fesPhi
        else:
            self.fesdPhi = H1(mesh, order=feOrder, dirichlet="out|in|upper|lower", complex=jomega)
            self.fesAPhi = self.fesA * self.fesPhi * self.fesdPhi
            
        self.gfAPhi = GridFunction(self.fesAPhi)
        # 行列計算用の線形形式などもここで定義可能

    def Z_calc(self, freq):
        print("frequency=", freq)
        import gc
        SI=self.SI
        Darwin=self.Darwin
        Augumented=self.Augumented
        # 関数内で fesAPhi = ... とせず、self.fesAPhi を使う
        fesAPhi = self.fesAPhi
        fesPhi=self.fesPhi
        gfAPhi = self.gfAPhi
        
        c = 299792458.
        mu = 4*math.pi*1e-7
        eps = 1/(c*c*mu)
        
        jomega=self.jomega
        mesh=self.mesh
        sigma = self.sigma
        
        if jomega:
            s=2j*math.pi*freq
        else:
            s=2*math.pi*freq

        self.s=s
        self.freq=freq
        
        if not Darwin or not Augumented:
            (A,phi), (N, psi) = fesAPhi.TnT()
            gfA, gfPhi=gfAPhi.components
        else:
            (A,phi, dphi), (N, psi, dpsi) = fesAPhi.TnT()
            gfA, gfPhi, gfdPhi=gfAPhi.components
    
        a=BilinearForm(fesAPhi)
        a += 1/mu*curl(A)*curl(N)*dx
        if not SI:
            a += sigma*s*(A+grad(phi))*(N+grad(psi))*dx("sig")
        else:
            Y=sqrt(self.sig/mu/s)
            a += (N.Trace()+grad(psi).Trace())* Y*( s*(A.Trace()+grad(phi).Trace()) ) *ds("conductorBND")
            
        if not Darwin:
            a += eps*s*s*(A+grad(phi))*(N+grad(psi))*dx
        else:
            a += eps*s*s*( grad(phi)*N + A*grad(psi) +grad(phi)*grad(psi))*dx
            if Augumented:
                a += -eps*s*s*( grad(dphi)*N + A*grad(dpsi)+ grad(dphi)*grad(dpsi))*dx

        with TaskManager():
            a.Assemble()

        #1V/m
        Vin=self.h/s
        #gfAPhi=GridFunction(fesAPhi)
        #gfA, gfPhi=gfAPhi.components
        #gfPhi.Set(Vin, BND, mesh.Boundaries("in|upper"))
        #Draw(1./Vin*gfPhi, mesh)
        gfPhi.Set(Vin, BND, mesh.Boundaries("in"))
        f=LinearForm(fesAPhi)
        r=f.vec-a.mat*gfAPhi.vec

        logout=True
        gfAPhi, log=solver.iccg_solve(fesAPhi, gfAPhi, a, r.Evaluate(), tol=1.e-20, max_iter=2000,
                         diviter=10, divfac=10., accel_factor=0,  complex=jomega, cpp_solver=cpp_solver, logplot=True, logout=logout, 
                        scaling=False)
        del a, f, r
        #gc.collect()
        #R0=1/crosssectionsel/self.sig
        #print("R0=", R0)

        A_PHI=A+grad(gfPhi)
        self.E=-s*(gfA+grad(gfPhi))
        J=-s*sigma*(gfA+grad(gfPhi))
        self.B=curl(gfA)

        if not SI:
            #fes = H1(mesh, order=1,  definedon="sig", dirichlet="out")
            fes = fesPhi
            u,v=fes.TnT()
            gf=GridFunction(fes)
            gf.Set(1,definedon= mesh.Boundaries("out|upper"))
            current=-Integrate(J*grad(gf)*dx("sig"), mesh)
            #current +=Integrate(eps*s*self.E*grad(gf)*dx, mesh)
        else:
            #fes = H1(mesh, order=1,  dirichlet="upper")
            fes = fesPhi
            u,v=fes.TnT()
            gf=GridFunction(fes)
            gf.Set(1,definedon= mesh.Boundaries("upper"))
            current=Integrate(Y*( -s*( gfA.Trace()+grad(gfPhi).Trace() ) ) *grad(gf)*ds("conductorBND"), mesh)
            current +=Integrate(eps*s*self.E*grad(gf)*dx, mesh)
        del gf
            
        gc.collect()
        #print("current=", -current)
        impedance=-1/current
        #print("impedance=", impedance)
        if logout:
            return impedance, log
        return impedance


    def Mesh(self, SI=False,**kwargs):
        import netgen.meshing as ngmeshing
        default_values = {"curveOrder":2, }
        default_values.update(kwargs)
        self.curveOrder=default_values["curveOrder"]
        SI=self.SI
        h=self.h
        r=self.b
        a=self.a
        curveOrder=self.curveOrder

        conductor = Cylinder((0,0,0), Z, r=0.01, h=h)
        conductor.faces.name="conductorBND"
        conductor.faces.Max(Z).name="in"
        conductor.faces.Min(Z).name="out"
        conductor.mat("sig")
        conductor.maxh=a/5
        crosssection = conductor.faces.Max(Z).mass

        outer = Cylinder((0,0,0), Z, r=0.1, h=h)
        outer.faces.name = "outer"

        air = outer-conductor
        air.faces.Max(Z).name="upper"
        air.faces.Min(Z).name="lower"
        #air.faces.col = (0,0,1,0.3)
        air.mat("air")
        air.maxh=r/5

        if not SI:
            model=Glue([conductor,air])
        else:
            model=air
            model.faces[3].name="conductorBND"           
        #model=Glue([conductor,air, integration_disc])
        #DrawGeo(model)
        geo =OCCGeometry(model)
        with TaskManager():
            #mesh = Mesh(geo.GenerateMesh(meshsize.coarse, quad_dominated=True)).Curve(1)
            #mesh = Mesh(geo.GenerateMesh(meshsize.fine)).Curve(1)
            #mesh = Mesh(geo.GenerateMesh(meshsize.fine)).Curve(1) 
            mesh = Mesh(geo.GenerateMesh(maxh=a/2)).Curve(1)

        self.sigma = mesh.MaterialCF({ "sig" : self.sig }, default=0)
        self.mesh=mesh
        return model, mesh
        
    def Z_analytical(self):
        from sympy import symbols, sqrt, solve, simplify, diff, pi
        from sympy.functions.special.bessel import besselk, besseli
        """
        MathematicaのZfull[s]関数をSymPyを使用してPythonで再現します。
        同軸ケーブルなどの円筒座標系におけるラプラス領域のインピーダンスを計算する式です。
        """

        # 1. 記号の定義 (Mathematicaの変数をSymPyの記号として定義)
        s, mu, eps, sigma, a, b, es = symbols('s mu eps sigma a b es')

        # 2. パラメータ g, g0 の定義
        g = sqrt(mu * s * (s * eps + sigma))
        g0 = sqrt(mu * eps * s**2) # s^2 のルートは s になるが、SymPyではそのまま記述

        # 3. 磁位 Az1(r), Az2(r) の定義
        # 積分定数 A0, A1, A2 はシンボルとして定義
        A0, A1, A2, r = symbols('A0 A1 A2 r')

        # r 依存の Az1, Az2 の式
        # Az1: r <= a (内側の導体と誘電体)
        Az1_r = A0 * besseli(0, g * r)
        # Az2: r > a (誘電体と外側の導体)
        Az2_r = A1 * besseli(0, g0 * r) + A2 * besselk(0, g0 * r)

        # 4. Az1' と Az2' の定義 (rに関する微分)
        # Mathematica の D[f[r], r] に相当
        dAz1_dr = diff(Az1_r, r)
        dAz2_dr = diff(Az2_r, r)

        # 5. 境界条件に基づく連立方程式の定義
        # Mathematica の Solve[{...}, {A0, A1, A2}] に相当

        # 境界条件 1: r = a で Az の連続性
        eq1 = Az1_r.subs(r, a) - Az2_r.subs(r, a)

        # 境界条件 2: r = a で Az' の連続性 (接線方向の E と H の連続性から導出)
        eq2 = dAz1_dr.subs(r, a) - dAz2_dr.subs(r, a)

        # 境界条件 3: r = b での Az の条件 (-s Az2[b] = es)
        # これは r=b での Ez の条件 Ez(b) = es に相当
        # Ez = -s Az だから、Ez(b) = -s Az2(b) = es
        eq3 = -s * Az2_r.subs(r, b) - es

        # 6. A0, A1, A2 について連立方程式を解く
        # 解は辞書のリストとして返される
        solution = solve([eq1, eq2, eq3], (A0, A1, A2))
    
        # solve の結果が空でないことを確認
        if not solution:
            return "SymPy solve could not find a solution for A0, A1, A2."

        # 7. Ht2(r) の定義 (H_phi の定義)
        # Ht2[r] = - Az2'[r]/mu (Phi方向の磁界 H_phi)
        Ht2_r = - dAz2_dr / mu
    
        # 8. 最終的なインピーダンス Z の式を計算
        # Z = -s Az2[b] / (Ht2[b] * 2 * Pi * b)
        # Z = Ez[b] / (H_phi[b] * 2 * Pi * b)
        # Ez[b] = es
    
        # 分子の計算: -s Az2[b] は境界条件3から es と等しい
        numerator = es

        # 分母の計算: Ht2[b] * 2 * Pi * b
        # Ht2[b] は r=b での Ht2_r に A1, A2 の解を代入
        Ht2_b = Ht2_r.subs(r, b).subs(solution)
        denominator = Ht2_b * 2 * pi * b

        # 9. インピーダンス Z の計算と簡略化
        Z_full_s = simplify(numerator / denominator)

        return Z_full_s
    """
    def Z_numerical(self, s_value, xp=None):
        import numpy as np
        import scipy.special as sp
        import math

        a, b, sigma = self.a, self.b, self.sig
        c = 299792458.
        mu = 4.0e-7 * math.pi
        eps = 1.0 / (c * c * mu)
        s = complex(s_value)
        if s == 0: return 0

        # パラメータ計算
        g = np.sqrt(mu * s * (s * eps + sigma))
        g0 = s / c  # 誘電体中の伝搬定数 (j*omega/c)
        ga, g0a, g0b = g * a, g0 * a, g0 * b

        # --- 改善ポイント: スケーリング済みベッセル関数を利用 ---
        # I1(ga)/I0(ga) は高周波で1に収束する。iveを使えばオーバーフローしない。
        ratio_I = sp.ive(1, ga) / sp.ive(0, ga)

        # 誘電体領域 (r=a, r=b) のベッセル関数
        # 誘電体側は引数が小さいため、通常の iv, kv で安定
        I0g0a = sp.iv(0, g0a); K0g0a = sp.kv(0, g0a)
        I1g0a = sp.iv(1, g0a); K1g0a = sp.kv(1, g0a)
        I0g0b = sp.iv(0, g0b); K0g0b = sp.kv(0, g0b)
        K1g0b = sp.kv(1, g0b)

        # 4. 正規化された連立方程式 M_scaled * [A0_tilde, A1, A2] = V
        # A0_tilde = A0 * I0(ga) と置くことで、巨大な I0(ga) を行列から消去
        M1 = [1.0, -I0g0a, -K0g0a]
        M2 = [g * ratio_I, -g0 * I1g0a, g0 * K1g0a]
        M3 = [0.0, I0g0b, K0g0b]

        M = np.array([M1, M2, M3], dtype=complex)
        es = 1.0
        V = np.array([0.0, 0.0, -es / s], dtype=complex)

        try:
            sol = np.linalg.solve(M, V)
            A1, A2 = sol[1], sol[2]

            # Ht2(b) の計算
            Az2_prime_b = A1 * g0 * sp.iv(1, g0b) - A2 * g0 * sp.kv(1, g0b)
            Ht2_b = - Az2_prime_b / mu
    """
    def Z_numerical(self, s_value, xp=None):
        import numpy as np
        import scipy.special as sp
        import math

        a_param, b_param, sigma = self.a, self.b, self.sig
        c = 299792458.
        mu = 4.0e-7 * math.pi
        eps = 1.0 / (c * c * mu)
        s = complex(s_value)
        if s == 0: return 0

        # 1. パラメータ計算
        g = np.sqrt(mu * s * (s * eps + sigma))
        g0 = s / c  
        ga, g0a, g0b = g * a_param, g0 * a_param, g0 * b_param

        # 2. スケーリング済みベッセル関数で導体部分の比率を計算 (2e8 Hz以上の対策)
        ratio_I = sp.ive(1, ga) / sp.ive(0, ga)

        # 3. 誘電体領域のベッセル関数 (r=a, r=b)
        I0g0a, K0g0a = sp.iv(0, g0a), sp.kv(0, g0a)
        I1g0a, K1g0a = sp.iv(1, g0a), sp.kv(1, g0a)
        I0g0b, K0g0b = sp.iv(0, g0b), sp.kv(0, g0b)
        K1g0b = sp.kv(1, g0b)

        # 4. 行列構築と求解
        M = np.array([
            [1.0, -I0g0a, -K0g0a],
            [g * ratio_I, -g0 * I1g0a, g0 * K1g0a],
            [0.0, I0g0b, K0g0b]
        ], dtype=complex)

        es = 1.0
        V = np.array([0.0, 0.0, -es / s], dtype=complex)

        try:
            sol = np.linalg.solve(M, V)
            A1, A2 = sol[1], sol[2]

            # --- 5. 任意地点 xp での電磁界計算 (誘電体中 a <= r <= b) ---
            E_out = []
            B_out = []
            if xp is not None:
                for x in xp:
                    A0_tilde = sol[0]
                    if x < a_param:
                        # --- 導体内部 (r < a) ---
                        # Az(r) = A0 * I0(gx) = A0_tilde * (I0(gx) / I0(ga))
                        # 比率で計算することでオーバーフローを防止
                        ratio_Az = sp.ive(0, g*x) / sp.ive(0, g*a_param) * np.exp(np.real(g)*(x - a_param))
                        Ez_r = -s * A0_tilde * ratio_Az
                        
                        # dAz/dr = A0 * g * I1(gx) = A0_tilde * g * (I1(gx) / I0(ga))
                        ratio_dAz = sp.ive(1, g*x) / sp.ive(0, g*a_param) * np.exp(np.real(g)*(x - a_param))
                        dAz_dr = A0_tilde * g * ratio_dAz
                        B_phi = -dAz_dr # B = curl A より B_phi = -dAz/dr

                    else:
                        # Ez(r) = -s * (A1 * I0(g0*r) + A2 * K0(g0*r))
                        Az_r = A1 * sp.iv(0, g0*x) + A2 * sp.kv(0, g0*x)
                        Ez_r = -s * Az_r

                    
                        # H_phi(r) = - (1/mu) * dAz/dr
                        # dAz/dr = g0 * (A1 * I1(g0*r) - A2 * K1(g0*r))
                        dAz_dr = g0 * (A1 * sp.iv(1, g0*x) - A2 * sp.kv(1, g0*x))
                        B_phi = -dAz_dr
                        # B = mu * H
                    
                    E_out.append(Ez_r)
                    B_out.append(B_phi)

            # 6. インピーダンス計算用の磁界 (r=b)
            Az2_prime_b = g0 * (A1 * sp.iv(1, g0b) - A2 * sp.kv(1, g0b))
            Ht2_b = - Az2_prime_b / mu
            impedance = es / (Ht2_b * 2.0 * math.pi * b_param)
            
            if xp is not None:
                return impedance, E_out, B_out
            return impedance

        except (np.linalg.LinAlgError, RuntimeWarning):
            return np.nan
    
    def Zsi_numerical(self, s_value): 
        import math
        import numpy as np
        import scipy.special as sp
        
        a = self.a
        b = self.b
        sigma = self.sig
        c = 299792458.
        mu = 4.0e-7 * math.pi
        eps = 1.0 / (c * c * mu)
        
        s = complex(s_value)
        if s == 0: return 0

        # 1. パラメータ計算 (sとの積でg0を定義し、枝切りの問題を回避)
        g0 = s * np.sqrt(mu * eps) 
        Zs = np.sqrt(s * mu / sigma)  

        g0a = g0 * a
        g0b = g0 * b

        # 2. ベッセル関数の計算 (kveを使用して安定化)
        I0g0a = sp.iv(0, g0a); I1g0a = sp.iv(1, g0a)
        I0g0b = sp.iv(0, g0b); I1g0b = sp.iv(1, g0b)
        
        K0g0a = sp.kve(0, g0a) * np.exp(-g0a); K1g0a = sp.kve(1, g0a) * np.exp(-g0a)
        K0g0b = sp.kve(0, g0b) * np.exp(-g0b); K1g0b = sp.kve(1, g0b) * np.exp(-g0b)

        # 3. 連立方程式の構築 (符号を修正: Ez = Zs * Ht)
        # 方程式1 (r=a): -s*Az2(a) + (Zs/mu)*Az2'(a) = 0
        # Az2'(a) = g0 * (A1*I1g0a - A2*K1g0a)
        # -> A1 * (-s*I0g0a + Zs*g0/mu*I1g0a) + A2 * (-s*K0g0a - Zs*g0/mu*K1g0a) = 0
        
        coeff_A1_eq1 = -s * I0g0a + (Zs * g0 / mu) * I1g0a # 符号修正
        coeff_A2_eq1 = -s * K0g0a - (Zs * g0 / mu) * K1g0a # 符号修正
    
        # 方程式2 (r=b): -s*Az2(b) = es
        coeff_A1_eq2 = -s * I0g0b
        coeff_A2_eq2 = -s * K0g0b

        M = np.array([[coeff_A1_eq1, coeff_A2_eq1],
                      [coeff_A1_eq2, coeff_A2_eq2]], dtype=complex)
        es = 1.0
        V = np.array([0, es], dtype=complex)

        try:
            A_sol = np.linalg.solve(M, V)
            A1, A2 = A_sol[0], A_sol[1]

            # Ht2(b) = -Az2'(b) / mu
            Ht2_b = -(g0 / mu) * (A1 * I1g0b - A2 * K1g0b)
        
            # インピーダンス Z = es / (Ht2_b * 2 * pi * b)
            return es / (Ht2_b * 2.0 * math.pi * b)

        except np.linalg.LinAlgError:
            return np.nan

    def PlotField(self):
        import matplotlib.pylab as plt
        import numpy as np
        mesh=self.mesh
        h=self.h
        b=self.b
        a=self.a
        n=1009
        if self.SI:
            d=(b-a)/n
            x=a
        else:
            d=b/n
            x=0
        y=0
        z=h/2
        xp=[]
        ezr=[]
        ezi=[]
        byr=[]
        byi=[]
        for n in range(n+1):
            pnt=mesh(x,y,z)
            xp.append(x)
            ezr.append(-self.E(pnt)[2].real)
            ezi.append(-self.E(pnt)[2].imag)
            byr.append(-self.B(pnt)[1].real)
            byi.append(-self.B(pnt)[1].imag)
            x=x+d

        impedance, E, B = self.Z_numerical(self.s, xp)
        Er=np.real(E)
        Ei=np.imag(E)
        Br=np.real(B)
        Bi=np.imag(B)
        

        plt.plot(xp, ezr,   linestyle='-',  color='Red',  linewidth=2,  label='Re[Edarwin]') 
        plt.plot(xp, Er,    linestyle=':',  color='Red',  linewidth=2,  label='Re[Ez]') 
        plt.xlabel("x")  # Add x-axis label
        plt.plot(xp, ezi,  linestyle='-',  color='black',  linewidth=2,  label='Im[Edarwin]')
        plt.plot(xp, Ei,    linestyle=':',  color='black',  linewidth=2,  label='Im[Ez]')
        plt.xlabel("x")  # Add x-axis label
        plt.ylabel("Ez")  # Add y-axis label
        plt.xlim(0, b)
        if self.Darwin: 
            plt.title(f"E field  ( Darwin model ) Frequency {self.freq/1.e6} MHz " )
        else:
             plt.title(f"E field  ( Full Maxwell ) Frequency {self.freq/1.e6} MHz " )           
        plt.legend()
        #xplt.ylim(0, 0.2)
        plt.show() 

        plt.plot(xp, byr,   linestyle='-',  color='Red',  linewidth=2,  label='Re[Bdarwin]')   
        plt.plot(xp, Br,    linestyle=':',  color='Red',  linewidth=2,  label='Re[Bt]') 
        plt.xlabel("x")  # Add x-axis label
        plt.plot(xp, byi ,  linestyle='-',  color='black',  linewidth=2,  label='Im[Bdarwin]')
        plt.plot(xp, Bi ,    linestyle=':',  color='black',  linewidth=2,  label='Im[Bt]')
        plt.xlabel("x")  # Add x-axis label
        plt.ylabel("By")  # Add y-axis label
        plt.xlim(0, b)
        #plt.ylim(-1e-5, 7e-5)
        if self.Darwin: 
            plt.title(f"B field  ( Darwin model ) Frequency {self.freq/1.e6} MHz " )
        else:
             plt.title(f"B field  ( Full Maxwell ) Frequency {self.freq/1.e6} MHz " )      
        plt.legend()
        plt.show()  
