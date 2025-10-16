import math
from netgen.occ import *
from ngsolve import *
from Prob7Coil import Prob7Coil
sys.path.append(r'..\COIL\include')
from EMPY_COIL import *
sys.path.append(r'..\bin\Release') 
from EMPY_Field import *
from ngsolve.webgui import Draw
from netgen.webgui import Draw as DrawGeo 
class Team7Model():
    def __init__(self,  **kwargs):
        default_values = {"name":"Team Workshop Problem #7 Model",
                          "sigma": 3.526e7, 
                          "msize": meshsize.moderate,
                          "boxsize": 1.0,
                          "rKelvin":0
                         }
        default_values.update(kwargs)
        sigma=default_values["sigma"]
        msize=default_values["msize"]
        boxsize=default_values["boxsize"]
        rKelvin=default_values["rKelvin"]  
        self.rKelvin=rKelvin
        
        xl=-1.353
        xu=1.637
        yl=-1.357
        yu=1.647
        zl=-0.3
        zu=0.449
        self.curveOrder=3

        pc0=gp_Pnt(0,0,0)
        pc1=gp_Pnt(294.e-3,294.e-3, 19.e-3)
        conductor=Box(pc0,pc1)
        center=gp_Pnt(0.5*(pc0[0]+pc1[0]),0.5*(pc0[1]+pc1[1]),0.5*(pc0[2]+pc1[2]))
        ph0=(18e-3,18e-3,0)
        ph1=(126.e-3,126.e-3, 19.e-3)
        hole=Box(ph0,ph1)
        hole.mat("hole")
        hole.faces.col = (0,0,0,0.0)
        conductor=conductor-hole
        conductor.mat("conductor")
        conductor.faces.name="conductor_boundary"
        conductor.faces.col = (1,0,0,1)
        
        #d=15e-3
        #d=30e-3
        d=22.5e-3
        pt0=pc0-gp_Vec(d,d,d)
        pt1=pc1+gp_Vec(d,d,d)
        total=Box(pt0,pt1)
        total.faces.name="total_boundary"
        #total.faces.col = (1,1,0,0.5)

        air=total-conductor-hole
        air.mat("air")
        #air.faces.col = (1,1,0,0.5)

        if rKelvin ==0:
            if boxsize ==0:
                p0=gp_Pnt(xl,yl,zl)
                p1=gp_Pnt(xu,yu,zu)
            else:
                d=boxsize
                p0=center-gp_Vec(d,d,d)
                p1=center+gp_Vec(d,d,d)
            reduced=Box(p0,p1)
        else:
            reduced=Sphere(center,rKelvin)
    
        reduced.faces.name="reduced_boundary"
        reduced=reduced-total
        reduced.mat("reduced_region")
        reduced.faces.col = (0,0,1,0.2)

        if rKelvin !=0:
            rk=rKelvin
            kcenter=rk*2.5+center[0]
            self.kcenter=kcenter
            external_domain = Sphere(Pnt(kcenter,0,0.), r=rk)
            external_domain.mat("Kelvin")
            external_domain.faces.col = (0,0,1,0.2)
            external_domain.faces[0].Identify(reduced.faces[0], "ud0",  IdentificationType.PERIODIC)
            geo=Glue([conductor, hole, air, reduced, external_domain])
        else:
            geo=Glue([conductor, hole, air, reduced])

        self.geo=geo

        with TaskManager():
            #mesh = Mesh(geo.GenerateMesh(meshsize.coarse, quad_dominated=True)).Curve(1)
            #mesh = Mesh(geo.GenerateMesh(meshsize.fine)).Curve(1)
            #mesh = Mesh(geo.GenerateMesh(meshsize.moderate)).Curve(curveOrder) 
            mesh = Mesh(OCCGeometry(geo).GenerateMesh(msize)).Curve(self.curveOrder) 
        self.mesh=mesh

        sigma_d={"conductor":sigma, "hole":0,  "air":0, "reduced_region":0, "Kelvin":0, 'default':0}
        rho_d={"conductor":1./sigma, "hole":0,  "air":0, "reduced_region":0, "Kelvin":0,'default':0}
        mu = 4*math.pi*1e-7
        mu_d={"conductor":mu, "hole":mu,  "air":mu, "reduced_region":mu, "Kelvin":mu,'default':mu}
        self.Sigma = CoefficientFunction([sigma_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        self.Rho = CoefficientFunction([rho_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        self.Mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性
        
        self.conductor_boundary="conductor_boundary"
        self.conductive_region="conductor"
        self.symmetric_plane="None"
        self.reduced_region="reduced_region"
        self.air_region="hole|air|reduced_region" 
        self.total_region="conductor|hole|air"
        self.total_boundary="total_boundary"
        self.reduced_boundary="reduced_boundary"
        self.total_air_region="air"
        self.Bn0_boundary=""
        self.Ht0_boundary=""
        
        #self.coil=EMPY_UNIF(0,0,1,0)
        self.coil=Prob7Coil()
        self.coil.geo.faces.col = (0,1,0,1)
        if self.coil.geo :
            self.geo=geo-self.coil.geo
            self.geo=Glue([self.geo,self.coil.geo])

 

        self.conductor=conductor

    def ReducedField(self, Bv, zero):
        Bs_dic = { "reduced_region":Bv,  "air": zero, "total": zero, 'default':zero,
                  "conductor":zero, "hole":zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in self.mesh.GetMaterials()])
        return Bs
    def AppliedField(self):
        Bv=Bfield(self.coil.field)
        mesh = Mesh(OCCGeometry(self.conductor).GenerateMesh(meshsize.coarse)).Curve(self.curveOrder) 
        with TaskManager():
            Draw(Bv, mesh, order=3)

    def PlotResults(self, mesh, BField, JField, freq):
        self.PlotBField(mesh, BField, freq)
        self.PlotJField(mesh, JField, freq)

    def PlotBField(self, mesh, BField, freq):
        import numpy as np
        import matplotlib.pylab as plt
        
        xi_msm = np.array([0, 18, 36, 54, 72, 90, 108, 126, 144, 162, 180, 198, 216, 234, 252, 270, 288]) * 1e-3
        xi_sim = np.linspace(0,0.288,500)
        Bz_A1B1_50Hz_ref = np.array([ -4.9  -1.16j, -17.88 +2.48j, -22.13 +4.15j, -20.19 +4.j,   -15.67 +3.07j,
           0.36 +2.31j,  43.64 +1.89j,  78.11 +4.97j,  71.55+12.61j , 60.44+14.15j,
          53.91+13.04j,  52.62+12.4j,  53.81+12.05j, 56.91+12.27j,  59.24+12.66j,
          52.78 +9.96j,  27.61 +2.26j])
        Bz_A1B1_200Hz_ref = np.array([ -3.63  -1.38j, -18.46 +1.20j, -23.62 +2.15j, -21.59 +1.73j,   -16.09 +1.10j,
           0.23 +0.27j,  44.35 -2.28j,  75.53 -1.40j,  63.42+4.17j , 53.20+ 3.94j,
          48.66+ 4.86j,  47.31+ 4.09j,  48.31+ 3.69j, 51.26+ 4.60j,  53.61+ 3.48j,
          46.11 +4.10j,  24.96 +0.98j])
               
        Bz_A2B2_50Hz_ref = np.array([-1.83-1.63j, -8.5-0.6j, -13.6-0.43j, -15.21+0.11j, -14.48+1.26j, -5.62+3.4j,
             28.77+6.53j, 60.34+10.25j, 61.84+11.83j, 56.64+11.83j, 53.4+11.01j, 52.36+10.58j, 53.93+10.8j, 56.82+10.54j, 
             59.48+10.62j, 52.08+9.03j, 26.56+1.79j])

        Bz_A2B2_200Hz_ref = np.array([-0.86-1.35j, -7.00-0.71j, -11.58-0.81j, -13.36+0.67j, -13.77+0.15j, -6.74+1.39j,
             24.63+2.67j, 53.19+ 3.00j, 54.89+4.01j, 50.72+3.80j, 48.03+4.00j, 47.13+3.02j, 48.25+2.20j, 51.35+2.78j, 
             53.35+1.58j, 45.37+1.37j, 24.01+0.93j])

        if freq==50:
            Bz_A1B1_ref=Bz_A1B1_50Hz_ref
            Bz_A2B2_ref=Bz_A2B2_50Hz_ref
        elif freq==200:
            Bz_A1B1_ref=Bz_A1B1_200Hz_ref
            Bz_A2B2_ref=Bz_A2B2_200Hz_ref
            
            
        Bz_A1B1_sim = np.array([BField[2](mesh(x, 0.072, 0.034)) for x in xi_sim])
        plt.figure(1)
        plt.clf()
        plt.title("Evaluate B.z on Line A1 - B1")
        plt.plot(xi_sim, 1e4*Bz_A1B1_sim.real, "b", label="sim Bz.real")
        plt.plot(xi_sim, -1e4*Bz_A1B1_sim.imag, "r", label="sim Bz.imag")
        plt.plot(xi_msm, Bz_A1B1_ref.real, "x", label="ref Bz.real")
        plt.plot(xi_msm, Bz_A1B1_ref.imag, "x", label="ref Bz.imag") 
        plt.xlabel("x")
        plt.ylabel("B in Gauss = $10^{-4}$ T")
        plt.legend()
        plt.show() 

        Bz_A2B2_sim = np.array([BField[2](mesh(x, 0.144, 0.034)) for x in xi_sim])
        plt.figure(2)
        plt.clf()
        plt.title("Evaluate B.z on Line A2 - B2")
        plt.plot(xi_sim, 1e4 * Bz_A2B2_sim.real, "b", label="sim Bz.real")
        plt.plot(xi_sim, -1e4 * Bz_A2B2_sim.imag, "r", label="sim Bz.imag")
        plt.plot(xi_msm, Bz_A2B2_ref.real, "b--x", label="ref Bz.real")
        plt.plot(xi_msm, Bz_A2B2_ref.imag, "r--x", label="ref Bz.imag")
        plt.xlabel("x")
        plt.ylabel("B in Gauss = $10^{-4}$ T")
        plt.legend()
        plt.show()

    def PlotJField(self, mesh, JField, freq):
        import numpy as np
        import matplotlib.pylab as plt

        xi_msm = np.array([0, 18, 36, 54, 72, 90, 108, 126, 144, 162, 180, 198, 216, 234, 252, 270, 288]) * 1e-3
        xi_sim = np.linspace(0,0.288,500)

        Jy_A3B3_50Hz_ref = np.array([0.249-0.629j,  0.685-0.873j,  0+0j,  0+0j,  0+0j,  0+0j,  0+0j,  -0.015-0.593j,
             -0.103-0.249j,  -0.061-0.101j,  -0.004-0.001j,  0.051+0.087j,  0.095+0.182j,  0.135+0.322j,  
             0.104+0.555j,  -0.321+0.822j,  -0.687+0.855j,])
        
        Jy_A3B3_200Hz_ref = np.array([0.427-0.623j,  0.794-0.755j,  0+0j,  0+0j,  0+0j,  0+0j,  0+0j,  1.401-1.304j,
             -0.035-0.229j,  0.005-0.041j,  -0.011-0.014j,  0.007-0.002j,  0.027-0.000j,  0.042+0.008j,  
             0.043+0.033j,  0.050+0.116j,  -0.321+0.893j,])

        Jy_A4B4_50Hz_ref = np.array([0.461-0.662j,  0.621-0.664j,  0+0j,  0+0j,  0+0j,  0+0j,  0+0j,  1.573-1.027j,
             0.556-0.757j,  0.237-0.364j,  0.097-0.149j,  -0.034+0.015j,  -0.157+0.154j,  -0.305+0.311j,
             -0.478+0.508j,  -0.66+0.747j,  -1.217+1.034j])

        Jy_A4B4_200Hz_ref = np.array([1.057-0.915j,  1.597-1.036j,  0+0j,  0+0j,  0+0j,  0+0j,  0+0j,  4.163-2.328j,
             1.143-1.193j,  0.672-0.613j,  0.307-0.259j,  -0.050+0.0061j,  -0.370+0.334j,  -0.749+0.674j,
             -1.205+1.064j,  -1.575+1.494j,  -2.583+2.331j])

        if freq==50:
            Jy_A3B3_ref=Jy_A3B3_50Hz_ref
            Jy_A4B4_ref=Jy_A4B4_50Hz_ref
        elif freq==200:
            Jy_A3B3_ref=Jy_A3B3_200Hz_ref
            Jy_A4B4_ref=Jy_A4B4_200Hz_ref
        
        Jy_A3B3_sim = np.array([1e-6*1j*JField[1](mesh(x, 0.072, 0.019-1e-5)) for x in xi_sim])
        plt.figure(3)
        plt.clf()
        plt.title("Evaluate J.y on Line A3 - B3")
        plt.plot(xi_sim, -Jy_A3B3_sim.real, "b", label="sim Jy.real")
        plt.plot(xi_sim, -Jy_A3B3_sim.imag, "r", label="sim Jy.imag")
        plt.plot(xi_msm, Jy_A4B4_ref.real, "bx", label="ref Jy.real")
        plt.plot(xi_msm, Jy_A4B4_ref.imag, "rx", label="ref Jy.imag")
        plt.xlabel("x")
        plt.ylabel("J in $10^{6}$ A/mm$^2$")
        #plt.ylim([-1.5,1.5])
        plt.legend()
        plt.show()

        Jy_A4B4_sim = np.array([1e-6*1j*JField[1](mesh(x, 0.072, 0.000+1e-5)) for x in xi_sim])
        plt.figure(4)
        plt.clf()
        plt.title("Evaluate J.y on Line A4 - B4")
        plt.plot(xi_sim, -Jy_A4B4_sim.real, "b", label="sim Jy.real")
        plt.plot(xi_sim, -Jy_A4B4_sim.imag, "r", label="sim Jy.imag")
        plt.plot(xi_msm, Jy_A3B3_ref.real, "bx", label="ref Jy.real")
        plt.plot(xi_msm, Jy_A3B3_ref.imag, "rx", label="ref Jy.imag")
        plt.xlabel("x")
        plt.ylabel("J in $10^{6}$ A/mm$^2$")
        #plt.ylim([-1.2,1.2])
        plt.legend()
        plt.show()
        plt.show()

