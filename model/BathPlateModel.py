# -*- coding: utf-8 -*-
from netgen.occ import *
from ngsolve import TaskManager
import math
from ngsolve import *
import sys
sys.path.append(r'..\COIL\include')
from Prob3Coil import Prob3Coil

from ngsolve import *
from netgen.occ import *
from ngsolve.webgui import Draw
from netgen.webgui import Draw as DrawGeo    

class MeasureFace():
    num=0
    faces=[]
    mesh=None
    def __init__(self, face, from_side, to_side):
        self.face=face
        self.from_side=from_side
        self.to_side=to_side
        MeasureFace.num=MeasureFace.num+1
        MeasureFace.faces.append(self)
    def Printm(self):
        print("    ", self.face, "  ",self.from_side, "  ", self.to_side) 
 
    @classmethod
    def Print(cls):
        print("class MeasureFace: num=", cls.num, " faces")
        for n in range(cls.num):
            cls.faces[n].Printm()
        
    @classmethod
    def SetMesh(cls, mesh):
        cls.mesh=mesh
    @classmethod
    def GetMesh(cls):
        return cls.mesh

    @classmethod
    def CalcFluxes(cls, field):
        fluxes=[]
        for n in range(cls.num):
            face=cls.faces[n]
            fluxes.append(face.CalcFlux(field) )
        return fluxes
               
    def CalcFlux(self, field):        
        mesh=self.mesh
        fes = H1(mesh, order=1,  definedon=self.from_side, dirichlet=self.face)
        u,v=fes.TnT()
        gf=GridFunction(fes)
        gf.Set(1,definedon= mesh.Boundaries(self.face))
        int1=Integrate(field*grad(gf)*dx(self.from_side), self.mesh)
        #print("int1=", int1)
        fes = H1(mesh, order=1,  definedon=self.to_side, dirichlet=self.face)
        u,v=fes.TnT()
        gf=GridFunction(fes)
        gf.Set(1,definedon= mesh.Boundaries(self.face))
        int2=Integrate(field*grad(gf)*dx(self.to_side), mesh)
        #print("int2=", int2)
        flux=(int1-int2)/2
        #print("int=", int)
        return flux
    
class HolePot():
    num=0
    holes=[]
    def __init__(self, name, lower, upper, interface):
        self.name=name
        self.lower=lower
        self.upper=upper
        self.interface=interface
        HolePot.num=HolePot.num+1
        HolePot.holes.append(self)
        
    def Printm(self):
        print("    ", "neme:", self.name, "  lower:", self.lower, "  upper:", self.upper, "  interface:", self.interface)
        
    @classmethod
    def Print(cls): 
        print("class HolePot: num=", cls.num, " holes")
        for n in range(cls.num):
            cls.holes[n].Printm()
    
    def SetInterface(self, mesh, gfT, order):
        fes=H1(mesh, order=order, definedon=self.name,dirichlet=self.lower+"|"+self.upper)
        gf=GridFunction(fes)
        gf.Set(1, definedon=mesh.Boundaries(self.upper)) 
        SetBoundaryValue(gf, 1, gfT, self.interface)
        self.gf=gf
        return gf
    
    @classmethod   
    def SetInterfaces(cls,mesh, gfT, order=1):
        for n in range(cls.num):
            cls.holes[n].SetInterface(mesh, gfT, order)
        return gfT
    
    @classmethod
    def GetGfs(cls):
        return gfs
    
class BathPlateModel():
    def __init__(self,  **kwargs):
        self.mu=None
        self.sigma=None
        self.outer_boundary=None
        self.conductiveDomain=None
        self.conductor_boundary=None
        self.coil=None
        self.SetMesh(**kwargs)
        
    def GetGeometry(self):
        return self.model
    
    def GetMesh(self):
        return self.mesh

    def Print(self):
        geo=self.GetGeometry()
        for s in geo.solids:
            print("name:",s.name, "  mass:", s.mass, "  center:", s.center)
        print("conductive_region:", self.conductive_region)
        print("Conductor boundary:", self.conductor_boundary)
        print(self.sigma_d)
        print(self.mu_d)
        HolePot.Print()
        MeasureFace.Print()

    def SetMesh(self, **kwargs):
        from netgen.occ import Box
        from netgen.occ import Z
        from netgen.occ import X
        from netgen.occ import Y
        from netgen.occ import Glue
        from netgen.occ import OCCGeometry
        from ngsolve.comp import Mesh
        from netgen.meshing import meshsize

        default_values = {"name":"Bath Plate Model",
                          "holes":2,
                          "outerBox":"sphere",
                          "wz":6.35e-3,
                          "boxx":255e-3, 
                          "boxy":280e-3,      
                          "boxz":400e-3, 
                          "div_thick":1,
                          "sigma": 0.3278e8, 
                          "msize": meshsize.moderate,
                          "rKelvin":0
                         }
        default_values.update(kwargs)
        holes=default_values["holes"]
        outerBox=default_values["outerBox"]
        bx=default_values["boxx"]
        by=default_values["boxy"]
        bz=default_values["boxz"]
        wz=default_values["wz"]
        div_thick=default_values["div_thick"]
        msize=default_values["msize"]
        rKelvin=default_values["rKelvin"]

        sig=default_values["sigma"]
        curveOrder=3
        self.curveOrder=curveOrder
        mu = 4*math.pi*1e-7
        self.up_down=2.*bz

        if rKelvin !=0:
            outerBox="Sphere"
            bx=rKelvin
            self.rKelvin=rKelvin

        print(default_values["name"], "   holes:", holes, "   Thickness:", wz)
        print("boxx= ", bx," boxy= ", by, " boxz= ",bz)
        print("div_thick= ", div_thick)
            

        wx=110e-3
        wy=60e-3
        wz=default_values["wz"]
        hx=30e-3
        hy=40e-3
        sx=20e-3
        conductor=Box((-wx/2,-wy/2,-wz/2),(wx/2,wy/2,wz/2))
        if holes !=0:
            hole1=Box((-hx/2+sx,-hy/2,-wz/2),(hx/2+sx,hy/2,wz/2))
            conductor=conductor-hole1
            if holes==2:
                hole2=Box((-hx/2-sx,-hy/2,-wz/2),(hx/2-sx,hy/2,wz/2))
                conductor=conductor-hole2
        conductor.mat("conductor")
        conductor.faces.name="conductor_boundary"
        
        d=7.5e-3
        total=Box((-wx/2-d,-wy/2-d,-wz/2-d),(wx/2+d,wy/2+d,wz/2+d))
        total.faces.name="total_boundary"
        if outerBox=="box":
            outer=Box((-bx,-by,-bz),(bx,by,bz))
            curveOrder=1
        else:outer=Sphere(Pnt(0,0,0),bx)
        outer.faces.name="reduced_boundary"
        reduced=outer#-total
        reduced.mat("reduced_region")
        air=total #-conductor
        air.mat("air")

        if rKelvin !=0:
            rk=rKelvin
            kcenter=rk*2.5
            self.kcenter=kcenter
            external_domain = Sphere(Pnt(kcenter,0,0.), r=rk)
            external_domain.mat("Kelvin")

            external_domain.faces[0].Identify(outer.faces[0], "ud0",  IdentificationType.PERIODIC)
            model=Glue([conductor,air, reduced, external_domain])
        else:
            model=Glue([conductor,air, reduced])
            
        self.conductor=conductor
        self.model=model

        geo =OCCGeometry(model)
        with TaskManager():
            #mesh = Mesh(geo.GenerateMesh(meshsize.coarse, quad_dominated=True)).Curve(1)
            #mesh = Mesh(geo.GenerateMesh(meshsize.fine)).Curve(1)
            #mesh = Mesh(geo.GenerateMesh(meshsize.moderate)).Curve(curveOrder) 
            mesh = Mesh(geo.GenerateMesh(msize)).Curve(curveOrder) 
        self.mesh=mesh
        
        #self.sigmaDomain="to_side|from_side"
        self.sigma_d={"conductor":sig,  "air":0, "reduced_region":0, "Kelvin":0,'default':0}
        self.rho_d={"conductor":1./sig,  "air":0, "reduced_region":0, "Kelvin":0, 'default':0}
        self.mu_d={"conductor":mu,  "air":mu, "reduced_region":mu, "Kelvin":mu, 'default':mu}
        self.Sigma = CoefficientFunction([self.sigma_d[mat] for mat in mesh.GetMaterials()]) 
        self.Mu = CoefficientFunction([self.mu_d[mat] for mat in mesh.GetMaterials()])  
        self.Rho = CoefficientFunction([self.rho_d[mat] for mat in mesh.GetMaterials()]) 
        #self.mu=4.e-7*math.pi
        self.conductor_boundary="conductor_boundary"
        self.conductive_region="conductor"
        #MeasureFace.SetMesh(mesh)
        self.Print()

        self.symmetric_plane="None"
        self.reduced_region="reduced_region"
        self.air_region="air|reduced_region" 
        self.total_region="conductor|air|"
        self.total_boundary="total_boundary"
        self.reduced_boundary="reduced_boundary"
        self.total_air_region="air"
        self.Bn0_boundary=""
        self.Ht0_boundary=""

        coil=Prob3Coil()
        self.coil=coil
        self.model=self.model-coil.geo
        self.model=Glue([self.model,coil.geo])

        if self.coil.geo :
            self.geo=model-self.coil.geo
            self.geo=Glue([self.geo,self.coil.geo])

        plotGeo=conductor*Box((-1,0,-1), (1,1,1))
        
        self.plotMesh = Mesh(OCCGeometry(plotGeo).GenerateMesh(maxh=1.e-3))

    def ReducedField(self, Bv, zero):
        Bs_dic = { "reduced_region":Bv, "total": zero, 'default':zero,
                  "conductor":zero, "air":zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in self.mesh.GetMaterials()])
        return Bs
 
    def GetGeometry(self):
        return self.model
    
    def GetMesh(self):
        return self.mesh

    def PlotBFieldonLine(self, mesh, BField):
        import matplotlib.pylab as plt
        x0=0.055
        y0=0
        z0=6.35e-3/2+0.5e-3
        dx=2*x0/100
        x=x0
        y=y0
        z=z0
        xp=[]
        ypreal=[]
        ypimag=[]
        for n in range(101):
            pnt=mesh(x,y,z)
            xp.append(x)
            ypreal.append(BField(pnt)[2].real)
            ypimag.append(BField(pnt)[2].imag)
            x=x-dx

        plt.plot(xp, ypreal )  
        plt.xlabel("x")  # Add x-axis label
        plt.plot(xp, ypimag ) 
        plt.xlabel("x")  # Add x-axis label
        plt.ylabel("Bz")  # Add y-axis label
        plt.show()  

    def PlotJ(self, JField, phase=0):
        import numpy as np
        import matplotlib.pyplot as plt
        eps=0.01e-3
        z=6.35e-3/2-eps
        nxgrid=55
        nygrid=30
        x_grid=np.linspace(-55e-3,55e-3, nxgrid)
        y_grid=np.linspace(-30e-3,30e-3, nygrid)
        X, Y= np.meshgrid(x_grid, y_grid)
        Jx, Jy, Jz=self.GetJ(X, Y, z, JField, phase)

        scale=10e-3
        plt.figure(figsize=(110e-3/scale, 60e-3/scale))
        #plt.quiver(X, Y, Jx, Jy, np.sqrt(Jx**2 + Jy**2 + Jz**2), cmap='viridis')
        #plt.pcolormesh(X, Y, np.sqrt(Jx**2 + Jy**2 + Jz**2), cmap='jet')
        #plt.quiver(X, Y, Jx, Jy, np.sqrt(Jx**2 + Jy**2 + Jz**2), cmap='jet')
        # 流線（Streamlines）の描画
        plt.streamplot(
            X, Y, Jx, Jy,
            color=np.sqrt(Jx**2 + Jy**2 + Jz**2), # マグニチュードで流線の色を変化させる
            cmap='jet',    # カラーマップを指定
            linewidth=1.5,     # 流線の幅
            density=1.5,       # 流線の密度（値が大きいほど線が増えます）
            arrowsize=1.5      # 流線上の方向を示す矢印のサイズ
        )
        plt.colorbar(label='Curent density [A/m^2]')
        plt.title(f'Magnetic Field, Phase={phase} deg')
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.axis('equal')
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    def GetJ(self, X, Y, z, JField, phase):
        import numpy as np
        mesh=self.mesh
        ex=exp(1j*phase*math.pi/180.)
        N=X.shape[0]
        M=X.shape[1]
        Jx=np.zeros((N, M))
        Jy=np.zeros((N, M))
        Jz=np.zeros((N, M))
        for n in range(N):
            for m in range(M):
                x=X[n][m]
                y=Y[n][m]
                mip = mesh(x,y,z)
                vx=JField(mip)[0]*ex
                vy=JField(mip)[1]*ex
                vz=JField(mip)[2]*ex
                Jx[n][m]=vx.real
                Jy[n][m]=vy.real
                Jz[n][m]=vz.real
        return Jx, Jy, Jz