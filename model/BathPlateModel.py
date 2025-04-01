# -*- coding: utf-8 -*-
from netgen.occ import *
from ngsolve import TaskManager
import math
from ngsolve import *
import sys
sys.path.append('..\COIL\include')
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
                          "msize": meshsize.moderate
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
        print(default_values["name"], "   holes:", holes, "   Thickness:", wz)
        print("boxx= ", bx," boxy= ", by, " boxz= ",bz)
        print("div_thick= ", div_thick)
        sig=default_values["sigma"]
        curveOrder=3
        mu = 4*math.pi*1e-7
        self.up_down=2.*bz

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
        reduced=outer-total
        reduced.mat("reduced_region")
        air=total-conductor
        air.mat("air")

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
        self.sigma_d={"conductor":sig,  "air":0, "reduced_region":0, 'default':0}
        self.mu_d={"conductor":mu,  "air":mu, "reduced_region":mu, 'default':mu}
        self.Sigma = CoefficientFunction([self.sigma_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        self.Mu = CoefficientFunction([self.mu_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
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

        coil=Prob3Coil()
        self.coil=coil
        self.model=self.model-coil.geo
        self.model=Glue([self.model,coil.geo])

    def ReducedField(self, Bv, zero):
        Bs_dic = { "reduced_region":Bv, "total": zero, 'default':zero,
                  "conductor":zero, "air":zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in self.mesh.GetMaterials()])
        return Bs
 
    def GetGeometry(self):
        return self.model
    
    def GetMesh(self):
        return self.mesh

