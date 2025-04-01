# -*- coding: utf-8 -*-
from netgen.occ import *
from ngsolve import TaskManager
import math
from ngsolve import *
import sys
sys.path.append('include')
from EMPY_COIL import EMPY_LOOP

def Prob3Coil():
    rin=20e-3
    rout=40e-3
    h=20e-3
    thick=6.35e-3
    z=thick/2+15e-3
    x=20e-3
    radius=(rin+rout)/2
    radialWidth=rout-rin
    axialWidth=h
    zp=z+h/2
    center=gp_Pnt(20e-3,0, zp)
    J=1260  #(AT)
    
    loop=EMPY_LOOP(J, radius, center, radialWidth, axialWidth)
    return loop
    
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
        print("conductiveDomain:", self.conductiveDomain)
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
                          "outer":"box",
                          "wz":6.35e-3,
                          "boxx":255e-3, 
                          "boxy":280e-3,      
                          "boxz":400e-3, 
                          "div_thick":1,
                          "sigma": 0.3278e8 
                         }
        default_values.update(kwargs)
        holes=default_values["holes"]
        outerBox=default_values["outer"]
        bx=default_values["boxx"]
        by=default_values["boxy"]
        bz=default_values["boxz"]
        wz=default_values["wz"]
        div_thick=default_values["div_thick"]
        print(default_values["name"], "   holes:", holes, "   Thickness:", wz)
        print("boxx= ", bx," boxy= ", by, " boxz= ",bz)
        print("div_thick= ", div_thick)
        sig=default_values["sigma"]
        curveOrder=2
        mu = 4*math.pi*1e-7
        self.up_down=2.*bz

        wx=110e-3
        wy=60e-3
        wz=default_values["wz"]
        hx=30e-3
        hy=40e-3
        sx=20e-3
        conductor=Box((-wx/2,-wy/2,-wz/2),(wx/2,wy/2,wz/2))
        conductor.maxh=wz/div_thick
        conductor.faces.name="conductorBND"
        if holes==0: conductor.mat("conductor")
        d=5.e-3
        total=Box((-wx/2-d,-wy/2-d,-wz/2-d),(wx/2+d,wy/2+d,wz/2+d))
        total.faces.name="totalBoundary"
        outer=Sphere(Pnt(0,0,0),bx)
        if outerBox=="box": outer=Box((-bx,-by,-bz),(bx,by,bz))
        outer.faces.name="outer"
        reduced=outer-total
        reduced.mat("reduced")
        air=total-conductor
        air.mat("air")
        #reduced.faces.Max(Z).name="outer"
        #reduced.faces.Min(Z).name="outer"
        #reduced.faces.Max(X).name="outer"
        #reduced.faces.Min(X).name="outer"
        #reduced.faces.Max(Y).name="outer"
        #reduced.faces.Min(Y).name="outer"
 
        plane =Box((-bx, 0, -bz), (bx, by, bz))
    
        if holes==0:
            model=Glue([conductor,air, reduced])
        else:
            hole1=Box((-hx/2+sx,-hy/2,-wz/2),(hx/2+sx,hy/2,wz/2))
            hole1.mat("hole")
            hole1.faces.Min(Z).name="hole_lower"
            hole1.faces.Max(Z).name="hole_upper"
            hole1.faces.Min(X).name="interface"
            hole1.faces.Max(X).name="interface"
            hole1.faces.Min(Y).name="interface"
            hole1.faces.Max(Y).name="interface"
            HolePot("hole", "hole_lower", "hole_upper","interface")
            #hole1.maxh=wz/div_thick
            conductor= conductor-hole1
            if holes==0: conductor.mat("conductor")
            #conductor.maxh=wz/div_thick
            con2=conductor-plane
            con1=conductor-con2
            #con1.maxh=wz/div_thick
            #con2.maxh=wz/div_thick
            con2.mat("from_side")
            con1.mat("to_side")

            MeasureFace("cut1", "from_side", "to_side")
            if holes==1:
                #conductor=Glue([con1,con2])
                model=Glue([hole1,con1,con2, air, reduced])
                model.faces[14].name="cut1"
            elif holes==2:
                hole2=Box((-hx/2-sx,-hy/2,-wz/2),(hx/2-sx,hy/2,wz/2))
                hole2.mat("hole2")
                hole2.faces.Min(Z).name="hole_lower2"
                hole2.faces.Max(Z).name="hole_upper2"
                hole2.faces.Min(X).name="interface2"
                hole2.faces.Max(X).name="interface2"
                hole2.faces.Min(Y).name="interface2"
                hole2.faces.Max(Y).name="interface2"
                #hole2.maxh=wz/div_thick
                con1=con1-hole2
                con2=con2-hole2
                model=Glue([hole1,hole2, con1, con2, air, reduced])

                HolePot("hole2", "hole_lower2", "hole_upper2","interface2")
                #hole2.maxh=wz/div_thick
                #conductor= conductor-hole2
                #conductor.mat("sig")
                #conductor.maxh=wz/div_thick
                #model=Glue([hole1, hole2, conductor,air])
                model.faces[22].name="cut1"
                model.faces[20].name="cut2"
                model.faces[26].name="cut3"
                MeasureFace("cut2", "from_side", "to_side")
                MeasureFace("cut3", "from_side", "to_side")
                
        if holes==0:
            self.conductiveDomain="conductor"
            self.totalDomain="conductor|air"
        else: 
            self.conductiveDomain="to_side|from_side"
            self.totalDomain="to_side|from_side|air"
        self.airDomain="hole1|hole2|air|reduced"   
            
        self.model=model
        #conductor.faces.name="conductorBND"

        geo =OCCGeometry(model)
        self.geo=geo
        with TaskManager():
            #mesh = Mesh(geo.GenerateMesh(meshsize.coarse, quad_dominated=True)).Curve(1)
            #mesh = Mesh(geo.GenerateMesh(meshsize.fine)).Curve(1)
            mesh = Mesh(geo.GenerateMesh(meshsize.coarse)).Curve(curveOrder) 
        self.mesh=mesh
        
        #self.sigmaDomain="to_side|from_side"
        self.sigma_d={"conductor":sig, "to_side":sig, "from_side":sig,  "air":0, "reduced":0,  "hole":0, "hole2":0,'default':0}
        self.mu_d={"conductor":mu, "to_side":mu, "from_side":mu,  "air":mu, "reduced":mu, "hole":mu, "hole2":mu,'default':mu}
        self.sigma = CoefficientFunction([self.sigma_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        self.Mu = CoefficientFunction([self.mu_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        #self.mu=4.e-7*math.pi
        self.outer_boundary="outer"
        self.conductor_boundary="conductorBND"
        for n in range(HolePot.num):
            self.conductor_boundary=self.conductor_boundary+"|"+HolePot.holes[n].interface
        MeasureFace.SetMesh(mesh)
        self.Print()
        
        self.symmetric_plane="XXX"
        self.reduced_region="reduced"
        self.total_region="conductor|to_side|from_side|air|hole|hole2"
        self.total_boundary="totalBoundary"
        self.reduced_boundary="outer"

        coil=Prob3Coil()
        self.coil=coil
        self.model=self.model-coil.geo
        self.model=Glue([self.model,coil.geo])
 
    def GetGeometry(self):
        return self.model
    
    def GetMesh(self):
        return self.mesh

