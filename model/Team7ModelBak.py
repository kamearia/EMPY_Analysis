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
        #self.coil=Prob7Coil()
        
        self.coil=Prob7Coil()
        if self.coil.geo :
            self.coil.geo.faces.col = (0,1,0,1)
            #self.geo=geo-self.coil.geo
            self.geo=Glue([self.coil.geo, self.geo])

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

    def PlotBFieldonLine(self, mesh, BField):
        import matplotlib.pylab as plt
        x0=0.
        y0=72.e-3
        z0=34.e-3
        N=100
        dx=288e-3/N
        #dx=18.e-3
        x=x0
        y=y0
        z=z0
        xp=[]
        ypreal=[]
        ypimag=[]
        #for n in range(17):
        for n in range(N+1):
            pnt=mesh(x,y,z)
            xp.append(x)
            ypreal.append(BField(pnt)[2].real)
            ypimag.append(-BField(pnt)[2].imag)
            x=x+dx

        plt.plot(xp, ypreal )  
        plt.xlabel("x")  # Add x-axis label
        plt.plot(xp, ypimag ) 
        plt.xlabel("x")  # Add x-axis label
        plt.ylabel("Bz")  # Add y-axis label
        plt.show() 

