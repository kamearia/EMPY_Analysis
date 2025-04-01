import math
from netgen.occ import *
from ngsolve import *
from Prob7Coil import Prob7Coil
class Team7Model():
    def __init__(self,  **kwargs):
        default_values = {"name":"Team Workshop Problem #7 Model",
                          "sigma": 3.526e7, 
                          "msize": meshsize.moderate,
                          "boxsize": 1.0
                         }
        default_values.update(kwargs)
        sigma=default_values["sigma"]
        msize=default_values["msize"]
        boxsize=default_values["boxsize"]
        xl=-1.353
        xu=1.637
        yl=-1.357
        yu=1.647
        zl=-0.3
        zu=0.449

        pc0=gp_Pnt(0,0,0)
        pc1=gp_Pnt(294.e-3,294.e-3, 19.e-3)
        conductor=Box(pc0,pc1)
        center=gp_Pnt(0.5*(pc0[0]+pc1[0]),0.5*(pc0[1]+pc1[1]),0.5*(pc0[2]+pc1[2]))
        ph0=(18e-3,18e-3,0)
        ph1=(126.e-3,126.e-3, 19.e-3)
        hole=Box(ph0,ph1)
        hole.mat("hole")
        conductor=conductor-hole
        conductor.mat("conductor")
        conductor.faces.name="conductor_boundary"
        
        #d=15e-3
        d=30e-3
        pt0=pc0-gp_Vec(d,d,d)
        pt1=pc1+gp_Vec(d,d,d)
        total=Box(pt0,pt1)
        total.faces.name="total_boundary"

        air=total-conductor-hole
        air.mat("air")

        if boxsize ==0:
            p0=gp_Pnt(xl,yl,zl)
            p1=gp_Pnt(xu,yu,zu)
        else:
            d=boxsize
            p0=center-gp_Vec(d,d,d)
            p1=center+gp_Vec(d,d,d)
        reduced=Box(p0,p1)
        reduced.faces.name="reduced_boundary"
        reduced=reduced-total
        reduced.mat("reduced_region")

        geo=Glue([conductor, hole, air, reduced])

        with TaskManager():
            #mesh = Mesh(geo.GenerateMesh(meshsize.coarse, quad_dominated=True)).Curve(1)
            #mesh = Mesh(geo.GenerateMesh(meshsize.fine)).Curve(1)
            #mesh = Mesh(geo.GenerateMesh(meshsize.moderate)).Curve(curveOrder) 
            mesh = Mesh(OCCGeometry(geo).GenerateMesh(msize)).Curve(1) 
        self.mesh=mesh

        self.sigma_d={"conductor":sigma, "hole":0,  "air":0, "reduced_region":0, 'default':0}
        mu = 4*math.pi*1e-7
        self.mu_d={"conductor":mu, "hole":mu,  "air":mu, "reduced_region":mu, 'default':mu}
        self.Sigma = CoefficientFunction([self.sigma_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        self.Mu = CoefficientFunction([self.mu_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性
        
        self.conductor_boundary="conductor_boundary"
        self.conductive_region="conductor"
        self.symmetric_plane="None"
        self.reduced_region="reduced_region"
        self.air_region="hole|air|reduced_region" 
        self.total_region="conductor|hole|air"
        self.total_boundary="total_boundary"
        self.reduced_boundary="reduced_boundary"
        self.total_air_region="air"

        self.coil=Prob7Coil()
        self.geo=geo-self.coil.geo
        self.geo=Glue([self.geo,self.coil.geo])

    def ReducedField(self, Bv, zero):
        Bs_dic = { "reduced_region":Bv,  "air": zero, "total": zero, 'default':zero,
                  "conductor":zero, "hole":zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in self.mesh.GetMaterials()])
        return Bs