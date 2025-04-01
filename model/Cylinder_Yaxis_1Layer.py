from netgen.meshing import *
from netgen.csg import *
from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw
from netgen.webgui import Draw as DrawGeo
import math


class Cylinder_Yaxis_1Layer():
    def __init__(self,  **kwargs):
        default_values = {"name":"Cylinder_Yaxis_1Layer",
                          "curveOrder":3,
                          "sigma": 0,
                          "ndiv":5,
                          "muiron": 10
                         }
        default_values.update(kwargs)
        sigma=default_values["sigma"]
        ndiv=default_values["ndiv"]
        curveOrder=default_values["curveOrder"]
        muiron=default_values["muiron"]

        iron = Cylinder((0,0,-0.05), Z, r=0.1, h=0.1)
        iron.faces.name="iron_boundary"
        iron.mat("iron")

        total = Cylinder((0,0,-0.05), Z, r=0.2, h=0.1)
        total.faces.name = "total_boundary"

        reduced=Cylinder((0,0,-0.05), Z, r=1.0, h=0.1)
        reduced.faces.name = "reduced_boundary"

        reduced = reduced-total
        total = total-iron
        reduced.mat("reduced")
        total.mat("total")

        iron.faces.Min(Z).name="iron_bottom"
        iron.faces.Max(Z).name="iron_top"
        iron.faces.Min(Z).Identify(iron.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)
        iron.maxh=0.1/ndiv


        total.faces.Min(Z).name="total_bottom"
        total.faces.Max(Z).name="total_top"
        total.faces.Min(Z).Identify(total.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)
        total.maxh=0.2/ndiv

        reduced.faces.Min(Z).name="reduced_bottom"
        reduced.faces.Max(Z).name="reduced_top"
        reduced.faces.Min(Z).Identify(reduced.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)
        reduced.maxh=1.0/ndiv

        geo =OCCGeometry(Glue([iron, total, reduced]))
        ngmesh = geo.GenerateMesh( quad_dominated=False)
        ngmesh.ZRefine("bot-top", [])
        mesh = Mesh(ngmesh).Curve(curveOrder)
        mu=4.e-7*math.pi
        mu_d={"iron":mu*muiron, "total":mu, "reduced":mu,'default':mu}
        self.Mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()]) 
        if sigma !=0:
            sigma_d={"iron":sigma, "total":0, "reduced":0,'default':0}
            self.Sigma = CoefficientFunction([sigma_d[mat] for mat in mesh.GetMaterials()])
        self.conductive_region="iron"
        self.symmetric_plane="iron_top|iron_bottom|total_top|total_bottom|reduced_top|reduced_bottom"
        self.total_region="iron|total"
        self.reduced_region="reduced"
        self.air_region="total|reduced"
        self.total_boundary="total_boundary"
        self.reduced_boundary="reduced_boundary"
        self.conductor_boundary="iron_boundary"
        self.mesh=mesh
        
        #Draw(mesh)

    def ReducedField(self, Bv, zero):
        Bs_dic = {"iron":zero, "total": zero, "reduced":Bv, 'default':zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in self.mesh.GetMaterials()])
        return Bs
