from netgen.meshing import *
from netgen.csg import *
from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw

#import math
class CubeMesh():

    def __init__(self,  **kwargs):  
        self.SetMesh(**kwargs)

    def SetMesh(self, **kwargs):
        default_values = {"name":"Cube",
                          "mur":1000, 
                          "msize": meshsize.moderate,
                          "ndiv":5,
                          "type":0
                         }
        
        default_values.update(kwargs)
        self.name=default_values["name"]
        self.mur=default_values["mur"]
        self.msize=default_values["msize"]
        ndiv=default_values["ndiv"]
        type=default_values["type"]
        
        curveOrder=1
        iron = Box((0,0,0),(1, 1, 1))
        iron.faces.Min(X).name="Bn0"
        iron.faces.Min(Y).name="Bn0"
        iron.faces.Min(Z).name="Ht0"
        if type==1:
            iron.faces.Max(X).name="A_Omega_boundary"
            iron.faces.Max(Y).name="A_Omega_boundary"
            iron.faces.Max(Z).name="A_Omega_boundary"
        iron.mat("iron")
        #iron.maxh=1.0/ndiv
 
        A_domain = Box((0,0,0),(1.5, 1.5, 1.5))
        A_domain=A_domain-iron
        #A_domain.faces.Min(Z).Identify(A_domain.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)

        A_domain.faces.Min(X).name="Bn0"
        A_domain.faces.Min(Y).name="Bn0"
        A_domain.faces.Min(Z).name="Ht0"
        if type==0:
            A_domain.faces.Max(X).name="A_Omega_boundary"
            A_domain.faces.Max(Y).name="A_Omega_boundary"
            A_domain.faces.Max(Z).name="A_Omega_boundary"
            A_domain.mat("A_domain")
        elif type==1:
            A_domain.mat("Omega_domain")
        #A_domain.maxh=1.0/ndiv

        Omega_domain = Box((0,0,0),(10,10,10))
        Omega_domain=Omega_domain-A_domain-iron
        Omega_domain.faces.Min(X).name="Bn0"
        Omega_domain.faces.Min(Y).name="Bn0"
        Omega_domain.faces.Min(Z).name="Ht0"
        Omega_domain.faces.Max(Y).name="Omega0"
        Omega_domain.faces.Max(X).name="Omega0"
        Omega_domain.faces.Max(Z).name="Omega0"
        #Omega_domain.faces.Min(Z).Identify(Omega_domain.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)
        Omega_domain.mat("Omega_domain")
        #Omega_domain.maxh=5.0/ndiv


        geo=Glue([iron, A_domain, Omega_domain])
        occgeo =OCCGeometry(geo)
        #ngmesh = occgeo.GenerateMesh(self.msize, quad_dominated=False)
        #ngmesh = occgeo.GenerateMesh(grading=0.05)
        ngmesh = occgeo.GenerateMesh(self.msize)
        #ngmesh.ZRefine("bot-top", [])
        mesh = Mesh(ngmesh)
        #Draw(mesh)
        #for n in range(0, ndiv):
        #    mesh.Refine()

        self.geo=geo
        self.mesh=mesh
        if type==0:
            self.reduced_region="Omega_domain"
            self.total_region="iron|A_domain"
        elif type==1:
            self.reduced_region="Omega_domain"
            self.total_region="iron"           
        self.total_boundary="A_Omega_boundary"
        self.reduced_boundary="Omega0"
        self.Bn0_boundary="Bn0"
        self.Ht0_boundary="Ht0"

        import math
        mur=self.mur
        mu0=4.e-7*math.pi
        mu=mu0*mur
        mu_d={"iron":mu,  "A_domain":mu0, "Omega_domain":mu0, 'default':mu0}
        self.Mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        
        print("nv=", mesh.nv, " nedge=", mesh.nedge, " nfacet=", mesh.nfacet, " ne=",mesh.ne)

    def Print(self):
        geo=self.geo
        mesh=self.mesh
        print("Model:", self.name, "mur=", self.mur)
        print("nv=", mesh.nv, " nedge=", mesh.nedge, " nfacet=", mesh.nfacet, " ne=",mesh.ne)
        for s in geo.solids:
            print("name:",s.name, "  mass:", s.mass, "  center:", s.center)

class SphereMesh():
    def __init__(self,  **kwargs):  
        self.SetMesh(**kwargs)

    def SetMesh(self, **kwargs):
        default_values = {"name":"Spere",
                          "mur":1000, 
                          "curveOrder":3, 
                          "msize": meshsize.moderate,
                          "ndiv":5,
                          "type":0
                         }
        
        default_values.update(kwargs)
        self.name=default_values["name"]
        self.mur=default_values["mur"]
        self.curveOrder=default_values["curveOrder"]
        self.msize=default_values["msize"]
        ndiv=default_values["ndiv"]
        type=default_values["type"]
        
        curveOrder=3
        iron = Sphere( (0,0,0), 1)
        if type==1:
            iron.faces[0].name="A_Omega_boundary"
        iron.mat("iron")
        iron.maxh=0.5/ndiv
 
        A_domain =  Sphere( (0,0,0), 1.2)
        if type==0:
            A_domain.faces[0].name="A_Omega_boundary"
        A_domain=A_domain-iron
        #A_domain.faces.Min(Z).Identify(A_domain.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)

        if type==0:
            A_domain.mat("A_domain")
        elif type==1:
            A_domain.mat("Omega_domain")
        A_domain.maxh=0.6/ndiv

        Omega_domain = Sphere( (0,0,0), 5.)
        Omega_domain.faces[0].name="Omega0"
        Omega_domain=Omega_domain-A_domain-iron
        #Omega_domain.faces.Min(Z).Identify(Omega_domain.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)
        Omega_domain.mat("Omega_domain")
        Omega_domain.maxh=2.5/ndiv

        geo=Glue([iron, A_domain, Omega_domain])
        occgeo =OCCGeometry(geo)
        
        ngmesh = occgeo.GenerateMesh(grading=0.2)
        #ngmesh = occgeo.GenerateMesh(self.msize, quad_dominated=False)
        #ngmesh = occgeo.GenerateMesh(quad_dominated=False).Curve(curveOrder)
        #ngmesh.ZRefine("bot-top", [])
        mesh = Mesh(ngmesh).Curve(self.curveOrder)
        #Draw(mesh)
        #for n in range(0, ndiv):
        #    mesh.Refine()

        self.geo=geo
        self.mesh=mesh
        if type==0:
            self.reduced_region="Omega_domain"
            self.total_region="iron|A_domain"
        elif type==1:
            self.reduced_region="Omega_domain"
            self.total_region="iron"           
        self.total_boundary="A_Omega_boundary"
        self.reduced_boundary="Omega0"
        self.Bn0_boundary=""
        self.Ht0_boundary=""

        import math
        mur=self.mur
        mu0=4.e-7*math.pi
        mu=mu0*mur
        mu_d={"iron":mu,  "A_domain":mu0, "Omega_domain":mu0, 'default':mu0}
        self.Mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        
        print("nv=", mesh.nv, " nedge=", mesh.nedge, " nfacet=", mesh.nfacet, " ne=",mesh.ne)

    def Print(self):
        geo=self.geo
        mesh=self.mesh
        print("Model:", self.name, "mur=", self.mur, "curveOrder=", self.curveOrder)
        print("nv=", mesh.nv, " nedge=", mesh.nedge, " nfacet=", mesh.nfacet, " ne=",mesh.ne)
        for s in geo.solids:
            print("name:",s.name, "  mass:", s.mass, "  center:", s.center)

