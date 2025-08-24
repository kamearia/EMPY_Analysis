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
                          "msize": meshsize.moderate
                         }
        
        default_values.update(kwargs)
        self.name=default_values["name"]
        self.mur=default_values["mur"]
        self.msize=default_values["msize"]
        
        curveOrder=1
        iron = Box((0,0,0),(1, 1, 1))
        iron.faces.Min(X).name="Bn0"
        iron.faces.Min(Y).name="Bn0"
        iron.faces.Min(Z).name="Ht0"
        iron.mat("iron")
        #iron.maxh=0.1
 
        A_domain = Box((0,0,0),(1.2, 1.2,1.2))
        A_domain=A_domain-iron
        #A_domain.faces.Min(Z).Identify(A_domain.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)

        A_domain.faces.Min(X).name="Bn0"
        A_domain.faces.Min(Y).name="Bn0"
        A_domain.faces.Min(Z).name="Ht0"
        A_domain.faces.Max(X).name="A_Omega_boundary"
        A_domain.faces.Max(Y).name="A_Omega_boundary"
        A_domain.faces.Max(Z).name="A_Omega_boundary"
        A_domain.mat("A_domain")
        #A_domain.maxh=0.1

        Omega_domain = Box((0,0,0),(5,5,5))
        Omega_domain=Omega_domain-A_domain-iron
        Omega_domain.faces.Min(X).name="Bn0"
        Omega_domain.faces.Min(Y).name="Bn0"
        Omega_domain.faces.Min(Z).name="Ht0"
        Omega_domain.faces.Max(Y).name="Omega0"
        Omega_domain.faces.Max(X).name="Omega0"
        Omega_domain.faces.Max(Z).name="Omega0"
        #Omega_domain.faces.Min(Z).Identify(Omega_domain.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)
        Omega_domain.mat("Omega_domain")
        #Omega_domain.maxh=0.25


        geo=Glue([iron, A_domain, Omega_domain])
        occgeo =OCCGeometry(geo)
        ngmesh = occgeo.GenerateMesh(self.msize, quad_dominated=False)
        #ngmesh.ZRefine("bot-top", [])
        mesh = Mesh(ngmesh).Curve(curveOrder)
        #Draw(mesh)

        self.geo=geo
        self.mesh=mesh
        self.reduced_region="Omega_domain"
        self.total_region="iron|A_domain"
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

    def Print(self):
        geo=self.geo
        mesh=self.mesh
        print("Model:", self.name, "mur=", self.mur)
        print("nv=", mesh.nv, " nedge=", mesh.nedge, " nfacet=", mesh.nfacet, " ne=",mesh.ne)
        for s in geo.solids:
            print("name:",s.name, "  mass:", s.mass, "  center:", s.center)