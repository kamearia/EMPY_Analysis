from netgen.meshing import *
from netgen.csg import *
from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw
import sys
sys.path.append(r'..\COIL\include')
from EMPY_COIL import *
sys.path.append(r'..\bin\Release') 
from EMPY_Field import *

from netgen.meshing import MeshingParameters
#from netgen.meshing import ngsglobals # 追記
#from netgen.libngpy import ngsglobals
#from netgen import ngsglobals 

#import math
class SphereMesh():
    def __init__(self,  **kwargs):  
        self.SetMesh(**kwargs)

    def SetMesh(self, **kwargs):
        default_values = {"name":"Cube",
                          "mur":1000, 
                          "sigma": 5.e7,
                          "msize": meshsize.moderate,
                          "ndiv":5,
                          "type":0,
                          "curveOrder":1,
                          "rOuter": 5,
                          "rKelvin":5
                         }
        
        default_values.update(kwargs)
        self.name=default_values["name"]
        self.mur=default_values["mur"]
        self.sigma=default_values["sigma"]
        self.msize=default_values["msize"]
        ndiv=default_values["ndiv"]
        type=default_values["type"]
        curveOrder=default_values["curveOrder"]   
        rKelvin=default_values["rKelvin"] 
        rOuter=default_values["rOuter"]  
        self.curveOrder=curveOrder
        self.rKelvin=rKelvin
        if rKelvin !=0: rOuter=rKelvin

        self.radius=1.
        r= self.radius
        iron = Sphere(Pnt(0,0,0.0), r=r)*Box((0,0,0), (r,r,r))
        iron.faces.Min(X).name="Bn0"
        iron.faces.Min(Y).name="Bn0"
        iron.faces.Min(Z).name="Ht0"
        if type==1:
            iron.face[0].name="A_Omega_boundary"
        iron.mat("iron")
        #iron.maxh=1.0/ndiv

        A_domain =  Sphere(Pnt(0,0,0.0), r=1.5)*Box((0,0,0), (1.5,1.5,1.5))
        A_domain.faces.Min(X).name="Bn0"
        A_domain.faces.Min(Y).name="Bn0"
        A_domain.faces.Min(Z).name="Ht0"
        """
        A_domain =  Box((0,0,0), (1.5,1.5,1.5))
        A_domain.faces.Min(X).name="Bn0"
        A_domain.faces.Min(Y).name="Bn0"
        A_domain.faces.Min(Z).name="Ht0"
        """
        if type==0:
            A_domain.faces[0].name="A_Omega_boundary"
            """
            A_domain.faces.Max(X).name="A_Omega_boundary"
            A_domain.faces.Max(Y).name="A_Omega_boundary"
            A_domain.faces.Max(Z).name="A_Omega_boundary"
            """
            A_domain.mat("A_domain")
        elif type==1:
            A_domain.mat("Omega_domain")
        #A_domain.maxh=1.0/ndiv
        

        Omega_domain=Sphere(Pnt(0,0,0.0), r=rOuter)*Box((0,0,0), (rOuter,rOuter,rOuter))
        #Omega_domain.faces[0].Identify(Omega_domain.faces[4], "ud0",  IdentificationType.PERIODIC)
        Omega_domain.faces.Min(X).name="Bn0"
        Omega_domain.faces.Min(Y).name="Bn0"
        Omega_domain.faces.Min(Z).name="Ht0"
        Omega_domain.faces[0].name="Omega0"
        #Omega_domain.faces.Max(Y).name="Omega0"
        #Omega_domain.faces.Max(X).name="Omega0"
        #Omega_domain.faces.Max(Z).name="Omega0"
        #Omega_domain.faces.Min(Z).Identify(Omega_domain.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)
        Omega_domain.mat("Omega_domain")
        #Omega_domain.maxh=rk/5

        if rKelvin !=0:
            rk=rKelvin
            kcenter=rk*2
            self.kcenter=kcenter
            external_domain = Sphere(Pnt(kcenter,0,0.), r=rk)*Box((kcenter,0,0), (kcenter+rk,rk,rk))
            external_domain.faces.Min(X).name="Bn0"
            external_domain.faces.Min(Y).name="Bn0"
            external_domain.faces.Min(Z).name="Ht0"
            external_domain.mat("Kelvin")
            #external_domain.maxh=rk/5
            external_domain.faces[0].Identify(Omega_domain.faces[0], "ud0",  IdentificationType.PERIODIC)
            geo=Glue([iron, A_domain, Omega_domain, external_domain])
        else:
             geo=Glue([iron, A_domain, Omega_domain]) 
            
        occgeo =OCCGeometry(geo)
        
        # 最小角度を0.1ラジアン（約5.7度）に設定
        #ngsglobals.min_angle = 0.1 
        
        # OCCGeometryから MeshingParameters オブジェクトを取得し、設定を適用
        #mp = MeshingParameters(grading=0.2, optimize3d="m")
        #mp = MeshingParameters()
        # 最小角度を0.15ラジアン（約8.6度）に設定
        # 3Dテトラヘドラの品質を保証するための Netgen パラメータ
        #mp.minimal_angle = 0.15 
        
        # メッシュ生成後に最適化を有効にする
        #mp.optimize = True 
        
        # 3D最適化を明示的に有効にする（optimize=Trueに包括されていることが多いが、念のため）
        #mp.optimize3d = True
        """
        # 🚀 修正点 2: maxh (メッシュサイズ) も MeshingParameters に設定する
        # self.msize が meshsize.moderate のような Enum の場合、.value で数値を取得
        if isinstance(self.msize, meshsize):
            maxh = self.msize.value
        else:
            maxh = self.msize
        """
        #mp = MeshingParameters(maxh=5, grading=0.2, optimize3d="m")
        #mp = MeshingParameters(grading=0.2, optimize3d="m")
      
        """
        #ngmesh = occgeo.GenerateMesh(self.msize, quad_dominated=False)
        #ngmesh = occgeo.GenerateMesh(grading=0.05)
        MIN_ANGLE = 0.15
        """
        #ngmesh = occgeo.GenerateMesh(self.msize, meshing_parameters=mp)#, minimal_angle=MIN_ANGLE, optimize3d=True,)
        ngmesh = occgeo.GenerateMesh(self.msize)
        #ngmesh = occgeo.GenerateMesh(mp=mp)
        #ngmesh.ZRefine("bot-top", [])
        print("curveOrder=", curveOrder)
        mesh = Mesh(ngmesh).Curve(curveOrder)
        #mesh.Optimize()

        self.geo=geo
        self.mesh=mesh
        if type==0:
            self.reduced_region="Omega_domain"
            self.total_region="iron|A_domain"
        elif type==1:
            self.reduced_region="A_domain|Omega_domain"
            self.total_region="iron"
        self.conductive_region="iron"
        self.total_boundary="A_Omega_boundary"
        self.reduced_boundary="Omega0"
        self.Bn0_boundary="Bn0"
        self.Ht0_boundary="Ht0"

        self.coil=EMPY_UNIF(0,0,1,0)
        
        import math
        mur=self.mur
        mu0=4.e-7*math.pi
        mu=mu0*mur
        mu_d={"iron":mu,  "A_domain":mu0, "Omega_domain":mu0,"Kelvin":mu0, 'default':mu0}
        self.Mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        sigma_d={"iron":self.sigma,  "A_domain":0, "Omega_domain":0,"Kelvin":0, 'default':0}
        self.Sigma = CoefficientFunction([sigma_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        rho_d={"iron":1/self.sigma,  "A_domain":0, "Omega_domain":0,"Kelvin":0, 'default':0}
        self.Rho = CoefficientFunction([rho_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
               
        print("nv=", mesh.nv, " nedge=", mesh.nedge, " nfacet=", mesh.nfacet, " ne=",mesh.ne)

    def Print(self):
        geo=self.geo
        mesh=self.mesh
        print("Model:", self.name, "mur=", self.mur)
        print("nv=", mesh.nv, " nedge=", mesh.nedge, " nfacet=", mesh.nfacet, " ne=",mesh.ne)
        for s in geo.solids:
            print("name:",s.name, "  mass:", s.mass, "  center:", s.center)