import math
from ngsolve import *
from ngsolve.webgui import Draw
import numpy as np
import sys
sys.path.append('..\include')
from MatrixSolver import MatrixSolver as solver 
sys.path.append(r'..\bin\Release') 
from EMPY_Field import *
from Static_Method import Static_Method

class A_ReducedOmega_Method(Static_Method):
    def __init__(self,  model, coil,  **kwargs):  
        super().__init__(model, coil,  **kwargs)
        #self.Calc(model, coil,  **kwargs)
        
    def Calc(self, model, coil,  **kwargs):
        default_values = {"feOrder":1,
                          "boundaryCD":"Bn0",
                          "Kelvin": "off"
                         }
        default_values.update(kwargs)
        self.feOrder=default_values["feOrder"]
        boundaryCD=default_values["boundaryCD"]
        Kelvin=default_values["Kelvin"]
        self.Kelvin=Kelvin

        feOrder=self.feOrder
        self.mesh=model.mesh
        mesh=self.mesh
      
        import time
        start_time = time.perf_counter()

        total_region=model.total_region
        reduced_region=model.reduced_region
        total_boundary=model.total_boundary
        reduced_boundary=model.reduced_boundary
        Bn0_boundary=model.Bn0_boundary
        Ht0_boundary=model.Ht0_boundary
        self.Mu=model.Mu
        Mu=self.Mu
        self.total_region=total_region
        self.reduced_region=reduced_region       
        
        #coil=UNIF(0,0,1,0)
        Ov=Ofield(coil)
        Av=Afield(coil)
        Bv=Bfield(coil)
        mu0=4.e-7*math.pi
        Hs=Bv/mu0
        zero=(0,0,0)

        Bs_dic = {"iron":zero, "A_domain":zero, "Omega_domain":Bv, "Kelvin":zero, 'default':zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in mesh.GetMaterials()])

        fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, definedon=total_region, dirichlet=Bn0_boundary) 

        if Kelvin=="off":
            if boundaryCD=="Ht0":
                fesOmega=H1(mesh, order=feOrder, definedon=reduced_region, dirichlet=Ht0_boundary+"|"+reduced_boundary)
            elif boundaryCD=="Bn0":
                fesOmega=H1(mesh, order=feOrder, definedon=reduced_region, dirichlet=Ht0_boundary)
        else:
            fesOmega=H1(mesh, order=feOrder, definedon=reduced_region+"|Kelvin", dirichlet=Ht0_boundary)
            fesOmega=Periodic(fesOmega)

        fes=fesA*fesOmega   
        self.fes=fes
        (A, Omega),(N, o) = fes.TnT()
        normal = specialcf.normal(mesh.dim)
        a= BilinearForm(fes)
        a +=1/Mu*curl(A)*curl(N)*dx(total_region)
        a += -Mu*grad(Omega)*grad(o)*dx(reduced_region)
        if Kelvin=="on":
            rs=self.model.rKelvin
            xs=2*rs
            r=sqrt((x-xs)*(x-xs)+y*y+z*z)
            fac=rs*rs/r/r
            a += -Mu*fac*(grad(Omega)*grad(o))*dx("Kelvin") 

        a += -( Omega*(curl(N).Trace()*normal) + o*(curl(A).Trace()*normal) ) *ds(total_boundary)
        #a += ( Cross( grad(Omega).Trace(), N.Trace() )*normal + Cross(grad(o).Trace(), A.Trace() )*normal) *ds(total_boundary)

        f = LinearForm(fes)
        f +=-o*(Bv*normal)  *ds(total_boundary)
        f += Cross(N.Trace(),Hs)*normal*ds(total_boundary)
        with TaskManager():
            f.Assemble()

        gf=self.Solve(fes, a, f)

        gfA, gfOmega=gf.components

        fesOt=HCurl(mesh, order=feOrder, type1=True, nograds=True, definedon=total_region) 
        fesOr=H1(mesh, order=feOrder, definedon=reduced_region+"|Kelvin")
        Ot=GridFunction(fesOt)
        Or=GridFunction(fesOr)

        Ot.Set(gfA,VOL, definedon=total_region)
        Or.Set(gfOmega,VOL, definedon=reduced_region+"|Kelvin")
        print("*** Omega ***  ")
        Draw(Or, mesh)
        
        self.Bt=curl(gfA)
        self.Br=mu0*grad(gfOmega)
        BField=self.Bt+self.Br+Bs

        print("feOrder=", feOrder,"  ", "ndof=",fes.ndof,"  ")
        self.CalcResult(model, BField)

        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f"経過時間: {elapsed_time:.4f} 秒  ")

    def CalcError(self):
        mesh=self.mesh
        feOrder=self.feOrder
        total_region=self.total_region
        reduced_region=self.reduced_region
        Mu=self.Mu

        fesHd = HCurl(mesh, order=feOrder-1, type1=True, nograds=True, definedon=total_region)
        Hd=GridFunction(fesHd, "flux")
        H = 1/Mu*self.Bt
        Hd.Set(H)
        err=Mu*(H-Hd)*(H-Hd)
        eta = Integrate(err, mesh, VOL, element_wise=True)

        #gfur=self.Or
        fsBdr = HDiv(mesh, order=feOrder-1, definedon=reduced_region)
        Bdr=GridFunction(fsBdr, "flux")
        Br = self.Br
        Bdr.Set(Br)
        err_r=1/Mu*(Br-Bdr)*(Br-Bdr)
        eta_r = Integrate(err_r, mesh, VOL, element_wise=True) 

        if self.Kelvin=="on":
            fkBdr = HDiv(mesh, order=feOrder-1, definedon="Kelvin")
            Bdr=GridFunction(fkBdr, "flux")
            
            rs=self.model.rKelvin
            xs=2*rs
            r=sqrt((x-xs)*(x-xs)+y*y+z*z)
            fac=rs*rs/r/r
            Br=Br*fac
            Bdr.Set(Br)
            err_rk=1/(Mu*fac)*(Br-Bdr)*(Br-Bdr)
            eta_rk =Integrate(err_rk, mesh, VOL, element_wise=True) 
            print(" maxerr = ", max(eta), max(eta_r), max(eta_rk) )
            eta_r = eta_r+eta_rk

        error=eta+eta_r
        maxerr = max(error)
        print ("ndof=", self.fes.ndof, " maxerr = ", maxerr)
        return maxerr, error        