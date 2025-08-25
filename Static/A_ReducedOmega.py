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
                         }
        default_values.update(kwargs)
        self.feOrder=default_values["feOrder"]
        boundaryCD=default_values["boundaryCD"]

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
        
        #coil=UNIF(0,0,1,0)
        Av=Afield(coil)
        Bv=Bfield(coil)
        mu0=4.e-7*math.pi
        Hs=Bv/mu0
        zero=(0,0,0)

        Bs_dic = {"iron":zero, "A_domain":zero, reduced_region:Bv, 'default':zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in mesh.GetMaterials()])

        fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, definedon=total_region, dirichlet=Bn0_boundary) 

        if boundaryCD=="Ht0":
            fesOmega=H1(mesh, order=feOrder, definedon=reduced_region, dirichlet=Ht0_boundary+"|"+reduced_boundary)
        elif boundaryCD=="Bn0":
            fesOmega=H1(mesh, order=feOrder, definedon=reduced_region, dirichlet=Ht0_boundary)
 
        fes=fesA*fesOmega    
        (A, Omega),(N, o) = fes.TnT()
        normal = specialcf.normal(mesh.dim)
        a= BilinearForm(fes)
        a +=1/Mu*curl(A)*curl(N)*dx(total_region)
        a += -Mu*grad(Omega)*grad(o)*dx(reduced_region)
        a += -( Omega*(curl(N).Trace()*normal) + o*(curl(A).Trace()*normal) ) *ds(total_boundary)

        f = LinearForm(fes)
        f +=-o*(Bv*normal)  *ds(total_boundary)
        f += Cross(N.Trace(),Hs)*normal*ds(total_boundary)
        with TaskManager():
            f.Assemble()

        gf=self.Solve(fes, a, f)
        """
        with TaskManager():
            a.Assemble()
        gf=GridFunction(fes)
        gf=solver.iccg_solve(fes, gf, a, f.vec.FV(), tol=1.e-8, max_iter=1000, accel_factor=0, divfac=100, diviter=100,
                     scaling=True, complex=False,logplot=True)
        """
        gfA, gfOmega=gf.components
        BField=mu0*grad(gfOmega)+curl(gfA)+Bs

        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        print("feOrder=", feOrder,"  ", "ndof=",fes.ndof,"  ")
        self.CalcResult(model, BField)
        """
        mip = mesh(0,0,0)
        print("center magnetic field = ", BField(mip),"  ")
        Wm=Integrate(BField*BField/Mu*dx("iron"), mesh)
        print("magnetic energy=", Wm,"  ")
        print(f"経過時間: {elapsed_time:.4f} 秒  ")

        from ngsolve.webgui import Draw
        #print("**** Omega field ****")
        #Draw (gfOmega, mesh, order=feOrder, deformation=False) 
        print("**** B field ****")
        Draw (BField, mesh, order=feOrder, min=0., max=5.0, deformation=False) 
        """
        print(f"経過時間: {elapsed_time:.4f} 秒  ")
    """
    def CalcResult(self):
        BField=self.BField
        mesh=self.mesh
        mip = self.mesh(0,0,0)
        print("center magnetic field = ", BField(mip),"  ")
        Wm=Integrate(BField*BField/self.Mu*dx("iron"), mesh)
        print("magnetic energy=", Wm,"  ")


        from ngsolve.webgui import Draw
        #print("**** Omega field ****")
        #Draw (gfOmega, mesh, order=feOrder, deformation=False) 
        print("**** B field ****")
        Draw (BField, mesh, order=self.feOrder, min=0., max=5.0, deformation=False) 
    """
        