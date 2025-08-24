import math
from ngsolve import *
from ngsolve.webgui import Draw
import numpy as np
import sys
sys.path.append('..\include')
from MatrixSolver import MatrixSolver as solver 
sys.path.append(r'..\bin\Release') 
from EMPY_Field import *
class A_ReducedA_Method():
    def __init__(self,  model, coil,  **kwargs):  
        self.Calc(model, coil, **kwargs)
        
    def Calc(self, model, coil, **kwargs):
        default_values = {"feOrder":1,
                          "boundaryCD":"Bn0", 
                         }
        default_values.update(kwargs)
        self.feOrder=default_values["feOrder"]
        boundaryCD=default_values["boundaryCD"]

        feOrder=self.feOrder
        mesh=model.mesh
        
        import time
        start_time = time.perf_counter()

        total_region=model.total_region
        reduced_region=model.reduced_region
        total_boundary=model.total_boundary
        reduced_boundary=model.reduced_boundary
        Bn0_boundary=model.Bn0_boundary
        Mu=model.Mu
        
        #coil=UNIF(0,0,1,0)
        Av=Afield(coil)
        Bv=Bfield(coil)
        mu=4.e-7*math.pi
        Hv=Bv/mu
        zero=(0,0,0)
        As_dic = {"iron":zero, "A_domain":zero, reduced_region:Av, 'default':zero}
        As = CoefficientFunction([As_dic[mat] for mat in mesh.GetMaterials()])
        Bs_dic = {"iron":zero, "A_domain":zero, reduced_region:Bv, 'default':zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in mesh.GetMaterials()])

        if boundaryCD=="Bn0":
            fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary+'|'+reduced_boundary)
        elif boundaryCD=="Ht0":
            fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary) 
    
        A,N = fesA.TnT() 
        gfA = GridFunction(fesA)
        normal = specialcf.normal(mesh.dim)
        a= BilinearForm(fesA)
        a +=1/Mu*curl(A)*curl(N)*dx

        
        # Calculate Dirichlet condition terms
        gfA.Set(Av, BND, mesh.Boundaries(total_boundary))
        f = LinearForm(fesA)
        f +=1/Mu*curl(gfA)*curl(N)*dx(reduced_region)
        with TaskManager():
            f.Assemble()    
        #remove components of the Dirichlet boundary
        fcut = np.array(f.vec.FV())[fesA.FreeDofs()]
        np.array(f.vec.FV(), copy=False)[fesA.FreeDofs()] = fcut

        # Add Neumann condition terms
        f += Cross(N.Trace(),Hv)*normal*ds(total_boundary)
        with TaskManager():
            f.Assemble()

        with TaskManager():
            a.Assemble()
        gfA = GridFunction(fesA)   #Clear gfA
        gfA=solver.iccg_solve(fesA, gfA, a, f.vec.FV(), tol=1.e-8, max_iter=1000, accel_factor=0, divfac=100, diviter=10,
                     scaling=True, complex=False,logplot=True)

        fesAt=HCurl(mesh, order=feOrder, definedon=total_region, dirichlet=Bn0_boundary, nograds=True)
        fesAr=HCurl(mesh, order=feOrder, definedon=reduced_region, dirichlet=Bn0_boundary, nograds=True)
        At=GridFunction(fesAt)
        Arr=GridFunction(fesAr)
        Axr=GridFunction(fesAr)
        At.Set(gfA,VOL, definedon=total_region)
        Arr.Set(gfA,VOL, definedon=reduced_region)
        Axr.Set(Av, BND, mesh.Boundaries(total_boundary))

        Bt=curl(At)
        Ar=Arr-Axr
        Br=curl(Arr)-curl(Axr)
        BField=Bt+Br+Bs

        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        print("feOrder=", feOrder,"  ", "ndof=",fesA.ndof,"  ")
        mip = mesh(0,0,0)
        print("center magnetic field = ", BField(mip),"  ")
        Wm=Integrate(BField*BField/Mu*dx("iron"), mesh)
        print("magnetic energy=", Wm,"  ")
        print(f"経過時間: {elapsed_time:.4f} 秒  ")

        print("**** B field ****")
        Draw (BField, mesh, order=feOrder, min=0., max=5.0, deformation=False) 
