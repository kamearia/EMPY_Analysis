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

class Omega_ReducedOmega_Method(Static_Method):
    def __init__(self,  model, coil,  **kwargs):  
        super().__init__(model, coil,  **kwargs)
        
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
        Ht0_boundary=model.Ht0_boundary
        self.Mu=model.Mu
        Mu=self.Mu

        #field=UNIF(0,0,1,0)
        Ov=Ofield(coil)
        Bv=Bfield(coil)
        mu0=4.e-7*math.pi
        Hv=Bv/mu0
        zero=(0,0,0)
        Bs_dic = {"iron":zero, "A_domain":zero, reduced_region:Bv, 'default':zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in mesh.GetMaterials()])
        Os_dic = {"iron":zero, "A_domain":zero, reduced_region:Ov, 'default':zero}
        Os = CoefficientFunction([Os_dic[mat] for mat in mesh.GetMaterials()])

        if boundaryCD=="Ht0":
            fes=H1(mesh, order=feOrder, dirichlet=reduced_boundary+"|"+Ht0_boundary)
        elif boundaryCD=="Bn0":
            fes=H1(mesh, order=feOrder, dirichlet=Ht0_boundary)
    
        omega,psi = fes.TnT() 
        gfOmega = GridFunction(fes)
        a= BilinearForm(fes)
        a +=Mu*(grad(omega)*grad(psi))*dx
        with TaskManager():
            a.Assemble()
        normal = specialcf.normal(mesh.dim)

        #surfaceOmega=HtoOmega(mesh, total_boundary, feOrder, Hv)
        #surfaceOmega=y
        # Calculate Dirichlet condition terms
        gfOmega.Set(Ov, BND, mesh.Boundaries(total_boundary))
        #gfOmega.Set(surfaceOmega, BND, mesh.Boundaries(total_boundary))

        f = LinearForm(fes)
        f +=Mu*grad(gfOmega)*grad(psi)*dx(reduced_region)
        with TaskManager():
            f.Assemble() 
        #remove components of the Dirichlet boundary
        fcut = np.array(f.vec.FV())[fes.FreeDofs()]
        np.array(f.vec.FV(), copy=False)[fes.FreeDofs()] = fcut

        # Add Neumann condition terms
        f += (normal*Bv)*psi*ds(total_boundary)
        with TaskManager():
            f.Assemble()
        gfOmega = GridFunction(fes)   #Clear gfOmega
        gfOmega=solver.iccg_solve(fes, gfOmega, a, f.vec.FV(), tol=1.e-8, max_iter=1000, accel_factor=0, divfac=100, diviter=10,
                     scaling=True, complex=False,logplot=True)

        fesOt=H1(mesh, order=feOrder, definedon=total_region)
        fesOr=H1(mesh, order=feOrder, definedon=reduced_region)
        Ot=GridFunction(fesOt)
        Orr=GridFunction(fesOr)
        Oxr=GridFunction(fesOr)

        Ot.Set(gfOmega,VOL, definedon=total_region)
        Orr.Set(gfOmega,VOL, definedon=reduced_region)
        #Oxr.Set(surfaceOmega, BND, mesh.Boundaries(total_boundary))
        Oxr.Set(Ov, BND, mesh.Boundaries(total_boundary))

        Bt=grad(Ot)*Mu
        Or=Orr-Oxr
        Br=(grad(Orr)-grad(Oxr))*mu0
        BField=Bt+Br+Bs

        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print("feOrder=", feOrder,"  ", "ndof=",fes.ndof,"  ")
        self.CalcResult(model, BField)
        """   
        mip = mesh(0,0,0)
        print("center magnetic field = ", BField(mip),"  ")
        Wm=Integrate(BField*BField/Mu*dx("iron"), mesh)
        print("magnetic energy=", Wm,"  ")


        #Draw ((Ot+Or+Os)*mu, mesh, order=3, min=0., max=1.0, deformation=False)  
        print("**** B field ****")
        Draw (BField, mesh, order=feOrder, min=0, max=5, deformation=False)
        """

        print(f"経過時間: {elapsed_time:.4f} 秒  ")