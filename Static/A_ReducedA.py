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

class A_ReducedA_Method(Static_Method):
    def __init__(self,  model, coil,  **kwargs): 
        super().__init__(model, coil,  **kwargs)
        
    def Calc(self, model, coil, **kwargs):
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
        self.Mu=model.Mu
        Mu=self.Mu

        self.total_region=total_region
        self.reduced_region=reduced_region
        
        #coil=UNIF(0,0,1,0)
        Av=Afield(coil)
        Bv=Bfield(coil)
        mu=4.e-7*math.pi
        Hv=Bv/mu
        zero=(0,0,0)

        Bs_dic = {"iron":zero, "A_domain":zero, "Omega_domain":Bv, "Kelvin":zero, 'default':zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in mesh.GetMaterials()])
        As_dic = {"iron":zero, "A_domain":zero, "Omega_domain":Av, "Kelvin":zero, 'default':zero}
        As = CoefficientFunction([As_dic[mat] for mat in mesh.GetMaterials()])

        if Kelvin=="off":
            if boundaryCD=="Bn0":
                fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary+'|'+reduced_boundary)
            elif boundaryCD=="Ht0":
                fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary) 

        else:
            fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary) 
            fesA=Periodic(fesA)
            
        self.fes=fesA
        A,N = fesA.TnT() 
        gfA = GridFunction(fesA)
        normal = specialcf.normal(mesh.dim)
        a= BilinearForm(fesA)
        a +=1/Mu*curl(A)*curl(N)*dx(total_region)
        a +=1/Mu*curl(A)*curl(N)*dx(reduced_region)

        if Kelvin=="on":
            rs=self.model.rKelvin
            xs=2*rs
            r=sqrt((x-xs)*(x-xs)+y*y+z*z)
            fac=rs*rs/r/r
            a +=1/(Mu*fac) *curl(A)*curl(N)*dx("Kelvin")

        
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

        gfA=self.Solve(fesA, a, f)
        """
        with TaskManager():
            a.Assemble()
        gfA = GridFunction(fesA)   #Clear gfA
        gfA=solver.iccg_solve(fesA, gfA, a, f.vec.FV(), tol=1.e-8, max_iter=1000, accel_factor=0, divfac=100, diviter=10,
                     scaling=True, complex=False,logplot=True)
        """

        fesAt=HCurl(mesh, order=feOrder, definedon=total_region, dirichlet=Bn0_boundary, nograds=True)
        fesAr=HCurl(mesh, order=feOrder, definedon=reduced_region+"|Kelvin", dirichlet=Bn0_boundary, nograds=True)
        At=GridFunction(fesAt)
        Arr=GridFunction(fesAr)
        Axr=GridFunction(fesAr)
        At.Set(gfA,VOL, definedon=total_region)
        Arr.Set(gfA,VOL, definedon=reduced_region)
        Axr.Set(Av, BND, mesh.Boundaries(total_boundary))

        Bt=curl(At)
        Ar=Arr-Axr
        print("*** Ar ***")
        Draw(Ar, mesh)
        
        Br=curl(Arr)-curl(Axr)
        BField=Bt+Br+Bs
        self.Br=Br
        self.Bt=Bt
        


        print("feOrder=", feOrder,"  ", "ndof=",fesA.ndof,"  ")
        self.CalcResult(model, BField)
        """
        mip = mesh(0,0,0)
        print("center magnetic field = ", BField(mip),"  ")
        Wm=Integrate(BField*BField/Mu*dx("iron"), mesh)
        print("magnetic energy=", Wm,"  ")


        print("**** B field ****")
        Draw (BField, mesh, order=feOrder, min=0., max=5.0, deformation=False) 
        """
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
        fsHdr = HCurl(mesh, order=feOrder-1, type1=True, nograds=True,  definedon=reduced_region)
        Hdr=GridFunction(fsHdr, "flux")
        Hr = 1/Mu*self.Br
        Hdr.Set(Hr)
        err_r=Mu*(Hr-Hdr)*(Hr-Hdr)
        eta_r = Integrate(err_r, mesh, VOL, element_wise=True) 

        if self.Kelvin=="on":
            fkBdr = HCurl(mesh, order=feOrder-1, type1=True, nograds=True,  definedon="Kelvin")
            Hdk=GridFunction(fkBdr, "flux")
            
            rs=self.model.rKelvin
            xs=2*rs
            r=sqrt((x-xs)*(x-xs)+y*y+z*z)
            fac=rs*rs/r/r
            Hk=1/(Mu*fac)*self.Br
            Hdk.Set(Hk)
            err_rk=(Mu*fac)*(Hk-Hdk)*(Hk-Hdk)
            eta_rk =Integrate(err_rk, mesh, VOL, element_wise=True) 
            print(" maxerr = ", max(eta), max(eta_r), max(eta_rk) )
            eta_r = eta_r+eta_rk

        error=eta+eta_r
        maxerr = max(error)
        print ("ndof=", self.fes.ndof, " maxerr = ", maxerr)
        return maxerr, error
