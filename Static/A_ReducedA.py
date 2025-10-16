import math
from ngsolve import *
from ngsolve.webgui import Draw
import numpy as np
import sys
sys.path.append(r'..\include')
from MatrixSolver import MatrixSolver as solver 
sys.path.append(r'..\bin\Release') 
from EMPY_Field import *
from Static_Method import Static_Method

class A_ReducedA_Method(Static_Method):
    def __init__(self,  model,  **kwargs): 
        super().__init__(model, **kwargs)
        
    def Calc(self, **kwargs):
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
        mesh=self.mesh
        
        import time
        start_time = time.perf_counter()

        model=self.model
        total_region=model.total_region
        reduced_region=model.reduced_region
        total_boundary=model.total_boundary
        reduced_boundary=model.reduced_boundary
        Bn0_boundary=model.Bn0_boundary
        self.Mu=model.Mu
        Mu=self.Mu
        if self.jomega:
            self.Sigma=model.Sigma

        self.total_region=total_region
        self.reduced_region=reduced_region
        
        #coil=UNIF(0,0,1,0)
        coil=model.coil.field
        Av=Afield(coil)
        Bv=Bfield(coil)
        mu0=4.e-7*math.pi
        Hv=Bv/mu0
        zero=(0,0,0)

        Bs_dic = {"iron":zero, "A_domain":zero, "Omega_domain":Bv, "Kelvin":zero, 'default':zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in mesh.GetMaterials()])
        As_dic = {"iron":zero, "A_domain":zero, "Omega_domain":Av, "Kelvin":zero, 'default':zero}
        As = CoefficientFunction([As_dic[mat] for mat in mesh.GetMaterials()])

        if Kelvin=="off":
            if boundaryCD=="Bn0":
                fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary+'|'+reduced_boundary, complex=self.jomega)
            elif boundaryCD=="Ht0":
                fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary, complex=self.jomega) 

        else:
            fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary, complex=self.jomega) 
            fesA=Periodic(fesA)

        self.fes=fesA
        A,N = fesA.TnT() 
        gfA = GridFunction(fesA)
        normal = specialcf.normal(mesh.dim)
        a= BilinearForm(fesA)
        a +=1/Mu*curl(A)*curl(N)*dx(total_region)
        a +=1/Mu*curl(A)*curl(N)*dx(reduced_region)

        if self.jomega:
            a +=self.s*model.Sigma*A*N*dx(model.conductive_region)

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
        """
        with TaskManager():
            f.Assemble()    
        #remove components of the Dirichlet boundary
        fcut = np.array(f.vec.FV())[fesA.FreeDofs()]
        np.array(f.vec.FV(), copy=False)[fesA.FreeDofs()] = fcut
        """
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

        fesAt=HCurl(mesh, order=feOrder, definedon=total_region, dirichlet=Bn0_boundary, nograds=True, complex=self.jomega)
        fesAr=HCurl(mesh, order=feOrder, definedon=reduced_region+"|Kelvin", dirichlet=Bn0_boundary, nograds=True, complex=self.jomega)
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
        Draw(BField, mesh)
        if self.jomega: 
            self.JField=-self.s*model.Sigma*At
        else:
            JField=0
    
        self.CalcResult(model, BField, self.JField)
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

        fesH = HCurl(mesh, order=feOrder-1, type1=True, nograds=True, definedon=total_region, complex=self.jomega)
        Hd=GridFunction(fesH, "flux")
        H = 1/Mu*self.Bt
        Hd.Set(H)
        dH=H-Hd
        if self.jomega:
            err=Mu*(dH.real*dH.real+dH.imag*dH.imag)/4
        else:
            err=Mu*dH*dH/2
        eta_t = Integrate(err, mesh, VOL, element_wise=True)

        if self.jomega:
            fesJ = HDiv(mesh, order=feOrder-1,  definedon=self.model.conductive_region, complex=self.jomega)
            Jd=GridFunction(fesJ)
            J=self.JField
            Jd.Set(J)
            dJ=J-Jd
            err=1/(self.Sigma*self.omega)*(dJ.real*dJ.real+dJ.imag*dJ.imag)/2
            eta_j = Integrate(err, mesh, VOL, element_wise=True)
            print("max(eta_t)=", max(eta_t), "  max(eta_j)=", max(eta_j)) 
            eta_t=eta_t+eta_j
            
        #gfur=self.Or
        fesH = HCurl(mesh, order=feOrder-1, type1=True, nograds=True,  definedon=reduced_region, complex=self.jomega)
        Hd=GridFunction(fesH, "flux")
        H = 1/Mu*self.Br
        Hd.Set(H)
        dH=H-Hd
        if self.jomega:
            err=Mu*(dH.real*dH.real+dH.imag*dH.imag)/4
        else:
            err=Mu*dH*dH/2
        eta_r = Integrate(err, mesh, VOL, element_wise=True) 

        if self.Kelvin=="on":
            fesH = HCurl(mesh, order=feOrder-1, type1=True, nograds=True,  definedon="Kelvin", complex=self.jomega)
            Hd=GridFunction(fesH, "flux")
            
            rs=self.model.rKelvin
            xs=2*rs
            r=sqrt((x-xs)*(x-xs)+y*y+z*z)
            fac=rs*rs/r/r
            H=1/(Mu*fac)*self.Br
            Hd.Set(H)
            dH=H-Hd
            if self.jomega:
                err=Mu*(dH.real*dH.real+dH.imag*dH.imag)/4
            else:
                err=Mu*dH*dH/2
            eta_k =Integrate(err, mesh, VOL, element_wise=True) 
            print(" maxerr = ", max(eta_t), max(eta_r), max(eta_k) )
            eta_r = eta_r+eta_k

        error=eta_t+eta_r
        maxerr = max(error)
        print ("ndof=", self.fes.ndof, " maxerr = ", maxerr)
        return maxerr, error
