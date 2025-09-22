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

class A_Phi_ReducedA_Method(Static_Method):
    def __init__(self,  model,   **kwargs): 
        super().__init__(model,  **kwargs)

    def Ht_regularization(self, H, mesh, boundary):
        normal = specialcf.normal(mesh.dim)
        fesu = H1(mesh, order=self.feOrder, definedon=mesh.Boundaries(self.total_boundary), complex=False)
        u, v= fesu.TnT()
        a = BilinearForm(fesu)
        a +=Cross(normal,grad(u).Trace())*Cross(normal, grad(v).Trace())*ds
        f=LinearForm(fesu)
        f +=-H*Cross(normal, grad(v).Trace())*ds
        with TaskManager():
            a.Assemble()
            f.Assemble()
        gfu=GridFunction(fesu)
        gfu=solver.iccg_solve(fesu, gfu, a, f.vec.FV(), tol=1.e-16, max_iter=200, accel_factor=0, complex=False, logplot=True) 
        hreg=H+Cross(normal,grad(gfu).Trace())
        return hreg   
        
    def Calc(self,  **kwargs):
        default_values = {"feOrder":1,
                          "boundaryCD":"Bn0", 
                          "Kelvin": "off",
                          "regularization": False,
                          "tol":1.e-8
                         }
        default_values.update(kwargs)
        self.feOrder=default_values["feOrder"]
        boundaryCD=default_values["boundaryCD"]
        Kelvin=default_values["Kelvin"]
        regularization=default_values["regularization"]
        tol=default_values["tol"]
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
        Ht0_boundary=model.Ht0_boundary
        self.Mu=model.Mu
        Mu=self.Mu
        if self.jomega:
            self.Sigma=model.Sigma

        self.total_region=total_region
        self.reduced_region=reduced_region
        self.total_boundary=total_boundary
        self.Ht0_boundary=Ht0_boundary
        
        #coil=UNIF(0,0,1,0)
        coil=model.coil.field
        #coil=model.coil
        Av=Afield(coil)
        Bv=Bfield(coil)
        mu0=4.e-7*math.pi
        Hv=Bv/mu0
        zero=(0,0,0)

        Bs_dic = {"iron":zero, "conductor":zero, "air":zero, "hole":zero, "A_domain":zero, "Omega_domain":Bv, "reduced_region":Bv, "Kelvin":zero, 'default':zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in mesh.GetMaterials()])
        As_dic = {"iron":zero, "conductor":zero, "air":zero, "hole":zero, "A_domain":zero, "Omega_domain":Av, "reduced_region":Av, "Kelvin":zero, 'default':zero}
        As = CoefficientFunction([As_dic[mat] for mat in mesh.GetMaterials()])

        if Kelvin=="off":
            if boundaryCD=="Bn0":
                fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary+'|'+reduced_boundary, complex=self.jomega)
            elif boundaryCD=="Ht0":
                fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary, complex=self.jomega) 

        else:
            fesA=HCurl(mesh, order=feOrder, type1=True, nograds=True, dirichlet=Bn0_boundary, complex=self.jomega) 
            fesA=Periodic(fesA)

        fesPhi=H1(mesh, order=feOrder, definedon=model.conductive_region, dirichlet=Bn0_boundary, complex=self.jomega) 
        fes=fesA*fesPhi
        (A,phi), (N, psi) = fes.TnT() 

        dofLimit=1.0e6
        if fes.ndof > dofLimit:
            print(" fespace Dof >  dofLimit, DOF=", fes.ndof)
            return 0


        self.fes=fes
        (A,phi), (N, psi) = fes.TnT() 
        gf = GridFunction(fes)
        gfA, gfPhi=gf.components

        normal = specialcf.normal(mesh.dim)
        a= BilinearForm(fes)
        a +=1/Mu*curl(A)*curl(N)*dx(total_region)
        a +=1/Mu*curl(A)*curl(N)*dx(reduced_region)

        if self.jomega:
            a +=self.s*model.Sigma*(A+grad(phi))*(N+grad(psi))*dx(model.conductive_region)

        if Kelvin=="on":
            rs=self.model.rKelvin
            xs=self.model.kcenter
            r=sqrt((x-xs)*(x-xs)+y*y+z*z)
            fac=rs*rs/r/r
            a +=1/(Mu*fac) *curl(A)*curl(N)*dx("Kelvin")

        
        # Calculate Dirichlet condition terms
        gfA.Set(Av, BND, mesh.Boundaries(total_boundary))
        #Draw(gfA, mesh)
        f = LinearForm(fes)
        f +=1/Mu*curl(gfA)*curl(N)*dx(reduced_region)
        
        with TaskManager():
            f.Assemble()    
        #remove components of the Dirichlet boundary
        fcut = np.array(f.vec.FV())[fes.FreeDofs()]
        np.array(f.vec.FV(), copy=False)[fes.FreeDofs()] = fcut
        
        # Add Neumann condition terms
        #Draw(Hv, mesh)
        if regularization:
            hst=self.Ht_regularization(Hv, mesh, total_boundary)
            f += Cross(N.Trace(),hst)*normal*ds(total_boundary)
        else:
            f += Cross(N.Trace(),Hv)*normal*ds(total_boundary)
        with TaskManager():
            f.Assemble()

        gf=self.Solve(fes, a, f, tol)
        gfA, gfPhi=gf.components
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
        #print("*** Ar ***")
        #Draw(Ar, mesh)
        
        Br=curl(Arr)-curl(Axr)
        #Draw(Br.real, mesh)
        #Draw(Br.imag, mesh)
        BField=Bt+Br+Bs
        self.BField=BField
        self.Br=Br
        self.Bt=Bt
        #Draw(BField.real, mesh)

        if self.jomega: 
            self.JField=-self.s*model.Sigma*(At+grad(gfPhi))

        else:
            JField=0
    
        #self.CalcResult(model, BField, self.JField)
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
        return 1


    def CalcError(self):
        mesh=self.mesh
        feOrder=self.feOrder
        total_region=self.total_region
        reduced_region=self.reduced_region
        Mu=self.Mu

        fesH = HCurl(mesh, order=feOrder-1, type1=True, nograds=True, definedon=total_region,  dirichlet=self.Ht0_boundary, complex=self.jomega)
        Hd=GridFunction(fesH, "flux")
        H = 1/Mu*self.Bt
        Hd.Set(H)
        dH=H-Hd
        if self.jomega:
            err=Mu*(dH.real*dH.real+dH.imag*dH.imag)/4
        else:
            err=Mu*dH*dH/2
        eta_t = Integrate(err, mesh, VOL, definedon=mesh.Materials(total_region), element_wise=True)
       
    
        if self.jomega:
            fesJ = HDiv(mesh, order=feOrder-1,  definedon=self.model.conductive_region, complex=self.jomega)
            Jd=GridFunction(fesJ)
            J=self.JField
            Jd.Set(J)
            dJ=J-Jd
            err=self.model.Rho/self.omega*(dJ.real*dJ.real+dJ.imag*dJ.imag)/2
            eta_j = Integrate(err, mesh, VOL, definedon=mesh.Materials(self.model.conductive_region), element_wise=True)
            print("max(eta_t)=", max(eta_t), "  max(eta_j)=", max(eta_j)) 
            eta_t=eta_t+eta_j
          
        #gfur=self.Or
        fesH = HCurl(mesh, order=feOrder-1, type1=True, nograds=True,  definedon=reduced_region, dirichlet=self.Ht0_boundary, complex=self.jomega)
        Hd=GridFunction(fesH, "flux")
        H = 1/Mu*self.Br
        Hd.Set(H)
        dH=H-Hd

        if self.jomega:
            err=Mu*(dH.real*dH.real+dH.imag*dH.imag)/4
        else:
            err=Mu*dH*dH/2
        eta_r = Integrate(err, mesh, VOL, definedon=mesh.Materials(reduced_region), element_wise=True) 
        """
        print("eta_r", eta_r)
        for n in range(len(eta_r)):
            print(n, eta_r[n])
        """

        if self.Kelvin=="on":
            fesH = HCurl(mesh, order=feOrder-1, type1=True, nograds=True,  definedon="Kelvin", dirichlet=self.Ht0_boundary, complex=self.jomega)
            Hd=GridFunction(fesH, "flux")
            
            rs=self.model.rKelvin
            xs=self.model.kcenter
            r=sqrt((x-xs)*(x-xs)+y*y+z*z)
            fac=rs*rs/r/r
            H=1/(Mu*fac)*self.Br
            Hd.Set(H)
            dH=H-Hd

            if self.jomega:
                err=Mu*(dH.real*dH.real+dH.imag*dH.imag)/4

            else:
                err=Mu*dH*dH/2

            eta_k =Integrate(err, mesh,  VOL, definedon=mesh.Materials("Kelvin"), element_wise=True) 
            print(" maxerr = ", max(eta_t), max(eta_r), max(eta_k) )
            eta_r = eta_r+eta_k

        error=eta_t+eta_r
        maxerr = max(error)
        print ("ndof=", self.fes.ndof, " maxerr = ", maxerr)
        return maxerr, error
