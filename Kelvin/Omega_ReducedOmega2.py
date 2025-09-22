import math
from ngsolve import *
from ngsolve.webgui import Draw
import numpy as np
import sys
sys.path.append('..\include')
from MatrixSolver import MatrixSolver as solver 
sys.path.append(r'..\bin\Release') 
from EMPY_Field import *
sys.path.append(r'..\Static') 
from Static_Method import Static_Method
from netgen.occ import *
class Omega_ReducedOmega_Method2(Static_Method):
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
        self.mesh=mesh
        self.rsurf=model.rsurf
        
        import time
        start_time = time.perf_counter()


        total_region=model.total_region
        reduced_region=model.reduced_region
        total_boundary=model.total_boundary
        reduced_boundary=model.reduced_boundary
        Bn0_boundary=model.Bn0_boundary
        Ht0_boundary=model.Ht0_boundary
        Kelvin_region=model.Kelvin_region

        R0=5
        z0=10
        mu0=4.e-7*math.pi
        mu=mu0*1000
        mu_d={"iron":mu,  "A_domain":mu0, "Omega_domain":mu0,"Kelvin":mu0, 'default':mu0}
        Mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])  # デフォルトの物性値
        self.Mu=Mu
        
        self.total_region=total_region
        self.reduced_region=reduced_region
        


        #field=UNIF(0,0,1,0)
        Ov=Ofield(coil)
        Bv=Bfield(coil)
        mu0=4.e-7*math.pi
        Hv=Bv/mu0
        zero=(0,0,0)
        Bs_dic = {"iron":zero, "A_domain":zero, "Omega_domain":Bv, "Kelvin":Bv,'default':zero}
        Bs = CoefficientFunction([Bs_dic[mat] for mat in mesh.GetMaterials()])
        Os_dic = {"iron":zero, "A_domain":zero, "Omega_domain":Ov, "Kelvin":Bv, 'default':zero}
        Os = CoefficientFunction([Os_dic[mat] for mat in mesh.GetMaterials()])

        if boundaryCD=="Ht0":
            fes=H1(mesh, order=feOrder, dirichlet=reduced_boundary+"|"+Ht0_boundary)
        elif boundaryCD=="Bn0":
            fes=H1(mesh, order=feOrder, dirichlet=Ht0_boundary)
  
        self.fes=fes
        omega,psi = fes.TnT() 
        gfOmega = GridFunction(fes)
        a= BilinearForm(fes)
        a +=Mu*(grad(omega)*grad(psi))*dx(total_region)
        a +=Mu*(grad(omega)*grad(psi))*dx(reduced_region)

        H=grad(omega)
        B=self.HtoB_Kelvin(H)
        
        a +=grad(psi)*B*dx(Kelvin_region)
        
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
        
        gfOmega=self.Solve(fes, a, f)
        """
        gfOmega = GridFunction(fes)   #Clear gfOmega
        gfOmega=solver.iccg_solve(fes, gfOmega, a, f.vec.FV(), tol=1.e-8, max_iter=1000, accel_factor=0, divfac=100, diviter=10,
                     scaling=True, complex=False,logplot=True)
        """
        fesOt=H1(mesh, order=feOrder, definedon=total_region)
        fesOr=H1(mesh, order=feOrder, definedon=reduced_region+'|'+'Kelvin')
        Ot=GridFunction(fesOt)
        Orr=GridFunction(fesOr)
        Oxr=GridFunction(fesOr)

        Ot.Set(gfOmega,VOL, definedon=total_region)
        Orr.Set(gfOmega,VOL, definedon=reduced_region+'|'+'Kelvin')
        #Oxr.Set(surfaceOmega, BND, mesh.Boundaries(total_boundary))
        Oxr.Set(Ov, BND, mesh.Boundaries(total_boundary))

        Bt=grad(Ot)*Mu
        self.Orr=Orr
        Or=Orr-Oxr
        """
        xp=0
        yp=0
        zp=1.5
        for n in range(100):
            mip = mesh(xp,yp,zp)
            print(zp, "  ", Or(mip))
            zp=zp+(self.rsurf*2-1.5)/100
        """            
        
        Br=(grad(Orr)-grad(Oxr))*mu0
        BField=Bt+Br+Bs

        print("feOrder=", feOrder,"  ", "ndof=",fes.ndof,"  ")
        #self.CalcResult(model, BField)

        #Draw (Ot, mesh, order=3)  
        Draw (Or, mesh, order=3) 
        
        self.Br=Br
        self.Bt=Bt
        #Draw ((Ot+Or+Os)*Mu, mesh, order=3, min=0., max=1.0, deformation=False)  
   
        mip = mesh(0,0,0)
        print("center magnetic field = ", BField(mip),"  ")
        Wm=Integrate(BField*BField/Mu*dx("iron"), mesh)
        print("magnetic energy=", Wm,"  ")


        #Draw ((Ot+Or+Os)*mu, mesh, order=3, min=0., max=1.0, deformation=False)  
        print("**** B field ****")
        Draw (BField, mesh, order=feOrder, min=0, max=5, deformation=False)
        

        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f"経過時間: {elapsed_time:.4f} 秒  ")

    def CalcError(self):
        mesh=self.mesh
        feOrder=self.feOrder
        total_region=self.total_region
        reduced_region=self.reduced_region
        Mu=self.Mu
        """
        #gfu=self.Ot
        fesBd = HDiv(mesh, order=feOrder-1, definedon=total_region)
        Bd=GridFunction(fesBd, "flux")
        B = self.Bt
        Bd.Set(B)
        err=1/Mu*(B-Bd)*(B-Bd)
        eta = Integrate(err, mesh, VOL, element_wise=True)

        #gfur=self.Or
        fsBdr = HDiv(mesh, order=feOrder-1, definedon=reduced_region)
        Bdr=GridFunction(fsBdr, "flux")
        Br = self.Br
        Bdr.Set(Br)
        err_r=1/Mu*(Br-Bdr)*(Br-Bdr)
        eta_r = Integrate(err_r, mesh, VOL, element_wise=True) 
        """
        fsK = HDiv(mesh, order=feOrder-1, definedon="Kelvin")
        Bdk=GridFunction(fsK, "flux")
        
        fesOr=H1(mesh, order=feOrder, definedon="Kelvin")
        Ok=GridFunction(fesOr)
        Ok.Set(self.Orr)
        Hk = grad(Ok)
        Bk=self.HtoB_Kelvin(Hk)
        Bdk.Set(Bk)
        error_k=(Bk-Bdk)*self.BtoH_Kelvin(Bk-Bdk)
        eta_k = Integrate(error_k, mesh,  VOL, element_wise=True) 
        print("eta_k=", eta_k)
        #error=eta+eta_r+eta_k
        error=eta_k
        maxerr = max(error)
        print ("ndof=", self.fes.ndof, " maxerr = ", maxerr)
        return maxerr, error


    def HtoB_Kelvin(self, H):
        r=sqrt(x*x+y*y+z*z)
        nvector=CoefficientFunction( (x,y,z) )/r
        Hn=(H*nvector)*nvector
        Ht=H-Hn

        rs=self.rsurf
        rinf=rs*2
        mu0=4.e-7*math.pi
        mun=(rinf-rs)*rs*r/(rinf-r)/(rinf-r)/(rinf-r)*mu0
        mut=(rinf-rs)*rs/(rinf-r)/(rinf-r)*mu0
        B=mun*Hn + mut*Ht
        return B

    def BtoH_Kelvin(self, B):
        r=sqrt(x*x+y*y+z*z)
        nvector=CoefficientFunction( (x,y,z) )/r
        Bn=(B*nvector)*nvector
        Bt=B-Bn

        rs=self.rsurf
        rinf=rs*2
        mu0=4.e-7*math.pi
        mun=(rinf-rs)*rs*r/(rinf-r)/(rinf-r)/(rinf-r)*mu0
        mut=(rinf-rs)*rs/(rinf-r)/(rinf-r)*mu0
        H=Bn/mun + Bt/mut
        return H
        

        