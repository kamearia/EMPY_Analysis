import numpy as np
import math
from ngsolve import *

#import include.MatrixSolver as solver
import os, sys
sys.path.append('/home/jupyter-ksugahar/亀有/include')
from EM_Model import EM_Model
from MatrixSolver import MatrixSolver as solver

class T_Omega_Method(EM_Model):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        default_values = {"Bz0": 1 ,     
                          "boxz":400e-3
                         }
        default_values.update(kwargs)
        self.Bz0=default_values["Bz0"]
        mu=4.e-7*math.pi
        self.Omega0=self.Bz0*default_values["boxz"]*2/mu
        print("Bz0= ",  self.Bz0, "  Omega0= ", self.Omega0)
        
    def LoopFields(self, domain, **kwargs):
        mesh=self.GetMesh()
        default_values = {"connected": 1}
        default_values.update(kwargs)
        connected=default_values["connected"]
    
        #smesh=self.surface_mesh( domain )
        smesh=self.surface_mesh_from_boundary("conductorBND")
        g=self.surface_genus(smesh, connected)
        print("genus= ", g)
        fes = HCurl(mesh, order=1, nograds=True, definedon=domain)
        u,v = fes.TnT()
        fesPhi = H1(mesh, order=1, definedon="air")
        phi,psi= fesPhi.TnT()
        
        loops=[]
        for k in range(g):
            gfu = GridFunction(fes)
            id=self.random_edge(smesh)
            edge_dofs = fes.GetDofNrs(id)[0]   
            
            #edge_dofs=427
            
            gfu.vec[:] = 0
            gfu.vec[edge_dofs] = 1
            print("edge DOF = ", edge_dofs)

            fes.FreeDofs().__setitem__(edge_dofs,False)

            a = BilinearForm(fes)
            a += curl(u)*curl(v)*dx
            f=LinearForm(fes)
            with TaskManager():
                a.Assemble()
                f.Assemble()
            fr=-a.mat*gfu.vec

            gfu=solver.iccg_solve(fes, gfu, a, fr.Evaluate(), tol=1.e-16, max_iter=200, accel_factor=1.1)

            gfPhi = GridFunction(fesPhi)
            a = BilinearForm(fesPhi)
            a += grad(phi)*grad(psi)*dx
            f=LinearForm(fesPhi)
            f += grad(psi)*gfu*dx
            with TaskManager():
                a.Assemble()
                f.Assemble()
    
            gfPhi=solver.iccg_solve(fesPhi, gfPhi, a, f.vec.FV(), tol=1.e-16, max_iter=200, accel_factor=1.0)  
            gfw=gfu-grad(gfPhi)
            a = BilinearForm(fes)
            a += u*v*dx
            f=LinearForm(fes)
            f += gfw*v*dx
            with TaskManager():
                a.Assemble()
                f.Assemble()
            gfw= GridFunction(fes)
            gfw=solver.iccg_solve(fes, gfw, a, f.vec.FV(), tol=1.e-16, max_iter=200, accel_factor=1.0) 
            
            gft=gfw
            for kd in range(len(loops)):
                prod=Integrate(gfw*loops[kd]*dx, mesh)
                for i in range(len(gft.vec)):
                    gft.vec[i] -=prod*loops[kd].vec[i]
      
            norm2=Integrate(gft*gft*dx, mesh)
            norm=sqrt(norm2)
            for i in range(len(gft.vec)):
                gft.vec[i]/=norm
              
            loops.append(gft)       
        return loops
        
    def SetLoopFields(self):
        self.loopFields=self.LoopFields("air" )
#        self.loopField=self.loopFields[0]
    def GetLoopFields(self):
        return self.loopFields
#    def GetLoopField(self, n):
#        return self.loopFields[n]
        
    def SetFeSpace(self):
        mesh=self.GetMesh()
        fesT = HCurl(mesh, order=self.feOrder, nograds=True, definedon=self.model.sigmaDomain,
                     dirichlet="conductorBND", complex=False)
        fesOmega = H1(mesh, order=self.feOrder, dirichlet="upper|lower",complex=False)
        fespace=fesT*fesOmega
        self.fespace=fespace

    def CalcRL(self, gfTOmega, loopField):
        gfT, gfOmega = gfTOmega.components
        mesh=gfTOmega.space.mesh
        s=self.GetLaplace_s()
        #mu=self.materials["air"].mu
        mu=self.model.mu
        #sig=self.materials["sig"].sig
        sigmaDomain=self.model.sigmaDomain
        sig=self.model.sigma
    
        R=Integrate(1/sig*curl(gfT)**2*dx(sigmaDomain) , mesh )
        L=Integrate(mu*(gfT+grad(gfOmega)+loopField)**2 *dx , mesh )
        return L,R
    
    def TroughCurrent(self, gf, loop):
        mesh=self.GetMesh()
        current=0.
        it=mesh.BBoundaries(loop).Elements()
        for e in it:
            start=e.vertices[0].nr
            end=e.vertices[1].nr
            dir=1
            if start> end: dir=-1
            for i in range(len(e.edges)): 
                k=gf.space.GetDofNrs(e.edges[i])
                current +=gf.vec[k[0]]*dir
        return current
    
    def LoopSelfTerm(self, loopField, gfTOmega):
        gfT, gfOmega = gfTOmega.components
        mesh=gfTOmega.space.mesh
        s=self.GetLaplace_s()
        mu=self.materials["air"].mu
        sig=self.materials["sig"].sig
        faf = Integrate(s*mu*gfT**2*dx("sig"), mesh)   
        faf +=Integrate(1./sig*curl(gfT)*curl(gfT)*dx("sig"), mesh)
        faf +=Integrate(s*mu*loopField*loopField*dx("air"), mesh) 
        return faf  
    
    def GetLoopCouplings(self, loopFields, fesTOmega, boundary):
        s=self.GetLaplace_s()
        mu=self.model.mu
        sig=self.model.sigma
        g=len(loopFields)
        fs=[]
        fafs=[]
        gfTs=[]
        fes=fesTOmega
        (T,omega),(W,psi) = fes.TnT()
        for k in range(g):
            loopField=loopFields[k]
            gfTOmega=GridFunction(fes)
            self.SetBoundaryValue(loopField, 1, gfTOmega, boundary )
            gfT, gfOmega = gfTOmega.components
            gfTs.append(gfTOmega)
            f=LinearForm(fes)
            f += 1./sig*curl(gfT)*curl(W)*dx("sig")
            f += s*mu*gfT*(W+grad(psi))*dx("sig")
            f += s*mu*loopField*grad(psi)*dx("air")
            with TaskManager():
                f.Assemble()
            fs.append(f)           
             
            mesh=self.GetMesh()
            tmp=[]
            for k2 in range(k+1):
                loopField2=loopFields[k2]
                gfT2, gfOmega2 = gfTs[k2].components
                #print("plot in GetLoopCouplings")
                #from ngsolve.webgui import Draw
                #Draw (gfT2, self.GetMesh(), order=3, min=0, max=2, deformation=False)
                faf =Integrate(1./sig*curl(gfT)*curl(gfT2)*dx("sig"), mesh) 
                faf +=Integrate(s*mu*gfT*gfT2*dx("sig"), mesh)
                faf +=Integrate(s*mu*loopField*loopField2*dx("air"), mesh) 
                tmp.append(faf)
                if k2<k: fafs[k2].append(faf)
            fafs.append(tmp)
            
        return fs, fafs
    
    def SetBoundaryValue(self, sr, fac, gf, boundaryId):
        mesh=sr.space.mesh
        for t in mesh.Boundaries(boundaryId).Elements():
            for e in t.edges:
                k=sr.space.GetDofNrs(e)
                k2=gf.space.GetDofNrs(e)
                gf.vec[k2[0]] = sr.vec[k[0]]*fac
    
    def SetBoundaryValue2(self, sr, fac, gf, boundaryId):
        nloop=len(sr)
        if nloop !=0:
            mesh=sr[0].space.mesh
            for t in mesh.Boundaries(boundaryId).Elements():
                for e in t.edges:
                    k2=gf.space.GetDofNrs(e)
                    sum=0
                    for n in range(nloop):
                        k=sr[n].space.GetDofNrs(e)
                        sum += sr[n].vec[k[0]]*fac[n]
                    gf.vec[k2[0]] = sum
        
    def CalcSourceTerm(self):
        a=self.systemMatrix
        gfTOmega= GridFunction(fesTOmega)
        gfT, gfOmega=gfTOmega.components
        gfOmega.Set(self.Omega0, definedon=self.mesh.Boundaries("upper")) 
        self.sourceTerm=-a.mat*gfTOmega.vec
        
    def CalcBzGiven(self):   
        fesTOmega=self.GetFeSpace()
        gfTOmega = GridFunction(fesTOmega)
        loopFields=self.loopFields
        c, faf=self.GetLoopCouplings(loopFields, fesTOmega, "conductorBND")
        #gfT, gfOmega=gfTOmega.components
        #print(" Draw(gfT) in CalcBzGiven")
        #mesh=self.GetMesh()
        #from ngsolve.webgui import Draw
        #self.SetBoundaryValue(loopFields[0], 1, gfTOmega,  "conductorBND" )
        #Draw(gfT, mesh)
        
        voltage=None
        a=self.GetSystemMatrix()
        gfTOmega= GridFunction(fesTOmega)
        gfT, gfOmega=gfTOmega.components
        mu=4.e-7*math.pi
        self.Omega0=self.Bz0*self.model.up_down/mu
        print("Bz0= ",  self.Bz0, "  Omega0= ", self.Omega0)
        gfOmega.Set(self.Omega0, definedon=self.GetMesh().Boundaries("upper")) 
        
        #from ngsolve.webgui import Draw
        #Draw (gfOmega, self.GetMesh(), order=3, min=0, max=2, deformation=False)
        
        sourceTerm=-a.mat*gfTOmega.vec

        x, amp= solver.SolveCoupled2(fesTOmega, a, c, faf, sourceTerm.Evaluate(),  voltage)
        a=None
        c=None
        np.array(gfTOmega.vec.FV(), copy=False)[fesTOmega.FreeDofs()] +=x
        x=None
       
        self.SetBoundaryValue2(self.loopFields, amp, gfTOmega, "conductorBND")  
    
        print("amp =", amp)
        self.Hfield=gfT+grad(gfOmega)
        for n in range(len(loopFields)):
            self.Hfield += amp[n]*loopFields[n]
        self.Bfield=self.Hfield*self.model.mu
        self.JField=curl(gfT)
        self.Omega=gfOmega
    
    def CalcCurretGiven(self):
        fesTOmega=self.GetFeSpace()
        gfTOmega = GridFunction(fesTOmega)
        cs, fafs=self.GetLoopCouplings(self.loopFields, fesTOmega, "conductorBND") 
        cs[0].vec.__imul__(-1.)
        
        print("faf=",fafs[0]) 
        
        gfT, gfOmega = gfTOmega.components
        self.SetBoundaryValue(self.loopFields[0], 1, gfTOmega,  "conductorBND" )
        gfTOmega=solver.iccg_solve(fesTOmega, gfTOmega, self.systemMatrix, cs[0].vec.FV(), 
                                   tol=1.e-16, max_iter=200, accel_factor=1.1)
        
        mulf_u=np.dot(gfTOmega.vec.FV(),cs[0].vec.FV())
        snum=self.model.GetSeriesNum()
        print("mulf_u=", -mulf_u/snum)
        
        L, R =self.CalcRL(gfTOmega, self.loopFields[0])
        # R+s*L=-cs[0]*gfTOmega + faf[0][0]

        R/=snum
        L/=snum
        s=self.GetLaplace_s()
        print("z=", R+s*L)
        
        mesh=self.GetMesh()
        normal = specialcf.normal(mesh.dim)
        #sum=Integrate(curl(gfT)*(0,0,1)*ds("in"), mesh)
        sum=Integrate(curl(gfT).Trace()*normal*ds("in"), mesh)
        print("Current=", sum)    
        sum=self.TroughCurrent(gfT, "inedge")
        print("line sum=", sum)

        self.Hfield=gfT+grad(gfOmega)+self.loopField
        self.JField=curl(gfT)

 
    def CalcVoltageGiven(self):   
        fesTOmega=self.GetFeSpace()
        gfTOmega = GridFunction(fesTOmega)

        loopFields=self.loopFields
        loopField=loopFields[0]
        cs, fafs=self.GetLoopCouplings(loopFields, fesTOmega, "conductorBND") 

        voltage=[1]
        voltage_part=[1*self.model.GetPortion()]
        a=self.GetSystemMatrix()
        x, current= solver.SolveCoupled2(fesTOmega, a, cs, fafs, None,  voltage_part)
        a=None
        c=None
        np.array(gfTOmega.vec.FV(), copy=False)[fesTOmega.FreeDofs()] =x[0]
        x=None
        print("voltage=", voltage[0], " current=", current[0])
        z=voltage[0]/current[0]
        print("z=", z)
        
        gfT, gfOmega = gfTOmega.components  
        mesh=self.GetMesh()
        #from ngsolve.webgui import Draw
        #Draw (gfT+grad(gfOmega), mesh, order=3, min=0, max=20, deformation=False)
        #Draw (curl(gfT), mesh, order=3, min=0, max=20, deformation=False)        
        normal = specialcf.normal(mesh.dim)
        sum=Integrate(curl(gfT)*normal*ds("in"), mesh)
        print("sum=", sum)
        #from ngsolve.webgui import Draw
        #Draw (curl(gfT), mesh, order=3, min=0, max=20, deformation=False)
        
        self.SetBoundaryValue(loopField, current[0], gfTOmega, "conductorBND")  
        mesh=self.GetMesh()
        normal = specialcf.normal(mesh.dim)
        sum=Integrate(curl(gfT)*normal*ds("in"), mesh)
        print("Current by surface integration=", sum)    
        sum=self.TroughCurrent(gfT, "inedge")
        print("Current by line integration=", sum)
        
        self.Hfield=(gfT+grad(gfOmega))/current[0]+loopField
        self.JField=curl(gfT)/current[0]
        #from ngsolve.webgui import Draw
        #Draw (gfT+grad(gfOmega), mesh, order=3, min=0, max=20, deformation=False)
        
    @classmethod
    def Exec_T_Omega_Method(cls, model, **kwargs):
        from ngsolve.webgui import Draw

#        system=EM_Model(**kwargs)
        system=T_Omega_Method(**kwargs)
        system.print()
        mesh=model.GetMesh()  
        system.SetMeshedModel(model)
        Draw(mesh)
        
        #system.SetLoopField()
        #loopField=system.LoopField("air" )
        system.loopFields=system.LoopFields("air" , connected=1)
        loopField=system.loopFields[0]
        current=-system.TroughCurrent(loopField, "inedge")
        loopField.vec.__imul__(1./current)
        system.loopField=loopField
        system.loopFields=[loopField]
        #loopField=system.GetLoopField()
        Draw (loopField, mesh, order=3, min=0, max=20, deformation=False) 
        
        #system.SetFeSpace()
        mesh=system.GetMesh()
        fesT = HCurl(mesh, order=system.feOrder, nograds=True, definedon="sig",
                     dirichlet="conductorBND", complex=False)
        fesOmega = H1(mesh, order=system.feOrder, dirichlet="upper|lower", complex=False)
        fespace=fesT*fesOmega
        system.fespace=fespace
        
        system.CalcSystemMatrix()
        a=system.GetSystemMatrix()
        if system.Bz0 !=0:
            system.CalcBzGiven()
        else: 
            if system.Egiven ==True:
                system.CalcVoltageGiven()
            else:
                system.CalcCurretGiven()

        Draw (system.Bfield, mesh, order=3, min=0, max=20, deformation=False)
        Draw (system.JField, mesh, order=3, min=0, max=1.e4, deformation=False)
        Draw (system.Omega, mesh,  deformation=False)
        return system
    
    @classmethod
    def CylinderCalc(cls, **kwargs):
        from ngsolve.webgui import Draw
        from include.CylinderModel import CylinderModel
        
#        system=EM_Model(**kwargs)
        system=T_Omega_Method(**kwargs)
        system.print()
        model=CylinderModel(**kwargs)
        mesh=model.GetMesh()  
        system.SetMeshedModel(model)
        Draw(mesh)
        
        #system.SetLoopField()
        #loopField=system.LoopField("air" )
        system.loopFields=system.LoopFields("air" , connected=1)
        loopField=system.loopFields[0]
        current=system.TroughCurrent(loopField, "inedge")
        loopField.vec.__imul__(1./current)
        system.loopField=loopField
        system.loopFields=[loopField]
        #loopField=system.GetLoopField()
        Draw (loopField, mesh, order=3, min=0, max=20, deformation=False) 
        
        #system.SetFeSpace()
        mesh=system.GetMesh()
        fesT = HCurl(mesh, order=system.feOrder, nograds=True, definedon="sig",
                     dirichlet="conductorBND", complex=False)
        fesOmega = H1(mesh, order=system.feOrder, complex=False)
        fespace=fesT*fesOmega
        system.fespace=fespace
        
        system.CalcSystemMatrix()
        a=system.GetSystemMatrix()
        if system.Egiven ==True:
            system.CalcVoltageGiven()
        else:
            system.CalcCurretGiven()
            

        Draw (system.Hfield, mesh, order=3, min=0, max=2, deformation=False)
        Draw (system.JField, mesh, order=3, min=20, max=30, deformation=False)     
        return system
    