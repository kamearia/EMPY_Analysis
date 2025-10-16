from ngsolve import *
import sys
sys.path.append(r'..\include')
from MatrixSolver import MatrixSolver as solver 
import math
from ngsolve.webgui import Draw

class Static_Method():
    def __init__(self,  model,  **kwargs):
        default_values = {"jomega":False,
                          "freq":0,
                          "dofLimit":1.e6,
                         }
        default_values.update(kwargs)

        self.jomega=default_values["jomega"]
        self.dofLimit=default_values["dofLimit"]
        if self.jomega:
            self.freq=default_values["freq"]
            self.omega=2.*math.pi*self.freq
            self.s=1.j*self.omega

        self.model=model
        self.mesh=model.mesh
        #self.Calc(model, coil,  **kwargs)


    #def Calc(self, model, coil,  **kwargs):
    #    return

    def CalcResult(self, model, plotBFieldonLine=False, drawFields=True, pltBField=True, useModelPlot=False):
        BField=self.BField
        JField=self.JField
        mesh=model.mesh
        mip = mesh(0,0,0)
        print("center magnetic field = ", BField(mip),"  ")

        if self.jomega==False:
            Wm=Integrate(BField*BField/self.Mu*dx("iron"), mesh)/2
            print(" Average magnetic energy in conductor=", Wm)
        else:
            Wm=Integrate(BField*BField/self.Mu/2.*dx(model.conductive_region), mesh)
            #WJ=Integrate((JField.real*JField.real+JField.imag*JField.imag)/model. Sigma*dx(model.conductive_region), mesh) /2
            WJ=Integrate((JField*JField)/model.Sigma*dx(model.conductive_region), mesh) 
            print(" Magnetic energy in conductor=", Wm, " Joule loss= ", WJ)
        
        from ngsolve.webgui import Draw
        #print("**** Omega field ****")
        #Draw (gfOmega, mesh, order=feOrder, deformation=False) 
        if drawFields:
            if self.jomega==False:
                clipping_plane= (0, 1, 0, 0)
                Draw (BField, mesh, order=self.feOrder, min=0., max=5.0, deformation=False, 
                      clipping = { "pnt" : (0.0,0.001,0), "vec" : (0,1,0) })
            else:
                if pltBField:
                    print("**** B field (real)****")
                    self.PlotBField(BField, 0)
                    print("**** B field (imag)****")
                    self.PlotBField(BField, -90)

                    print("**** J field (real)****")
                    self.PlotJphi(JField, 0, 0)
                    print("**** J field (imag)****")
                    self.PlotJphi(JField, -90, 0)
                else:
                    print("**** B field (real)****")
                    Draw (BField.real, mesh, order=self.feOrder, deformation=False) 
                    print("**** B field (imag)****")
                    Draw (BField.imag, mesh, order=self.feOrder, deformation=False) 
                    print("**** J field (real)****")
                    
                    if useModelPlot:
                        self.model.PlotJ(JField, 0)
                        self.model.PlotJ(JField, -90)
                    else:
                        Draw (JField.real, mesh, order=self.feOrder) #,  subdivision=mesh.Materials(model.total_region)) 
                        print("**** J field (imag)****")
                        Draw (JField.imag, mesh, order=self.feOrder) #, subdivision=mesh.Materials(model.total_region)) 


        if plotBFieldonLine:
            #self.PlotBFieldonLine(mesh, BField)
            self.model.PlotBFieldonLine(mesh, BField)
            
    def Solve(self, fes, a,f, tol=1.e-8):
        with TaskManager():
            a.Assemble()
        gf=GridFunction(fes)
        gf=solver.iccg_solve(fes, gf, a, f.vec.FV(), tol=tol, max_iter=1000, accel_factor=0, divfac=10, diviter=10,
                     scaling=True, complex=self.jomega,logplot=True)  
        return gf

    def Refine(self, maxerr, elmErrors, **kwargs):
        default_values = {"adaptive":True,
                          "errorRatio":0.25,
                         }
        default_values.update(kwargs)
        adaptiveRefine=default_values["adaptive"]
        errorRatio=default_values["errorRatio"]
        
        mesh=self.mesh
        if adaptiveRefine:
            error=elmErrors
            elms=0
            for el in mesh.Elements():
                criterion=error[el.nr] > errorRatio*maxerr
                if criterion==True:
                    mesh.SetRefinementFlag(el, True)
                    elms =elms+1
                else:
                    mesh.SetRefinementFlag(el, False)
            print("Number of selected elements to refime mesh =", elms)

        curveOrder=self.model.curveOrder
        mesh.Refine()
        """
        temp_filename = "temp_mesh_after_refine.vol"
        mesh.ngmesh.Save(temp_filename)
        mesh = Mesh(temp_filename)
        os.remove(temp_filename)
        """
        mesh.Curve(curveOrder)
        self.mesh=mesh
        self.model.mesh=mesh
        
        print("Refined mesh: nv=", mesh.nv, " nedge=", mesh.nedge, " nfacet=", mesh.nfacet, " ne=",mesh.ne) 
        Draw(mesh)
 
    def GetB(self, X, Z, BField, phase):
        import numpy as np
        mesh=self.model.mesh
        ex=exp(1j*phase*math.pi/180.)
        N=X.shape[0]
        M=X.shape[1]
        Bx=np.zeros((N, M))
        By=np.zeros((N, M))
        Bz=np.zeros((N, M))
        for n in range(N):
            for m in range(M):
                x=X[n][m]
                z=Z[n][m]
                mip = mesh(x,0,z)
                vx=BField(mip)[0]*ex
                vy=BField(mip)[1]*ex
                vz=BField(mip)[2]*ex
                Bx[n][m]=vx.real
                By[n][m]=vy.real
                Bz[n][m]=vz.real
        return Bx, By, Bz

    def GetJphi(self, X, Z, JField, phase, phi):
        import numpy as np
        mesh=self.model.mesh
        radphi=phi/180.*math.pi
        ex=exp(1j*phase*math.pi/180.)
        N=X.shape[0]
        M=X.shape[1]
        Jphi=np.zeros((N, M))
        for n in range(N):
            for m in range(M):
                x=X[n][m]*cos(radphi)
                y=X[n][m]*sin(radphi)
                z=Z[n][m]
                mip = mesh(x,y,z)
                vx=JField(mip)[0]*ex
                vy=JField(mip)[1]*ex
                Jphi[n][m]=-vx.real*sin(radphi)+vy.real*cos(radphi)
        return Jphi
        
    def PlotBField(self, BField, phase=0):
        import numpy as np
        import matplotlib.pyplot as plt
        a=self.model.radius
        f=self.freq
        ngrid=50

        mesh=self.model.mesh
        r_grid = np.linspace(0, a, ngrid)
        theta_grid = np.linspace(0, np.pi/2, ngrid)
        R, THETA = np.meshgrid(r_grid, theta_grid)
        X = R * np.sin(THETA)
        Z = R * np.cos(THETA)
        Bx, By, Bz=self.GetB(X, Z, BField, phase)
        # 磁場ベクトルのプロット
    
        plt.figure(figsize=(4, 4))
        plt.quiver(X, Z, Bx, Bz,
           np.sqrt(Bx**2 + By**2 + Bz**2), cmap='jet')
        plt.colorbar(label='Magnetic Field Magnitude [T]')
        plt.title(f'Magnetic Field inside the Conducting Sphere (f={f}Hz)')
        plt.xlabel('x [m]')
        plt.ylabel('z [m]')
        plt.axis('equal')
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    def PlotJphi(self, JField, phase=0, phi=0):
        import numpy as np
        import matplotlib.pyplot as plt
        a=0.999*self.model.radius
        f=self.freq
        ngrid=50
        mesh=self.model.mesh
        r_grid = np.linspace(0, a, ngrid)
        theta_grid = np.linspace(0, np.pi/2, ngrid)
        R, THETA = np.meshgrid(r_grid, theta_grid)
        X = R * np.sin(THETA)
        Z = R * np.cos(THETA)
        Jphi=self.GetJphi(X, Z, JField, phase, phi)
        # 磁場ベクトルのプロット
    
        plt.figure(figsize=(4, 4))
        plt.pcolormesh(X, Z, Jphi, cmap='coolwarm')
        plt.colorbar(label='Eddy Current Density J_phi [A/m^2]')
        plt.title(f'Eddy Current Density J_phi inside the Conducting Sphere (f={f}Hz)')
        plt.xlabel('x [m]')
        plt.ylabel('z [m]')
        plt.axis('equal')
        plt.grid(True)
        plt.tight_layout()
        plt.show()



        
   

        
        