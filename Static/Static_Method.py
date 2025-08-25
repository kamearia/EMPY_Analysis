from ngsolve import *
import sys
sys.path.append('..\include')
from MatrixSolver import MatrixSolver as solver 

class Static_Method():
    def __init__(self,  model, coil,  **kwargs):  
        self.Calc(model, coil,  **kwargs)

    #def Calc(self, model, coil,  **kwargs):
    #    return

    def CalcResult(self, model, BField):
        mesh=model.mesh
        mip = mesh(0,0,0)
        print("center magnetic field = ", BField(mip),"  ")
        Wm=Integrate(BField*BField/self.Mu*dx("iron"), mesh)
        print("magnetic energy=", Wm,"  ")

        from ngsolve.webgui import Draw
        #print("**** Omega field ****")
        #Draw (gfOmega, mesh, order=feOrder, deformation=False) 
        print("**** B field ****")
        Draw (BField, mesh, order=self.feOrder, min=0., max=5.0, deformation=False) 

    def Solve(self,fes, a,f):
        with TaskManager():
            a.Assemble()
        gf=GridFunction(fes)
        gf=solver.iccg_solve(fes, gf, a, f.vec.FV(), tol=1.e-8, max_iter=1000, accel_factor=0, divfac=100, diviter=10,
                     scaling=True, complex=False,logplot=True)  
        return gf