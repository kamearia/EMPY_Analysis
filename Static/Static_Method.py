from ngsolve import *
import sys
sys.path.append('..\include')
from MatrixSolver import MatrixSolver as solver 

class Static_Method():
    def __init__(self,  model, coil,  **kwargs):  
        self.model=model
        self.Calc(model, coil,  **kwargs)

    #def Calc(self, model, coil,  **kwargs):
    #    return

    def CalcResult(self, model, BField):
        mesh=model.mesh
        mip = mesh(0,0,0)
        print("center magnetic field = ", BField(mip),"  ")
        Wm=Integrate(BField*BField/self.Mu*dx("iron"), mesh)/2
        print("magnetic energy=", Wm,"  ")
        
        from ngsolve.webgui import Draw
        #print("**** Omega field ****")
        #Draw (gfOmega, mesh, order=feOrder, deformation=False) 
        print("**** B field ****")
        Draw (BField, mesh, order=self.feOrder, min=0., max=5.0, deformation=False) 

    def Solve(self, fes, a,f):
        with TaskManager():
            a.Assemble()
        gf=GridFunction(fes)
        gf=solver.iccg_solve(fes, gf, a, f.vec.FV(), tol=1.e-8, max_iter=1000, accel_factor=0, divfac=10, diviter=10,
                     scaling=True, complex=False,logplot=True)  
        return gf

    def Refine(self, maxerr, elmErrors):
        mesh=self.mesh
        error=elmErrors
        elms=0
        for el in mesh.Elements():
            criterion=error[el.nr] > 0.25*maxerr
            if criterion==True:
                mesh.SetRefinementFlag(el, True)
                elms =elms+1
            else:
                mesh.SetRefinementFlag(el, False)
        print("Number of selected elements to refime mesh =", elms)
        curveOrder=self.model.curveOrder
        mesh.Refine()
        mesh.Curve(curveOrder)
        print("Refined mesh: nv=", mesh.nv, " nedge=", mesh.nedge, " nfacet=", mesh.nfacet, " ne=",mesh.ne) 