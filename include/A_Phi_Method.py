from ngsolve import *
from MatrixSolver import MatrixSolver as solver 
def Ht_regularization(H, mesh, boundary, feOrder):
    normal = specialcf.normal(mesh.dim)
    fesu = H1(mesh, order=feOrder, definedon=mesh.Boundaries(boundary), complex=False)
    u, v= fesu.TnT()
    a = BilinearForm(fesu)
    a +=Cross(normal,grad(u).Trace())*Cross(normal, grad(v).Trace())*ds
    f=LinearForm(fesu)
    f +=-H*Cross(normal, grad(v).Trace())*ds
    with TaskManager():
        a.Assemble()
        f.Assemble()
    gfu=GridFunction(fesu)
    gfu=solver.iccg_solve(fesu, gfu, a, f.vec.FV(), tol=1.e-16, max_iter=200, accel_factor=0, 
                          complex=False, logplot=False) 
    hreg=H+Cross(normal,grad(gfu).Trace())
    return hreg