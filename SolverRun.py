# -*- coding: utf-8 -*-
from ngsolve import *
from ngsolve.webgui import Draw
from netgen.webgui import Draw as DrawGeo 
import math
import sys
#sys.path.append('/home/jupyter-ksugahar/‹T—L/include')
sys.path.append('include')
from MatrixSolver import MatrixSolver as solver
sys.path.append('C:\EMSolution\EMSolpy4\python\model')
from BathPlateModel import BathPlateModel
from BathPlateModel import MeasureFace
sys.path.append(r'..\bin\\Release') 
#sys.path.append(r'..\..\..\x64\Release') 
from EMPY_Field import *


curveOrder=1
feOrder=1
f=50
jomega=False
if jomega==True:
    s=2j*math.pi*f
else:
    s=2*math.pi*f   
Dirichlet=False
boxz=200e-3
B0=1 #(Tesla)
boxx=200e-3
boxy=200e-3
model=BathPlateModel(holes=2,boxx=boxx, boxy=boxy, boxz=boxz, div_thick=1)

DrawGeo(model.GetGeometry())
mesh=model.GetMesh()              
Draw(mesh)

if Dirichlet==False:
    fesA = HCurl(mesh, order=feOrder, nograds=True, complex=jomega)
else:
    fesA = HCurl(mesh, order=feOrder, nograds=True, dirichlet=model.outer_boundary, complex=jomega)
fesPhi = H1(mesh, order=feOrder, definedon=model.conductive_region, complex=jomega)
fesAPhi=fesA*fesPhi
(A,phi),(N,psi) = fesAPhi.TnT()


a = BilinearForm(fesAPhi)
a += s*model.Sigma*(A+grad(phi))*(N+grad(psi))*dx(model.conductive_region)
a += 1/model.Mu*curl(A)*curl(N)*dx
with TaskManager():
    a.Assemble()
gfAPhi = GridFunction(fesAPhi)
gfA, gfPhi = gfAPhi.components

#field=UNIF(0,0,B0,0)
if Dirichlet==False:
    normal = specialcf.normal(mesh.dim)
    mu=4.e-7*math.pi
    #h=EMPY_Field.UNIF(0.,0.,1.,0).B()
    #mip = mesh(0.,0.,0.)
    #print("h(mip)=",h(mip))
    h=(0.,0.,1.)
    sr = LinearForm(fesAPhi)
    sr += 1./mu*Cross(N.Trace(),h)*normal*ds(model.outer_boundary)
    #sr += Cross(N.Trace(),h)*normal*ds('upper|lower|xmax|xmin|ymax|ymin')
    print("point1")
    #with TaskManager():
    sr.Assemble()
    print(sr)
    print("point2")
    gfAPhi=solver.iccg_solve(fesAPhi, gfAPhi, a, sr.vec.FV(), tol=1.e-16, max_iter=200, accel_factor=1.1, complex=jomega)
else:
    Av=EMPY_Field.UNIF(0.,0.,1.,0).A()
    #gfA.Set(Av, BND, mesh.Boundaries('upper|lower|xmax|xmin|ymax|ymin'))
    gfA.Set(Av, BND, mesh.Boundaries(model.outer_boundary))
    print("point11")
    r=-a.mat*gfA.vec
    #print("point12")
    gfAPhi=solver.iccg_solve(fesAPhi, gfAPhi, a, r.Evaluate(), tol=1.e-16, max_iter=200, accel_factor=1.1, complex=jomega)

Bfield=curl(gfA)
Jfield=-model.Sigma*s*(gfA+grad(gfPhi))

from BathPlateModel import MeasureFace
flux=MeasureFace.CalcFluxes(Jfield)
print("Current fluxes=", flux)

Draw (Bfield, mesh, order=3, min=0.5, max=1.5, deformation=False)
Draw (Jfield, mesh, order=3, min=0, max=2.e8, deformation=False)
