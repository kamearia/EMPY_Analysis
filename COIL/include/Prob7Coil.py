from netgen.occ import *
import sys
sys.path.append(r'..\COIL\include')
from EMPY_COIL import *
#sys.path.append(r'C:\EMSolution\EMSolPy3\bin\Release') 
#import EMPY_Field

def Prob7Coil():
    current=2742
    a=12.5e-3
    b=50e-3
    rc=25e-3+a
    w1=50e-3
    w2=w1+rc
    w3=w2+a
    x0=294e-3 - w3
    y0=w3
    z0=19e-3+30e-3+b

    x=x0+w2
    y=y0-w1

    coils=EMPY_COIL()
    block=EMPY_BLOCK( current, gp_Pnt(x0+w2,y0-w1,z0), gp_Pnt(x0+w2,y0+w1,z0), 
                            gp_Vec(0,0,b), gp_Vec(a,0,0) )
    coils.Add(block)
    arc=EMPY_ARC( current, rc, gp_Pnt(x0+w1, y0+w1,z0), 2*a, 2*b,  0, 90)
    coils.Add(arc)
    block=EMPY_BLOCK( current, gp_Pnt(x0+w1,y0+w2,z0), gp_Pnt(x0-w1,y0+w2,z0),
                            gp_Vec(0,0,b), gp_Vec(0,a,0) )
    coils.Add(block)
    arc=EMPY_ARC( current, rc, gp_Pnt(x0-w1, y0+w1,z0), 2*a, 2*b,  90, 180)
    coils.Add(arc)

    block=EMPY_BLOCK( current, gp_Pnt(x0-w2,y0+w1,z0), gp_Pnt(x0-w2,y0-w1,z0),
                            gp_Vec(0,0,b), gp_Vec(-a,0,0) )
    coils.Add(block)
    arc=EMPY_ARC( current, rc, gp_Pnt(x0-w1, y0-w1,z0), 2*a, 2*b,  180, 270)
    coils.Add(arc)

    block=EMPY_BLOCK( current, gp_Pnt(x0-w1,y0-w2,z0), gp_Pnt(x0+w1,y0-w2,z0), 
                            gp_Vec(0,0,b), gp_Vec(0,-a,0) )
    coils.Add(block)
    arc=EMPY_ARC( current, rc, gp_Pnt(x0+w1, y0-w1,z0), 2*a, 2*b,  270, 360)
    coils.Add(arc)

    return coils