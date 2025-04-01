from netgen.occ import *
import sys
sys.path.append(r'..\COIL\include')
from EMPY_COIL import *

def Prob3Coil():
    rin=20e-3
    rout=40e-3
    h=20e-3
    thick=6.35e-3
    z=thick/2+15e-3
    x=20e-3
    radius=(rin+rout)/2
    radialWidth=rout-rin
    axialWidth=h
    zp=z+h/2
    center=gp_Pnt(20e-3,0, zp)
    J=1260  #(AT)  
    loop=EMPY_LOOP(J, radius, center, radialWidth, axialWidth)
    return loop