from netgen.occ import *
from ngsolve import *
import math
import sys
sys.path.append(r'..\bin\Release') 
import EMPY_Field

def toEM_Pnt(ng):
    return EMPY_Field.EM_3dPoint(ng.x, ng.y, ng.z)
def toEM_Vec(ng):
    return EMPY_Field.EM_3dVector(ng.x, ng.y, ng.z)

def DegtoRad(ang): return ang*math.pi/180

#***** BLOCK ******
class EMPY_BLOCK():
    def __init__(self, current, start_point, end_point, dx, dy):
        self.field=EMPY_Field.BLOCK(current, toEM_Pnt(start_point), toEM_Pnt(end_point),
                              toEM_Vec(dx), toEM_Vec(dy) )
        p1=start_point-dx-dy
        p2=end_point+dx+dy
        self.geo=Box(p1, p2)

class EMPY_LOOP():
    def __init__(self, current, radius, pz, radialWidth, axialWidth):
        if isinstance(pz, gp_Pnt):
            center=gp_Pnt(pz.x, pz.y, pz.z)
        else:
            center=gp_Pnt(0,0,pz)
            
        self.field=EMPY_Field.LOOP(current, radius, toEM_Pnt(center), radialWidth, axialWidth )
        a=radialWidth
        b=axialWidth
        angle=360
        f = WorkPlane(Axes((radius-a/2,0,-b/2), n=-Y, h=X)).Rectangle(a,b).Face()
        self.geo = f.Revolve(Axis((0,0,0),Z),angle)
        self.geo =self.geo.Move(gp_Vec(center.x, center.y, center.z))

class EMPY_ARC():
    def __init__(self, current, radius, pz, radialWidth, axialWidth, start_angle, end_angle):
        if isinstance(pz, gp_Pnt):
            center=gp_Pnt(pz.x, pz.y, pz.z)
        else:
            center=gp_Pnt(0,0,pz)

        self.field=EMPY_Field.ARC(current, radius, toEM_Pnt(center), radialWidth, axialWidth,
                                  DegtoRad(start_angle), DegtoRad(end_angle) )   
        a=radialWidth
        b=axialWidth
        f = WorkPlane(Axes((radius-a/2,0,-b/2), n=-Y, h=X)).Rectangle(a,b).Face()
        self.geo = f.Revolve(Axis((0,0,0),Z),end_angle-start_angle)
        self.geo =self.geo.Rotate(Axis((0,0,0),Z), start_angle)
        self.geo =self.geo.Move(gp_Vec(center.x, center.y, center.z))

class EMPY_COIL():
    def __init__(self):
        self.coils=[]
        self.field=EMPY_Field.EM_COIL()

    def Add(self, coil):
        self.coils.append(coil)
        self.field.Add(coil.field)
        if len(self.coils)==1: self.geo=coil.geo
        else : self.geo= self.geo+coil.geo

    

        