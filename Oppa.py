import math
import matplotlib
import numpy
from Drag import Drag

"""ATMOSPHERE CONDITIONS
"""
Alt = 300 #meters (1000 ft)
Temp = 286.21 #Kelvin
Press = 9.7773 #N/m^2
Dens = 1.1901
Default_Vinf = 45 #m/s

class Wing3d:
    Name = ""
    AR = 0
    OswaldE = 0
    WingSweep = 0
    Area = 0
    Span = 0
    Chord = 1
    A_zero_lift = 0
    A_Stall = 0
    TC = 0.12
    XC = 0.3
    Vinf = 0
    def __init__(self,name,Sweep,Area,Span,Azerolift,Astall,TC,XC) -> None:
        self.Name = name
        AR = (numpy.power(Span, 2))/Area
        OswaldE = Drag.calc_oswald_swept(AR,Sweep)
        WingSweep = Sweep
        Area = Area
        Span = Span
        A_zero_lift = Azerolift
        A_Stall = Astall
        self.TC = TC
        self.XC = XC

if __name__ == "__main__":
    #From CAD model, Area = 208.122 ft^2, 19.335 m^2
    #Span is 15ft, 9.144 m, total length is 4.77m
    NACA4312 = Wing3d("NACA4312",40,19.335,9.114,-6.3,17,0.12,0.3)
    dyn_vc = Drag.calc_dynamicviscosity(Temp)
    rey_nm = Drag.calc_reynolds(Dens,Default_Vinf,NACA4312.Chord,dyn_vc)
    print(NACA4312.Name + "\n Dynamic Viscocity: " + str(dyn_vc) + "\n Reynolds: "+ str(rey_nm))