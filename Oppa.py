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
    Interf_Factor = 1
    ##To be computed later
    dyn_vc = 0
    rey_nm = 0
    mach = 0
    FPSF = 0
    FormFactor = 0
    Swet = 0
    CD0_wing = 0
    airfoil = ""

    def __init__(self,name,Airfoil,Sweep,Area,Span,Chord,Azerolift,Astall,TC,XC) -> None:
        self.Name = name
        self.airfoil = Airfoil
        self.Chord = Chord
        self.AR = (numpy.power(Span, 2))/Area
        self.OswaldE = Drag.calc_oswald_swept(self.AR,Sweep)
        self.WingSweep = Sweep
        self.Area = Area
        self.Span = Span
        self.A_zero_lift = Azerolift
        self.A_Stall = Astall
        self.TC = TC
        self.XC = XC

    def compute(self):
        self.dyn_vc = Drag.calc_dynamicviscosity(Temp)
        self.rey_nm = Drag.calc_reynolds(Dens,Default_Vinf,NACA4312.Chord,self.dyn_vc)
        self.mach = Drag.calc_mach(Default_Vinf, Temp)
        self.FPSF = Drag.calc_flatplateskinfriction(self.rey_nm,self.mach)
        self.FormFactor = Drag.calc_formfactorwing(NACA4312.XC,NACA4312.TC,self.mach,NACA4312.WingSweep)
        self.Swet = Drag.calc_wettedareawing(NACA4312.TC,NACA4312.Area)
        self.CD0_wing = Drag.calc_cd0(self.FPSF,self.FormFactor,NACA4312.Interf_Factor,self.Swet,NACA4312.Area)

    def printStats(self):
        print("Stats for:" + self.Name )
        print("├Airfoil: " + str(self.airfoil))
        print("├Span: " + str(self.Span))
        print("├Area: " + str(self.Area))
        print("├Aspect Ratio: " + str(self.AR))
        print("├Chord: " + str(self.Chord))
        print("├Sweep: " + str(self.WingSweep))
        print("├Dyn Viscosity: " + str(self.dyn_vc))
        print("├Reynolds: " + str(self.rey_nm))
        print("├Mach: " + str(self.mach))
        print("├Oswald Eff: " + str(self.OswaldE))
        print("└DRAG BUILDUP")
        print(" ├Skin Friction: " + str(self.FPSF))
        print(" ├Form Factor: " + str(self.FormFactor))
        print(" ├Wetted Area: " + str(self.Swet))
        print(" ├Interf factor: " + str(self.Interf_Factor))
        print(" └CD0_Wing: " + str(self.CD0_wing))

if __name__ == "__main__":
    #From CAD model, wing Area = 208.122 ft^2, 19.335 m^2
    #Span is 15ft, 9.144 m, total length is 4.77m, ave chord 2.36 m
    NACA4312 = Wing3d("OppaStoppa","NACA4312",40,19.335,9.114,2.36,-4,17,0.12,0.3)
    NACA4312.compute()
    NACA4312.printStats()