import math
import numpy
from Drag import Drag

"""ATMOSPHERE CONDITIONS
"""
Alt = 3000 #meters (10000 ft)
Temp = 268 #Kelvin
Press = 70121 #N/m^2
Dens = 0.905 #kg/m^3  𝜌
Default_Vinf = 92.59 #m/s

class DragyComponent:
    def __init__(self) -> None:
        pass

    def compute():
        pass

    def getCD0() -> float:
        pass

    def printStats():
        pass

class Wing3d(DragyComponent):
    Name = ""
    AR = .0
    OswaldE = .0
    WingSweep = .0
    Area = .0
    AreaObs = .0
    Span = .0
    Chord = .1
    A_zero_lift = .0
    A_Stall = .0
    TC = 0.12
    XC = 0.3
    Vinf = .0
    Interf_Factor = 1
    ##To be computed later
    dyn_vc = .0
    rey_nm = .0
    mach = .0
    FPSF = .0
    FormFactor = .0
    Swet = .0
    CD0_wing = .0
    airfoil = ""

    def __init__(self,name,Airfoil,Sweep,Area,Span,Chord,Azerolift,Astall,TC,XC, AwingObscured) -> None:
        self.Name = name
        self.airfoil = Airfoil
        self.Chord = Chord
        self.AR = (numpy.power(Span, 2))/Area
        self.OswaldE = Drag.calc_oswald_any(AR=self.AR,ALE=Sweep)
        self.WingSweep = Sweep
        self.Area = Area
        self.Span = Span
        self.A_zero_lift = Azerolift
        self.A_Stall = Astall
        self.TC = TC
        self.XC = XC
        self.AreaObs = AwingObscured

    def compute(self):
        self.dyn_vc = Drag.calc_dynamicviscosity(Temp)
        self.rey_nm = Drag.calc_reynolds(Dens,Default_Vinf,self.Chord,self.dyn_vc)
        self.mach = Drag.calc_mach(Default_Vinf, Temp)
        self.FPSF = Drag.calc_flatplateskinfriction(self.rey_nm,self.mach)
        self.FormFactor = Drag.calc_formfactorwing(self.XC,self.TC,self.mach,self.WingSweep)
        self.Swet = Drag.calc_wetted_area_wing(self.TC,(self.Area-self.AreaObs))
        self.CD0_wing = Drag.calc_cd0(self.FPSF,self.FormFactor,self.Interf_Factor,self.Swet,self.Area)

    def getCD0(self):
        return(float(self.CD0_wing))

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


class Nacelle(DragyComponent):
    Name = ""
    Length = .0
    CrossSectionArea = .0
    InterfFactor = 1
    CD0C = .0
    FPSFC = .0
    FF = .0
    SWet = .0
    SWing = .0

    def __init__(self,name,len,CrossArea,InterFactor,SWet,SWing) -> None:
        self.Name = name
        self.Length = len
        self.CrossSectionArea = CrossArea
        self.InterfFactor = InterFactor
        self.SWet = SWet
        self.SWing =SWing

    def compute(self):
        self.dyn_vc = Drag.calc_dynamicviscosity(Temp)
        self.rey_nm = Drag.calc_reynolds(Dens,Default_Vinf,self.Length,self.dyn_vc)
        self.mach = Drag.calc_mach(Default_Vinf, Temp)
        self.FPSFC = Drag.calc_flatplateskinfriction(self.rey_nm,self.mach)
        self.FF = Drag.calc_formfactornacelle(self.Length,self.CrossSectionArea)
        self.CD0C = Drag.calc_cd0(self.FPSFC,self.FF,self.InterfFactor,self.SWet,self.SWing)

    def printStats(self):
        print("Stats for:" + self.Name )
        print("├Length: " + str(self.Length))
        print("├CrossArea: " + str(self.CrossSectionArea))
        print("├Dyn Viscosity: " + str(self.dyn_vc))
        print("├Reynolds: " + str(self.rey_nm))
        print("├Mach: " + str(self.mach))
        print("└DRAG BUILDUP")
        print(" ├Skin Friction: " + str(self.FPSFC))
        print(" ├Form Factor: " + str(self.FF))
        print(" ├Wetted Area: " + str(self.SWet))
        print(" ├Interf factor: " + str(self.InterfFactor))
        print(" └CD0_nac: " + str(self.CD0C))

    def getCD0(self):
        return(float(self.CD0C))
    
class FixedGear(DragyComponent):
    Name = ""
    CDC = 0
    S_Frontal = 0
    S_Wing = 0
    interfFactor = 0
    CDC0 = 0
    def __init__(self,name,CDC,SFrontal,SWing,InterfFactor) -> None:
        self.Name = name
        self.CDC = CDC
        self.S_Frontal= SFrontal
        self.S_Wing = SWing
        self.interfFactor = InterfFactor
        super().__init__()

    def compute(self):
        self.CDC0 = Drag.calc_cd0_gear(self.CDC,self.S_Frontal,self.S_Wing,self.interfFactor)
    
    def printStats(self):
        print("Stats for:" + self.Name )
        print("└DRAG BUILDUP")
        print(" ├Frontal Area: " + str(self.S_Frontal))
        print(" ├ref Wing Area: " + str(self.S_Wing))
        print(" ├CD Comp: " + str(self.CDC))
        print(" ├Interf factor: " + str(self.interfFactor))
        print(" └CD0_gear: " + str(self.CDC0))

    def getCD0(self) -> float:
        return(self.CDC0)
    
class Fuselage(DragyComponent):
    Name = ""

    CDC0 = 0
    def __init__(self,name,CD0) -> None:
        self.Name = name
        self.CDC0 = CD0
        super().__init__()

    def compute(self):
        pass
    
    def printStats(self):
        print("Stats for:" + self.Name )
        print("└DRAG BUILDUP")
        print(" └CD0_fuse: " + str(self.CDC0))

    def getCD0(self) -> float:
        return(self.CDC0)


if __name__ == "__main__":
    #Components of DC-3
    MainFuselage = Fuselage("DC3 Fuselage",0.0046) #Not Computing, just placeholding
    MainWing = Wing3d("MainWing","NACA2215",10.5,91.7,29,3.16,-1.85,14,0.15,0.3,(4.32*2.4384)) #Specs of Wing
    HorizTail = Wing3d("HorizTail","NA",25,19.5096,1,2.578608,0,0,0.08,0.5,0) #Specs of the horizontal tail, ALL UNITS CONVERTED TO METRIC, meters or m^2
    VertTail = Wing3d("VertTail","NA",33,9.0116,1,3.56616,0,0,0.09,0.50,0) #Specs of the Vertical tail, ALL UNITS CONVERTED TO METRIC, meters or m^2
    LNacelle = Nacelle("LeftNacelle",4.84632,1.5886,1.5,15.0967,91.7) #Specs of nacelles (metric)
    RNacelle = Nacelle("RightNacelle",4.84632,1.5886,1.5,15.0967,91.7)
    TailGear = FixedGear("Tail Gear",0.25,0.196129,91.7,1.2) #Specs of gear (metric)

    HorizTail.Interf_Factor = 1.04
    VertTail.Interf_Factor = 1.04

    #Compute CD0's for all the components
    MainWing.compute()
    HorizTail.compute()
    VertTail.compute()
    LNacelle.compute()
    RNacelle.compute()
    TailGear.compute()

    #Print all of the computed Values
    MainFuselage.printStats()
    MainWing.printStats()
    HorizTail.printStats()
    VertTail.printStats()
    LNacelle.printStats()
    RNacelle.printStats()
    TailGear.printStats()

    CD0ALL = MainFuselage.getCD0() + MainWing.getCD0() + HorizTail.getCD0() + VertTail.getCD0() + LNacelle.getCD0() + RNacelle.getCD0() + TailGear.getCD0()
    print("Overall CD0: " + str(CD0ALL))

    """
        Some Info for Power requred calculation
    """
    K = 0.05
    mass = 11430 
    W = 112128.3

    M1 = 0.5 * (Dens * math.pow(Default_Vinf,2) * MainWing.Area * CD0ALL)
    M2 = (2 * K * math.pow(W,2))/(Dens*math.pow(Default_Vinf,2)*MainWing.Area)
    ThrustRequred = M1 + M2

    H1 = 0.5*(Dens * (math.pow(Default_Vinf, 3)) * MainWing.Area * CD0ALL)
    H2 = (2 * K * math.pow(W,2))/ (Dens * Default_Vinf * MainWing.Area)
    PowerRequred = H1+H2

    print("Thrust Requred for level flight:" + str(ThrustRequred) + " N")
    print("Power Requred for level flight:" + str(PowerRequred) + " Watts")