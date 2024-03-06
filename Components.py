import math
import numpy
from Drag import Drag

"""ATMOSPHERE CONDITIONS

Alt = 300 #meters
Temp = 286.21 #Kelvin
Press = 9.77774/4 #N/m^2
Dens = 1.19 #kg/m^3  𝜌
Default_Vinf = 76 #m/s
"""

class Atmosphere:
    Altitude = 0
    Temp = 0
    Pressure = 0
    Density = 0 #𝜌
    Vinfinity = 0
    
    def __init__(self,Altitude, Temp, Pressure, Density, Vinf) -> None:
        """Variables for the atmosphere

        Args:
            Altitude (_type_): Altitude (meters)
            Temp (_type_): Temperature (kelvin)
            Pressure (_type_): Pressure (N/m^2)
            Density (_type_): Desnity of atmosphere (𝜌, kg/m^3)
            Vinf (_type_): Velocity of incoming airstream (m/s)
        """        
        self.Altitude = Altitude
        self.Temp = Temp
        self.Pressure = Pressure
        self.Density = Density
        self.Vinfinity = Vinf
    
    def printAtmos(self):
        print("Atmosphere State:\n ALT: " + str(self.Altitude) +
               " m\n Temp: " + str(self.Temp) + 
               " k\n Pressure: " + str(self.Pressure)+ 
               " N/m^2\n Density: " + str(self.Density)+ 
               " kg/m^3\n V_inf: " + str(self.Vinfinity) + "m/s")

class DragyComponent:
    Name = ""
    Atmos = Atmosphere(0,0,0,0,0)
    def __init__(self,Atmosphere) -> None:
        self.Atmos = Atmosphere
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
    Span = 1.0
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

    def __init__(self,name,Airfoil,Sweep,Area,Span,Chord,Azerolift,Astall,TC,XC, AwingObscured,atmos) -> None:
        """Creates a new 3D wing.

        Args:
            name: Name to be usef for this wing.
            Airfoil: Airfoil used for wing.
            Sweep: Sweep angle in degrees.
            Area: Planform area of wing in square meters.
            Span: Wingspan in meters.
            Chord: Chord length in meters.
            Azerolift: Angle of zero lift in degrees.
            Astall: Stall angle of wing in degrees.
            TC: Thickness Ratio.
            XC: Chordwise location of the maximum thickness.
            AwingObscured: Planform area that is obsucred in square meters.
            atmos (Atmosphere): Atmosphere to be used in calculations
        """
        self.Name = name
        self.airfoil = Airfoil
        self.Chord = Chord
        self.WingSweep = Sweep
        self.AR = (numpy.power(Span, 2))/Area
        self.OswaldE = Drag.calc_oswald_any(self.AR,self.WingSweep)
        self.WingSweep = Sweep
        self.Area = Area
        self.Span = Span
        self.A_zero_lift = Azerolift
        self.A_Stall = Astall
        self.TC = TC
        self.XC = XC
        self.AreaObs = AwingObscured
        super().__init__(Atmosphere=atmos)


    def compute(self):
        self.dyn_vc = Drag.calc_dynamicviscosity(self.Atmos.Temp)
        self.rey_nm = Drag.calc_reynolds(self.Atmos.Density,self.Atmos.Vinfinity,self.Chord,self.dyn_vc)
        self.mach = Drag.calc_mach(self.Atmos.Vinfinity, self.Atmos.Temp)
        self.FPSF = Drag.calc_flatplateskinfriction(self.rey_nm,self.mach)
        self.FormFactor = Drag.calc_formfactorwing(self.XC,self.TC,self.mach,self.WingSweep)
        self.Swet = Drag.calc_wetted_area_wing(self.TC,(self.Area-self.AreaObs))
        self.CD0_wing = Drag.calc_cd0(self.FPSF,self.FormFactor,self.Interf_Factor,self.Swet,self.Area)

    def getCD0(self):
        return(float(self.CD0_wing))

    def printStats(self):
        print("Stats for:" + self.Name )
        print("├Airfoil: " + str(self.airfoil))
        print("├Span: " + str(self.Span)  + " m")
        print("├Area: " + str(self.Area) + " m^2")
        print("├Aspect Ratio: " + str(self.AR))
        print("├Chord: " + str(self.Chord)+ " m")
        print("├Sweep: " + str(self.WingSweep)+ " deg")
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

    def __init__(self,name,len,CrossArea,InterFactor,SWet,SWing,atmos) -> None:
        """Creates a new Nacelle

        Args:
            name (string): Name of nacelle
            len (float): Length in meters
            CrossArea (float): Max cross-sectional area in square meters
            InterFactor (float): Interferance factor
            SWet (float): Wetted surface area in square meters
            SWing (float): Wing area in square meters
            atmos (Atmosphere): Atmosphere to be used in calculations
        """
        self.Name = name
        self.Length = len
        self.CrossSectionArea = CrossArea
        self.InterfFactor = InterFactor
        self.SWet = SWet
        self.SWing =SWing
        super().__init__(Atmosphere=atmos)


    def compute(self):
        self.dyn_vc = Drag.calc_dynamicviscosity(self.Atmos.Temp)
        self.rey_nm = Drag.calc_reynolds(self.Atmos.Density,self.Atmos.Vinfinity,self.Length,self.dyn_vc)
        self.mach = Drag.calc_mach(self.Atmos.Vinfinity, self.Atmos.Temp)
        self.FPSFC = Drag.calc_flatplateskinfriction(self.rey_nm,self.mach)
        self.FF = Drag.calc_formfactornacelle(self.Length,self.CrossSectionArea)
        self.CD0C = Drag.calc_cd0(self.FPSFC,self.FF,self.InterfFactor,self.SWet,self.SWing)

    def printStats(self):
        print("Stats for:" + self.Name )
        print("├Length: " + str(self.Length)+ " m")
        print("├CrossArea: " + str(self.CrossSectionArea)+" m^2")
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
    def __init__(self,name,CDC,SFrontal,SWing,InterfFactor,atmos) -> None:
        """Creates a new fixed gear

        Args:
            name (string): Name for gear.
            CDC (float): Estimated Cd for gear.
            SFrontal (float): Frontal area in square meters.
            SWing (float): Area of main wing in square meters
            InterfFactor (float): Interference factor
        """
        self.Name = name
        self.CDC = CDC
        self.S_Frontal= SFrontal
        self.S_Wing = SWing
        self.interfFactor = InterfFactor
        super().__init__(Atmosphere=atmos)

    def compute(self):
        self.CDC0 = Drag.calc_cd0_gear(self.CDC,self.S_Frontal,self.S_Wing,self.interfFactor)
    
    def printStats(self):
        print("Stats for:" + self.Name )
        print("└DRAG BUILDUP")
        print(" ├Frontal Area: " + str(self.S_Frontal) + " m^2")
        print(" ├ref Wing Area: " + str(self.S_Wing) + " m^2")
        print(" ├CD Comp: " + str(self.CDC))
        print(" ├Interf factor: " + str(self.interfFactor))
        print(" └CD0_gear: " + str(self.CDC0))

    def getCD0(self) -> float:
        return(self.CDC0)
    
class Fuselage(DragyComponent):
    Name = ""
    CDC0 = 0
    length = 1.0
    dyn_vc = 1.0
    FPSF = 1.0
    FormFactor = 1.0
    Swet = 1.0
    atop =1.0
    aside = 1.0
    Interf_Factor =1.0
    maxcross = 1.0
    refwing_area = 1.0
    dyn_vc = 1.0
    rey_nm = 1.0
    mach = 1.0
    def __init__(self,name,length,area_top,area_side,maxcross,interf_factor,refrence_wing_area,atmos) -> None:
        """Creates a new fuselage

        Args:
            name (string): Name for the fuselage.
            length (float): Length of the fuselage in meters.
            area_top (float): Area of the fuselage from top view in square meters.
            area_side (float): Area of the side of the fuslage in square meters.
            maxcross (float): Max cross-section area in square meters.
            interf_factor (float): interference factor.
            refrence_wing_area (float): Area of the main wing as refernce in square meters.
            atmos (Atmosphere): Atmosphere to be used in calculations
        """
        self.Name = name
        self.length = length
        self.atop = area_top
        self.aside = area_side
        self.Interf_Factor = interf_factor
        self.maxcross = maxcross
        self.refwing_area = refrence_wing_area
        super().__init__(Atmosphere=atmos)


    def compute(self):
        self.dyn_vc = Drag.calc_dynamicviscosity(self.Atmos.Temp)
        self.rey_nm = Drag.calc_reynolds(self.Atmos.Density,self.Atmos.Vinfinity,self.length,self.dyn_vc)
        self.mach = Drag.calc_mach(self.Atmos.Vinfinity, self.Atmos.Temp)
        self.FPSF = Drag.calc_flatplateskinfriction(self.rey_nm,self.mach)
        self.FormFactor = Drag.calc_formfactorfuse(self.length,self.maxcross)
        self.Swet = Drag.calc_wetted_area_fuse(self.atop,self.aside)
        self.CDC0 = Drag.calc_cd0(self.FPSF,self.FormFactor,self.Interf_Factor,self.Swet,self.refwing_area)

    
    def printStats(self):
        print("Stats for:" + self.Name )
        print("├Length: " + str(self.length)+ " m")
        print("├Top Area: " + str(self.atop)+" m^2")
        print("├Side Area: " + str(self.aside)+" m^2")
        print("├CrossArea: " + str(self.maxcross)+" m^2")
        print("├Dyn Viscosity: " + str(self.dyn_vc))
        print("├Reynolds: " + str(self.rey_nm))
        print("├Mach: " + str(self.mach))
        print("└DRAG BUILDUP")
        print(" ├Skin Friction: " + str(self.FPSF))
        print(" ├Form Factor: " + str(self.FormFactor))
        print(" ├Wetted Area: " + str(self.Swet) + " m^2")
        print(" ├Interf factor: " + str(self.Interf_Factor))
        print(" └CD0_fuse: " + str(self.CDC0))

    def getCD0(self) -> float:
        return(self.CDC0)

class Engine():
    BSFC = 0
    TSFC = 0
    MaxThrust = 0
    MaxPower =0
    Effic = 0
    Name = ""
    
    def __init__(self,Name,TSFC,BSFC,MAX_Thrust,MAX_Power,Effic) -> None:
        self.Name = Name
        self.TSFC = TSFC #Idealy in g/kN/s
        self.BSFC = BSFC
        self.MaxThrust = MAX_Thrust
        self.MaxPower = MAX_Power
        self.Effic = Effic

    def calc_fuelrate_from_thrust(self,Thrust):
        """_summary_

        Args:
            Thrust (_type_): Required Thrust in Newtons

        Returns:
            _type_: Fuel consumption(g) per second
        """        
        fr = self.TSFC * Thrust/1000
        return fr
    

if __name__ == "__main__":
    pass