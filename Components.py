import math
import numpy
from Atmosphere import Atmosphere
import Drag as DragClass
import pint

class DragyComponent:
    Name = ""
    Atmos: Atmosphere
    def __init__(self,Atmosphere) -> None:
        """DraggyComponents are parts the aircraft that contribute to its CD0 build up.

        Args:
            Atmosphere (_type_): Atmosphere object to be used in calcualtions.
        """        
        self.Atmos = Atmosphere
        self.ur = Atmosphere.ur
        self.Drag = DragClass.Drag(self.ur)

    def compute():
        """Computes the components CD0
        """        
        pass

    def getCD0() -> float:
        """Returns the components CD0

        Returns:
            float: CD0
        """        
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
    CLmax = 0
    CLmin = 0
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
        """Creates a new 3D wing. Units added during creation of wing, do not pass in quantites with demensions

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
        super().__init__(atmos)
        #print(self.Drag)
        self.Name = name
        self.airfoil = Airfoil
        self.Chord = Chord * self.ur.meters
        self.WingSweep = Sweep
        self.AR = (numpy.power(Span, 2))/Area
        self.OswaldE = self.Drag.calc_oswald_any(self.AR,self.WingSweep)
        self.WingSweep = Sweep
        self.Area = Area * self.ur.meters **2
        self.Span = Span * self.ur.meters
        self.A_zero_lift = Azerolift
        self.A_Stall = Astall
        self.TC = TC
        self.XC = XC
        self.AreaObs = AwingObscured * self.ur.meters **2


    def compute(self):
        self.AR = (numpy.power(self.Span, 2))/self.Area
        self.dyn_vc = self.Drag.calc_dynamicviscosity(self.Atmos.Temp)
        self.rey_nm = self.Drag.calc_reynolds(self.Atmos.Density,self.Atmos.Vinfinity,self.Chord,self.dyn_vc)
        self.mach = self.Drag.calc_mach(self.Atmos.Vinfinity, self.Atmos.Temp)
        self.FPSF = self.Drag.calc_flatplateskinfriction(self.rey_nm,self.mach)
        self.FormFactor = self.Drag.calc_formfactorwing(self.XC,self.TC,self.mach,self.WingSweep)
        self.Swet = self.Drag.calc_wetted_area_wing(self.TC,(self.Area-self.AreaObs))
        self.CD0_wing = self.Drag.calc_cd0(self.FPSF,self.FormFactor,self.Interf_Factor,self.Swet,self.Area)

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
        super().__init__(Atmosphere=atmos)
        self.Name = name
        self.Length = len * self.ur.meters
        self.CrossSectionArea = CrossArea * self.ur.meters **2
        self.InterfFactor = InterFactor
        self.SWet = SWet * self.ur.meters **2
        self.SWing =SWing* self.ur.meters **2

    def compute(self):
        self.dyn_vc = self.Drag.calc_dynamicviscosity(self.Atmos.Temp)
        print(f'Dynamic Viscosity: {self.dyn_vc}')
        self.rey_nm = self.Drag.calc_reynolds(self.Atmos.Density,self.Atmos.Vinfinity,self.Length,self.dyn_vc)
        print(f'Reynolds Num: {self.rey_nm}')
        self.mach = self.Drag.calc_mach(self.Atmos.Vinfinity, self.Atmos.Temp)
        print(f'Mach: {self.mach}')
        self.FPSFC = self.Drag.calc_flatplateskinfriction(self.rey_nm,self.mach)
        print(f'FlatPlateSF: {self.FPSFC}')
        self.FF = self.Drag.calc_formfactornacelle(self.Length,self.CrossSectionArea)
        self.CD0C = self.Drag.calc_cd0(self.FPSFC,self.FF,self.InterfFactor,self.SWet,self.SWing)

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
        super().__init__(Atmosphere=atmos)
        self.Name = name
        self.CDC = CDC
        self.S_Frontal= SFrontal * self.ur.meters **2
        self.S_Wing = SWing #* self.ur.meters **2 Already has units. (most of the time :())
        self.interfFactor = InterfFactor

    def compute(self):
        self.CDC0 = self.Drag.calc_cd0_gear(self.CDC,self.S_Frontal,self.S_Wing,self.interfFactor)
    
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
        super().__init__(Atmosphere=atmos)
        self.Name = name
        self.length = length * self.ur.meters
        self.atop = area_top * self.ur.meters**2
        self.aside = area_side * self.ur.meters**2
        self.Interf_Factor = interf_factor
        self.maxcross = maxcross * self.ur.meters**2
        self.refwing_area = refrence_wing_area #* self.ur.meters**2 already has units

    def compute(self):
        self.dyn_vc = self.Drag.calc_dynamicviscosity(self.Atmos.Temp)
        self.rey_nm = self.Drag.calc_reynolds(self.Atmos.Density,self.Atmos.Vinfinity,self.length,self.dyn_vc)
        self.mach = self.Drag.calc_mach(self.Atmos.Vinfinity, self.Atmos.Temp)
        self.FPSF = self.Drag.calc_flatplateskinfriction(self.rey_nm,self.mach)
        self.FormFactor = self.Drag.calc_formfactorfuse(self.length,self.maxcross)
        self.Swet = self.Drag.calc_wetted_area_fuse(self.atop,self.aside)
        self.CDC0 = self.Drag.calc_cd0(self.FPSF,self.FormFactor,self.Interf_Factor,self.Swet,self.refwing_area)

    
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
    
    def __init__(self,Name:str,TSFC,BSFC,MAX_Thrust,MAX_Power,Effic,UnitReg:pint.UnitRegistry) -> None:
        """The engine class is meant to be a Prop or jet engine, However much of the prop/reciprocating engine functionality is not yet implemented.

        Args:
            Name (str): Name of the engine being made
            TSFC (_type_): Thrust specific fuel consumption
            BSFC (_type_): Brake specifif fuel consumption
            MAX_Thrust (_type_): Max thrust (newtons)
            MAX_Power (_type_): Max power (watt)
            Effic (_type_): efficency for props
            UnitReg (pint.UnitRegistry): Unit regestry
        """        
        self.ur: pint.UnitRegistry = UnitReg
        self.Name = Name
        self.TSFC = TSFC #Idealy in g/kN/s
        self.BSFC = BSFC
        self.MaxThrust = MAX_Thrust * self.ur.newton
        self.MaxPower = MAX_Power * self.ur.watt
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
    
class ElectricMotor():
    def __init__(self,Name:str,MaxPower,Effic,UnitReg: pint.UnitRegistry) -> None:
        """PlaceHolder for the future Electric motor class

        Args:
            Name (str): Name
            MaxPower: Max motor power consumption
            Effic: Motor Effic
            UnitReg (pint.UnitRegistry): UnitReg
        """       
        self.ur = UnitReg
        self.Name = Name
        self.MaxPower = MaxPower * self.ur.watt
        self.effic = Effic

    def calc_power_consumption_from_output(self,Output_power):
        return (Output_power / self.effic)
    
class Battery():
    def __init__(self,EnergyDensity, mass,voltage,UnitReg: pint.UnitRegistry,currentEnergy = -1) -> None:
        self.ur = UnitReg
        self.density = EnergyDensity
        self.mass = mass* self.ur.kilogram
        self.voltage = voltage *self.ur.volt
        self.maxCapacity = EnergyDensity * mass * self.ur.watt * self.ur.hours
        self.MaxEnergy = self.maxCapacity * (3600 * self.ur.joules / (self.ur.watt * self.ur.hours))
        if currentEnergy != -1:
            self.CurrentEnergy = currentEnergy * self.ur.joule
        else:
            self.CurrentEnergy = self.MaxEnergy

    def discharge(self,current, voltage, time):
        self.CurrentEnergy -= (voltage * current * time)
        if self.CurrentEnergy < 0:
            print("Battery over drawn! We are now breaking the laws of physics!")
        else:
            print(f'Used {voltage * current * time} of energy. \n Current remaing energy: {self.CurrentEnergy}')

    def print_state(self):
        print("Battery Status")
        print(f' Current Energy stored {self.CurrentEnergy}. \n Max Energy {self.MaxEnergy}. \n Voltage {self.voltage}. \n State of Charge {self.CurrentEnergy/self.MaxEnergy}')

    def get_energy(self):
        return self.CurrentEnergy
    

if __name__ == "__main__":
    pass