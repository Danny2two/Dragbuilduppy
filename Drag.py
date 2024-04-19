import math
import numpy

import pint

class Drag :
    def __init__(self,ur: pint.UnitRegistry) -> None:
        """The Drag class serves as a hub for calculations around Cd0

        Args:
            ur (pint.UnitRegistry): Unitregestry, Should be a refrence to the Crafts Unit reg
        """        
        self.ur = ur

    def calc_reynolds(self,Density: float,V_inf: float,Chord,DynamicViscosity:float)->float:
        """Calculates the Reynolds number for an airfoil at given conditions.
        @Density: Density in kg/m^3
        @V_inf: Velocity of incoming flow in m/s
        @Chord: Chord length in meters, or characteristic length.
        @DynamicViscosity: μ value
        """
        #print(f'Density: {Density}')
        #print(f'Velocity: {V_inf}')
        #print(f'Chord: {Chord}')
        #print(f'Viscosity: {DynamicViscosity}')
        Re: float = (Density * V_inf * Chord)/DynamicViscosity
        #print(f'RE {Re}')
        return (Re.magnitude)

    def calc_dynamicviscosity(self,temperature: float)->float:
        """Calculates dynamic viscosity for a given temperature
        @temperature: Temperature in kelvin.
        """
        𝜇0 = 1.716e-5 #Constant
        𝑆𝜇 = 113 * self.ur.kelvin# Constant
        T0 = 273.15  * self.ur.kelvin
        #Dy: float = 𝜇0*((temperature/T0)^(3/2))*((T0+𝑆𝜇)/(temperature+𝑆𝜇))
        Dy: float = 𝜇0*(numpy.power(temperature/T0,1.5))*((T0+𝑆𝜇)/(temperature+𝑆𝜇))
        return(Dy)

    def calc_flatplateskinfriction(self,reynolds: float, mach: float)->float:
        """Calculates the Flat Plate Skin Friction Coefficient. TURBULENT FLOW ONLY
        @reynolds: reynolds number.
        @mach: mach number.
        """
        #print(f'MACH: {mach}')
        #print(f'Reynolds: {reynolds}')
        Cf = 0.455/(numpy.power((numpy.log10(reynolds)),(2.58)) * numpy.power(1 + 0.144*(numpy.power(mach,2)),(0.65)))
        return(Cf)

    def calc_formfactorwing(self,XC, TC, mach: float, sweepback)->float:
        """Calculates Form Factor for wings, tails, struts and pylons.
        @XC: (x/c)m is the chordwise location of the maximum thickness point.
        @TC: (t/c)m is the thickness ratio. 
        @mach: mach number.
        @sweepback: Sweepback of max thickness line (degrees).
        """
        FF = (1 + (0.6/XC)*TC + 100*numpy.power(TC,4))*((1.34*numpy.power(mach,(0.18)))*numpy.power((numpy.cos(numpy.deg2rad(sweepback))),(0.28)))
        return(FF)
    
    def calc_formfactorfuse(self,Len,A_max):
        """Calculates Form Factor for Fuselage
        @Len: Length of component
        @A_max: maximum cross-sectional area
        """
        f= Len/numpy.sqrt((4/math.pi)*A_max)
        FF = (0.9 + (5/numpy.power(f,1.5)) + (f/400))
        #print(f'Formfactor: {FF}')
        return(FF)
    
    def calc_formfactornacelle(self,Len,A_max):
        """Calculates Form Factor for Fuselage
        @Len: Length of component
        @A_max: maximum cross-sectional area
        """
        f= Len/numpy.sqrt((4/math.pi)*A_max)
        FF = 1 + (0.35/f)
        #print(str(f) + " | " + str(FF))
        return(FF)

    def calc_wetted_area_wing(self,TC, Area):
        """Estimates wetted area
        @TC: (t/c)m is the thickness ratio.
        @Aera: area of wing
        """
        if (TC > 0.05):
            Swet = (1.977+0.52*TC)*Area
        else:
            Swet = 2.003*Area
        return(Swet)
    
    def calc_wetted_area_fuse(self,Atop,Aside):
        """Estimates wetted area of fuselage
        @Atop: area of top view of fuse
        @Aside: side area
        """
        swet = 3.4*((Atop + Aside)/2)
        #print(f'Fuselage Swet: {swet}')
        return (swet)
    
    def calc_oswald_swept(self,AR, ALe):
        """Calculate the Oswald efficiency factor for a swept wing (ALE >= 30).
        @AR: Aspect Ratio
        @ALe: Angle of sweepback
        """
        e0 = 4.61 * (1-0.045*(numpy.power(AR,0.68))) * numpy.power(numpy.cos(numpy.deg2rad(ALe)),0.15) - 3.1
        return(e0)
    
    def calc_oswald_straight(self,AR):
        """Calculate the Oswald efficiency factor for a straight wing.
        @AR: Aspect Ratio
        """
        e0 = (1.78*(1-0.045*numpy.power(AR,0.68))) - 0.64
        return(e0)
    
    def calc_oswald_any(self,AR, ALE):
        """Calculate the Oswald efficiency factor for either swept or straight wings.
        @AR: Aspect Ratio
        @ALe: Angle of sweepback
        """
        if ALE == 0:
            e0 = self.calc_oswald_straight(AR)
        elif ALE >=30:
            e0 = self.calc_oswald_swept(AR,ALE)
        else:
            e0st = self.calc_oswald_straight(AR)
            e0sw = self.calc_oswald_swept(AR,30)

            e0 = e0st + ((ALE/30) * (e0sw-e0st))
        return (e0)
    
    def calc_mach(self,V_inf, Temp):
        """Calculates Mach

        Args:
            V_inf (_type_): Velocity
            Temp (_type_): Temperature (K)
        """
        ur = self.ur
        gasconst = 287 * ur.m**2 / (ur.s**2 * ur.kelvin)
        y = 1.4  #ratio
        #print(f'TEMP {Temp}')
        #print(f'Vel {V_inf}')

        M =V_inf/ numpy.sqrt(y * (gasconst * Temp))
        #print(f'Mach: {M}')
        return(M)

    def calc_cd0_gear(self,CD0, S_Frontal, S_wing, Q):
        """Calcualts Cd0 of landing gear

        Args:
            CD0 (_type_): _description_
            S_Frontal (_type_): Frontal surface area 
            S_wing (_type_): main wing surface area
            Q (_type_): Interf factor
        """        
        cd = CD0 * (S_Frontal/S_wing) * Q
        return(cd)

    def calc_cd0(self,Skin_Friction_Coefficient, Form_Factor, Interference_Factor, Swet, Swing):
        """Calculates Cd0 for components that have standard CDO buildup

        Args:
            Skin_Friction_Coefficient (_type_): Flatplate skin friction coeff
            Form_Factor (_type_): Formfactor 
            Interference_Factor (_type_): Interference factor 
            Swet (_type_): Wetted area
            Swing (_type_): Mainwing area
        """        
        Cd0 = (Skin_Friction_Coefficient * Form_Factor * Interference_Factor * (Swet/Swing))
        #print(f'CD0 Units: FPSF:{Skin_Friction_Coefficient}:, Formfactor: {Form_Factor}: InterfFactor: {Interference_Factor}: Swet: {Swet}: Swing: {Swing}')
        return(Cd0)
    
