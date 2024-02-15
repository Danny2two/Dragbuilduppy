import math
import numpy

class Drag :
    def calc_reynolds(Density: float,V_inf: float,Chord,DynamicViscosity:float)->float:
        """Calculates the Reynolds number for an airfoil at given conditions.
        @Density: Density in kg/m^3
        @V_inf: Velocity of incoming flow in m/s
        @Chord: Chord length in meters, or characteristic length.
        @DynamicViscosity: μ value
        """
        Re: float = (Density * V_inf * Chord)/DynamicViscosity

        return (Re)

    def calc_dynamicviscosity(temperature: float)->float:
        """Calculates dynamic viscosity for a given temperature
        @temperature: Temperature in kelvin.
        """
        𝜇0 = 1.716e-5
        𝑆𝜇 = 113
        T0 = 273.15
        #Dy: float = 𝜇0*((temperature/T0)^(3/2))*((T0+𝑆𝜇)/(temperature+𝑆𝜇))
        Dy: float = 𝜇0*(numpy.power(temperature/T0,1.5))*((T0+𝑆𝜇)/(temperature+𝑆𝜇))
        return(Dy)

    def calc_flatplateskinfriction(reynolds: float, mach: float)->float:
        """Calculates the Flat Plate Skin Friction Coefficient. TURBULENT FLOW ONLY
        @reynolds: reynolds number.
        @mach: mach number.
        """
        Cf = 0.455/(numpy.power((math.log10(reynolds)),(2.58)) * (1 + 0.144*(numpy.power(numpy.power(mach,2),(0.65)))))
        return(Cf)

    def calc_formfactorwing(XC, TC, mach: float, sweepback)->float:
        """Calculates Form Factor for wings, tails, struts and pylons.
        @XC: (x/c)m is the chordwise location of the maximum thickness point.
        @TC: (t/c)m is the thickness ratio. 
        @mach: mach number.
        @sweepback: Sweepback of max thickness line (degrees).
        """
        FF = (1 + (0.6/XC)*TC + 100*numpy.power(TC,4))*((1.34*numpy.power(mach,(0.18)))*numpy.power((numpy.cos(numpy.deg2rad(sweepback))),(0.28)))
        return(FF)
    
    def calc_formfactorfuse(Len,A_max):
        """Calculates Form Factor for Fuselage
        @Len: Length of component
        @A_max: maximum cross-sectional area
        """
        f= Len/numpy.sqrt((4/math.pi)*A_max)
        FF = (0.9 + (5/numpy.power(f,1.5)) + (f/400))
        return(FF)
    
    def calc_formfactornacelle(Len,A_max):
        """Calculates Form Factor for Fuselage
        @Len: Length of component
        @A_max: maximum cross-sectional area
        """
        f= Len/numpy.sqrt((4/math.pi)*A_max)
        FF = 1 + (0.35/f)
        #print(str(f) + " | " + str(FF))
        return(FF)

    def calc_wettedareawing(TC, Area):
        """Estimates wetted area
        @TC: (t/c)m is the thickness ratio.
        @Aera: area of wing
        """
        if (TC > 0.05):
            Swet = (1.977+0.52*TC)*Area
        else:
            Swet = 2.003*Area
        return(Swet)
    
    def calc_oswald_swept(AR, ALe):
        """Calculate the Oswald efficiency factor for a swept wing (ALE >= 30).
        @AR: Aspect Ratio
        @ALe: Angle of sweepback
        """
        e0 = 4.61 * (1-0.045*(numpy.power(AR,0.68))) * numpy.power(numpy.cos(numpy.deg2rad(ALe)),0.15) - 3.1
        return(e0)
    
    def calc_oswald_straight(AR):
        """Calculate the Oswald efficiency factor for a straight wing.
        @AR: Aspect Ratio
        """
        e0 = (1.78*(1-0.045*numpy.power(AR,0.68))) - 0.64
        return(e0)
    
    def calc_oswald_any(AR, ALE):
        """Calculate the Oswald efficiency factor for either swept or straight wings.
        @AR: Aspect Ratio
        @ALe: Angle of sweepback
        """
        if ALE == 0:
            e0 = Drag.calc_oswald_straight(AR)
        elif ALE >=30:
            e0 = Drag.calc_oswald_swept(AR,ALE)
        else:
            e0st = Drag.calc_oswald_straight(AR)
            e0sw = Drag.calc_oswald_swept(AR,30)

            e0 = e0st + ((ALE/30) * (e0sw-e0st))
        return (e0)
    
    def calc_mach(V_inf, Temp):
        """Calculate Mach
        @V_inf: Velocity in m/s
        @Temp: Temperature
        """
        gasconst = 287 #j/kg k
        y = 1.4
        a_inf = math.sqrt(y * gasconst * Temp)
        M = V_inf/a_inf
        return(M)

    def calc_cd0_gear(CD0, S_Frontal, S_wing, Q):
        cd = CD0 * (S_Frontal/S_wing) * Q
        return(cd)

    def calc_cd0(Skin_Friction_Coefficient, Form_Factor, Interference_Factor, Swet, Swing):
        Cd0 = (Skin_Friction_Coefficient * Form_Factor * Interference_Factor * (Swet/Swing))
        return(Cd0)
    