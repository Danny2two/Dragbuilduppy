from craft import Craft
from Components import *
import StFl
import Atmosphere
import Drag
import CraftStatistics
import LoadFactor

class Takeoff:
    def __init__(self, Craft:Craft,Alt,n_takeoff,climbangle) -> None:
        self.Craft = Craft
        self.CraftStats = CraftStatistics.CraftStatistics(Craft)
        self.LoadFactor = LoadFactor.LoadFactor(Craft)
        self.ur = Craft.ur
        self.ActiveAtmospere = Craft.Atmosphere
        self.Alt = Alt * ur.m
        self.g = 9.81 * ur.m /ur.second**2
        self.RollFricCoe = 0.04
        self.weight = Craft.weight_takeoff
        self.V_LO = 1.2 * StFl.calc_Vstall(self.ActiveAtmospere.dens_trop_alt(Alt),self.weight,Craft.mainwing.Area,Craft.CLmax).to("meter/second")
        self.n_to = n_takeoff
        self.height_Obs = 10.668 * ur.m #35 ft
        self.climbAngle = climbangle
    
    def GroundRoll(self):
        weight = self.weight
        wingArea = self.Craft.mainwing.Area
        dens = self.ActiveAtmospere.dens_trop_alt(self.Alt)
        thrust = self.CraftStats.get_ThrustAvailable_jet(self.Alt)
        Lift = StFl.calc_Lift(self.Craft.CLrolling,0.7 * self.V_LO,wingArea,dens).to("newton")
        Drag = StFl.calc_Drag(self.Craft.Cd0,0.7 * self.V_LO,wingArea,dens,self.Craft.K,self.Craft.CLrolling)

        RollDist =  (1.44 * numpy.power(weight,2))/(self.g * dens * wingArea * self.Craft.CLmax * (thrust - Drag - (self.RollFricCoe*(weight - Lift))))
        #print((self.g * dens * wingArea * self.Craft.CLmax * (thrust - Drag - (self.RollFricCoe*(weight - Lift)))))
        #print(f'GR: {RollDist}, w: {weight}, wa: {wingArea}, Dens: {dens}, thr: {thrust}, L: {Lift}, D: {Drag}, COE: {self.RollFricCoe}, g: {self.g}')
        print(f'Rolling Distance of {RollDist.to("meter")}')
        return (RollDist)
    
    def Transition(self):
        radius = self.LoadFactor.getRadiusPullUp(self.n_to,self.V_LO)
        Str = radius * numpy.sin(numpy.deg2rad(self.climbAngle))
        htr = radius - (radius * numpy.cos(numpy.deg2rad(self.climbAngle)))
        #print(Str,htr)
        print(f'Tranition distance of {Str.to("meter")} at VLO: {self.V_LO}')
        return([Str,htr])
    
    def AirDist(self,remainingObsH):
        Sa = remainingObsH / numpy.tan(numpy.deg2rad(self.climbAngle))
        #print(Sa)
        print(f'Air Distance of {Sa.to("meter")}')
        return Sa
    
    def DoTakeoff(self):
        ToD = 0
        remainingH = self.height_Obs

        ToD += self.GroundRoll()

        trans = self.Transition()
        ToD += trans[0]
        remainingH -= trans[1]

        ToD += self.AirDist(remainingH)
        return(ToD)
    
class Landing:
    def __init__(self, Craft: Craft, Alt,n_flare,ApproachAngle,weight) -> None:
        self.Craft = Craft
        self.CraftStats = CraftStatistics.CraftStatistics(Craft)
        self.LoadFactor = LoadFactor.LoadFactor(Craft)
        self.ur = Craft.ur
        self.Atmos = Craft.Atmosphere
        self.weight=weight
        self.Alt = Alt * ur.m
        self.n_flare = n_flare
        self.apprachA = ApproachAngle
        self.g = 9.81 * ur.m /ur.second**2
        self.RollFricCoe = 0.4
        self.obsH = 0.668 * ur.m #35 ft
        self.V_A = 1.3 * StFl.calc_Vstall(self.Atmos.dens_trop_alt(Alt),self.weight,Craft.mainwing.Area,Craft.CLmax).to("meter/second")
        self.V_F = 1.23 * StFl.calc_Vstall(self.Atmos.dens_trop_alt(Alt),self.weight,Craft.mainwing.Area,Craft.CLmax).to("meter/second")
        self.V_TD = 1.15 * StFl.calc_Vstall(self.Atmos.dens_trop_alt(Alt),self.weight,Craft.mainwing.Area,Craft.CLmax).to("meter/second")
        pass

    def Flare(self):
        radius = self.LoadFactor.getRadiusPullUp(self.n_flare,self.V_F)
        Str = radius * numpy.sin(numpy.deg2rad(self.apprachA))
        htr = radius - (radius * numpy.cos(numpy.deg2rad(self.apprachA)))
        #print(Str,htr)
        print(f'Flare distance of {Str.to("meter")} at V_flare: {self.V_F}')
        return([Str,htr])
    
    def Approach(self,ROH):
        Sa = ROH / numpy.tan(numpy.deg2rad(self.apprachA))
        print(f'Approach distance {Sa} as V_app: {self.V_A}')
        return Sa

    def GroundRoll(self):
        weight = self.weight
        wingArea = self.Craft.mainwing.Area
        dens = self.Atmos.dens_trop_alt(self.Alt)
        Lift = StFl.calc_Lift(self.Craft.CLrolling,0.7 * self.V_TD, wingArea,dens).to("newton")
        Drag = StFl.calc_Drag(self.Craft.Cd0,0.7 * self.V_TD,wingArea,dens,self.Craft.K,self.Craft.CLrolling)

        RollDist =  (1.69 * numpy.power(weight,2))/(self.g * dens * wingArea * self.Craft.CLmax * (Drag + (self.RollFricCoe*(weight - Lift))))
        #print((self.g * dens * wingArea * self.Craft.CLmax * (thrust - Drag - (self.RollFricCoe*(weight - Lift)))))
        #print(f'GR: {RollDist}, w: {weight}, wa: {wingArea}, Dens: {dens}, thr: {thrust}, L: {Lift}, D: {Drag}, COE: {self.RollFricCoe}, g: {self.g}')
        print(f'Rolling Distance of {RollDist.to("meter")} at V_TD: {self.V_TD}')
        return (RollDist)
    
    def DoLanding(self):
        ToD = 0
        remainingH = self.obsH

        ToD += self.GroundRoll()

        trans = self.Flare()
        ToD += trans[0]
        remainingH -= trans[1]

        ToD += self.Approach(remainingH)
        return(ToD)

if __name__ == "__main__":
    OppaStoppa = Craft("OppaStoppa")
    OppaStoppa.Atmosphere = Atmosphere.Atmosphere(300,286.21, 9.77774,1.19,170,OppaStoppa.ur)
    ur = OppaStoppa.ur
    atmo = OppaStoppa.Atmosphere
    OppaStoppa.weight_empty = 4450 * 9.81 * ur.newton
    OppaStoppa.weight_takeoff = 5225 * 9.81 * ur.newton

    OppaStoppa.CLmax = 1.45
    OppaStoppa.CLrolling = 0.43 #assuming AOA of 2.5 degrees

    """Defining the draggy components of our craft"""
    #Wing defined: NAME, AIRFOIL, SWEEP, AREA, SPAN, CHORD, angleZeroLift, AngleStall, TC, XC, AreaWIngObscured, atmosphere
    N4312Wing = Wing3d("Oppa Main Wing","NACA 4312",34.87,25.26,9.14,1.76,-4,17,0.12,0.3,6.56,atmo)

    MainWing = N4312Wing

    HorizontalTail = Wing3d("HT","NACA 0012",26.57,4.58,2.44,0.915,0,15,0.12,0.3,0.12,atmo)
    VerticalTail = Wing3d("VT","NACA 0012",26.57,4.58/2,2.44/2,0.915,0,15,0.12,0.3,0,atmo)

    #Fuselage defined: NAME, Length, AreaTop, AreaSide, maxCrossSectionArea, Interf, MainwingArea, atmosphere
    MainFuselage = Fuselage("Oppa Fuselage",7.51,8.11,5.24,1.00,1.0,MainWing.Area,atmo)

    #Gear defined: NAME, CD_component, FrontalArea, MainwingArea, Interf, atmosphere
    TailGear = FixedGear("Tail Gear",0.25,0.196129,MainWing.Area,1.2,atmo)

    """Defining our engines"""
    #Engine defined: NAME, TSFC, BSFC, MaxThrust, MaxPower, efficency.
    #Note that for a turbojet we dont really need BSFC or power
    WilliamsFJ33 = Engine("Willams FJ33",13.77,0,8210,0,0.9,ur)
    OppaStoppa.powertrain = [WilliamsFJ33,WilliamsFJ33]
    OppaStoppa.mainwing = MainWing

    OppaStoppa.dragcomponents = [MainWing,MainFuselage,HorizontalTail,VerticalTail,TailGear]
    OppaStoppa.PosStruturalNlimit = 5 #Set Structual N limits
    OppaStoppa.NegStruturalNlimit = 2
    OppaStoppa.compute_components()

    Takeoff = Takeoff(OppaStoppa, 300,1.15,5)
    print("Total Takeoff: " + str(Takeoff.DoTakeoff().to("meter")))
    print()
    Landing = Landing(OppaStoppa,300,1.1,1,OppaStoppa.weight_takeoff)
    print("Landing " + str(Landing.DoLanding().to("meter")))