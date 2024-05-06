from craft import Craft
from Components import *
import StFl
import Atmosphere
import Drag
import CraftStatistics
import LoadFactor
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

class Takeoff:
    def __init__(self, Craft:Craft,Alt,n_takeoff,climbangle,REPLACEWEIGHT = 0) -> None:
        """This class takes in an aircraft and conditons about takeoff and then is able to determine the length of takeoff

        Args:
            Craft (Craft): The craft to use for takeoff
            Alt (_type_): Altitide of takeoff (airport elevation)
            n_takeoff (_type_): The n limit for the aircraft
            climbangle (_type_): Climb angle of takeoff (degrees)
            REPLACEWEIGHT : Weight to use for takeoff. Defaults to 0, meaning the crafts takeoffwieght varable is used. 
        """        
        self.Craft = Craft
        self.CraftStats = CraftStatistics.CraftStatistics(Craft)
        self.LoadFactor = LoadFactor.LoadFactor(Craft)
        #print(Craft)
        #print(Craft.ur)
        self.ur = Craft.ur
        self.ActiveAtmospere = Craft.Atmosphere
        self.Alt = Alt * self.ur.m
        self.g = 9.81 * self.ur.m /self.ur.second**2
        self.RollFricCoe = 0.04
        if REPLACEWEIGHT == 0:
            self.weight = Craft.weight_takeoff
        else:
            self.weight = REPLACEWEIGHT
        self.V_LO = 1.2 * StFl.calc_Vstall(self.ActiveAtmospere.dens_trop_alt(Alt),self.weight,Craft.mainwing.Area,Craft.CLmax).to("meter/second")
        self.n_to = n_takeoff
        self.height_Obs = 10.668 * self.ur.m #35 ft
        self.climbAngle = climbangle
        self.totalDist = 0 * self.ur.m
    
    def GroundRoll(self):
        """Calculates ground roll distance

        Returns:
            _type_: Rolling distance in meters
        """        
        weight = self.weight
        wingArea = self.Craft.mainwing.Area
        dens = self.ActiveAtmospere.dens_trop_alt(self.Alt)
        thrust = self.CraftStats.get_thrustAvailable(self.Alt,self.V_LO)
        #print("THRUST: " + str(thrust.to("newton")))
        Lift = StFl.calc_Lift(self.Craft.CLrolling,0.7 * self.V_LO,wingArea,dens).to("newton")
        Drag = StFl.calc_Drag(self.Craft.Cd0,0.7 * self.V_LO,wingArea,dens,self.Craft.K,self.Craft.CLrolling)

        RollDist =  (1.44 * numpy.power(weight,2))/(self.g * dens * wingArea * self.Craft.CLmax * (thrust - Drag - (self.RollFricCoe*(weight - Lift))))
        #print((self.g * dens * wingArea * self.Craft.CLmax * (thrust - Drag - (self.RollFricCoe*(weight - Lift)))))
        #print(f'GR: {RollDist}, w: {weight}, wa: {wingArea}, Dens: {dens}, thr: {thrust}, L: {Lift}, D: {Drag}, COE: {self.RollFricCoe}, g: {self.g}')
        #print(f'Rolling Distance of {RollDist.to("meter")}')
        return (RollDist)
    
    def Transition(self):
        """Calculates the transiton period of the takeoff

        Returns:
            Array of transition distance and height [horizontal distance, vertical distance]
        """        
        radius = self.LoadFactor.getRadiusPullUp(self.n_to,self.V_LO)
        Str = radius * numpy.sin(numpy.deg2rad(self.climbAngle))
        htr = radius - (radius * numpy.cos(numpy.deg2rad(self.climbAngle)))
        #print(Str,htr)
        #print(f'Tranition distance of {Str.to("meter")} at VLO: {self.V_LO}')
        return([Str,htr])
    
    def AirDist(self,remainingObsH):
        """Calculates the in air section of the takeoff

        Args:
            remainingObsH (_type_): How much remaining height to clear the takeoff obsticle height 

        Returns:
            _type_: Horzontal distance of air section
        """        
        Sa = remainingObsH / numpy.tan(numpy.deg2rad(self.climbAngle))
        #print(Sa)
        #print(f'Air Distance of {Sa.to("meter")}')
        return Sa
    
    def DoTakeoff(self):
        """Performs the entire takeoff and prints the results.
        """        
        ToD = 0
        remainingH = self.height_Obs

        rd = self.GroundRoll().to("meter")
        print(f'Ground Roll Distance: {rd}')
        ToD += rd

        trans = self.Transition()
        print(f'Transition distance: {trans[0]}')
        ToD += trans[0]
        remainingH -= trans[1]

        aird = self.AirDist(remainingH)
        print(f'Air Distance: {aird}')
        ToD += aird
        self.totalDist = ToD
        return(ToD)
    
    def get_DistTraveled(self):
        ToD = 0
        remainingH = self.height_Obs

        rd = self.GroundRoll().to("meter")
        #print(f'Ground Roll Distance: {rd}')
        ToD += rd

        trans = self.Transition()
        #print(f'Transition distance: {trans[0]}')
        ToD += trans[0]
        remainingH -= trans[1]

        aird = self.AirDist(remainingH)
        #print(f'Air Distance: {aird}')
        ToD += aird
        
        linedist = [rd , numpy.sqrt((trans[0] + aird)**2 + self.height_Obs**2)]
        return(linedist)
    
    def Graph(self)->plt:
        """Graphs the trajectory of takeoff

        Returns:
            _type_: Figure containing graph
        """        
        fig, ax = plt.subplots()
        ax.set_aspect(10)

        ToD = 0
        remainingH = self.height_Obs

        GroundRoll = self.GroundRoll()

        trans = self.Transition()
        Transitiondist = trans[0]
        remainingH -= trans[1]
        #print(f'Height from tans{trans[1]}')
        #ax.hlines(trans[1],0,1000,colors=["orange"],linestyles="dashed")
        

        AirDist = self.AirDist(remainingH)

        totalDist = GroundRoll + Transitiondist + AirDist
        #ax.vlines(totalDist,0,50)
        airtime = numpy.linspace(GroundRoll+Transitiondist,totalDist,5)
        height = numpy.linspace(trans[1].to("meters"),self.height_Obs.to("meters"),5)
        ax.plot(airtime,height,linestyle="--",color="red",label="In air")

        #Circle for transition
        radius = self.LoadFactor.getRadiusPullUp(self.n_to,self.V_LO).magnitude
        ellipse = Arc((GroundRoll,radius),radius*2,radius*2,theta1=270,theta2=270+self.climbAngle,color='orange',linestyle="--",label="Transition")
        ax.add_patch(ellipse)
       
        #plot ground roll
        ax.hlines(0,0,GroundRoll,colors=["green"],linestyles="dashed",label = "Ground Roll")
        ax.set_xbound(0,totalDist)
        ax.set_ybound(-0.5,self.height_Obs.magnitude + 1)
        ax.set_xlabel("Distance from start (meters)")
        ax.set_ylabel("Height (meters)")
        ax.set_title("Takeoff Diagram")

        #Annotate for Vel
        bbox = dict(boxstyle="round", fc="0.8")
        arrowprops = dict(arrowstyle="->",connectionstyle="angle,angleA=15,angleB=-30,rad=10")
        VelLO = r'$V_{LO}$ ' + str(round(self.V_LO,1))
        ax.annotate(VelLO,(GroundRoll,0),arrowprops=arrowprops,bbox=bbox,xytext=(-45,15), textcoords='offset points')

        #InfoBox
        textstr = f'Takeoff Parameters \n Altitude: {self.Alt} \n ClimbAng: {self.climbAngle}° \n Obstacle: {self.height_Obs}'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.05, 0.8, textstr, transform=ax.transAxes, fontsize=10, bbox=props)

        return fig
    

    
class Landing:
    def __init__(self, Craft: Craft, Alt,n_flare,ApproachAngle,weight) -> None:
        """Takes in a craft and info about Landing to calculate its length

        Args:
            Craft (Craft): Craft to Land
            Alt (_type_): Altitude (Elevation of airport)
            n_flare (_type_): n value to use for flare
            ApproachAngle (_type_): Apprach angle (glide angle)
            weight (_type_): wight to use at takeoff (newtons)
        """        
        self.Craft = Craft
        self.CraftStats = CraftStatistics.CraftStatistics(Craft)
        self.LoadFactor = LoadFactor.LoadFactor(Craft)
        self.ur = Craft.ur
        self.Atmos = Craft.Atmosphere
        self.weight=weight
        self.Alt = Alt * self.ur.m
        self.n_flare = n_flare
        self.apprachA = ApproachAngle
        self.g = 9.81 * self.ur.m /self.ur.second**2
        self.RollFricCoe = 0.4
        self.obsH = 10.668 * self.ur.m #35 ft
        self.V_A = 1.3 * StFl.calc_Vstall(self.Atmos.dens_trop_alt(Alt),self.weight,Craft.mainwing.Area,Craft.CLmax).to("meter/second")
        self.V_F = 1.23 * StFl.calc_Vstall(self.Atmos.dens_trop_alt(Alt),self.weight,Craft.mainwing.Area,Craft.CLmax).to("meter/second")
        self.V_TD = 1.15 * StFl.calc_Vstall(self.Atmos.dens_trop_alt(Alt),self.weight,Craft.mainwing.Area,Craft.CLmax).to("meter/second")
        self.totalDist = 0 * self.ur.m
        pass

    def Flare(self):
        """Calcualtes flare distance

        Returns:
            The horizontal and vertical components of the flare [horiz, vertical]
        """        
        radius = self.LoadFactor.getRadiusPullUp(self.n_flare,self.V_F)
        Str = radius * numpy.sin(numpy.deg2rad(self.apprachA))
        htr = radius - (radius * numpy.cos(numpy.deg2rad(self.apprachA)))
        #print(Str,htr)
        #print(f'Flare distance of {Str.to("meter")} at V_flare: {self.V_F}')
        #print(f'Flare looses {htr}')
        return([Str,htr])
    
    def Approach(self,ROH):
        """Simulates approach 

        Args:
            ROH (_type_): How much of the landing obsticle height is not taken up by our flare (meter)

        Returns:
            _type_: Distance of Apprach
        """        
        Sa = ROH / numpy.tan(numpy.deg2rad(self.apprachA))
        #print(f'Approach distance {Sa} as V_app: {self.V_A}')
        return Sa

    def GroundRoll(self):
        """Simulates ground roll

        Returns:
            _type_: Ground roll distance. 
        """        
        weight = self.weight
        wingArea = self.Craft.mainwing.Area
        dens = self.Atmos.dens_trop_alt(self.Alt)
        Lift = StFl.calc_Lift(self.Craft.CLrolling,0.7 * self.V_TD, wingArea,dens).to("newton")
        Drag = StFl.calc_Drag(self.Craft.Cd0,0.7 * self.V_TD,wingArea,dens,self.Craft.K,self.Craft.CLrolling)

        RollDist =  (1.69 * numpy.power(weight,2))/(self.g * dens * wingArea * self.Craft.CLmax * (Drag + (self.RollFricCoe*(weight - Lift))))
        #print((self.g * dens * wingArea * self.Craft.CLmax * (thrust - Drag - (self.RollFricCoe*(weight - Lift)))))
        #print(f'GR: {RollDist}, w: {weight}, wa: {wingArea}, Dens: {dens}, thr: {thrust}, L: {Lift}, D: {Drag}, COE: {self.RollFricCoe}, g: {self.g}')
        #print(f'Rolling Distance of {RollDist.to("meter")} at V_TD: {self.V_TD}')
        return (RollDist)
    
    def DoLanding(self):
        """Performs all the parts of a landing and prints the results.
        """        
        ToD = 0
        remainingH = self.obsH

        gr = self.GroundRoll().to("meter")
        print(f'Landing ground roll dist: {gr}')
        ToD += gr

        trans = self.Flare()
        print(f'Transition Distance: {trans[0]}')
        ToD += trans[0]
        remainingH -= trans[1]

        approach = self.Approach(remainingH)
        print(f'Approach Distance: {approach}')
        ToD += approach
        self.totalDist = ToD
        return(ToD)
    
    def Graph(self)->plt:
        """Performs a landing and outputs a graph of the landing

        Returns:
            plt: matplotlib figure of landing
        """    
        #plt.rcParams['text.usetex'] = True   
        fig, ax = plt.subplots()
        
        ax.set_aspect(10)

        ToD = 0
        remainingH = self.obsH

        GroundRoll = self.GroundRoll()

        trans = self.Flare()
        flaredist = trans[0]
        remainingH -= trans[1]
        #print(f'Height from tans{trans[1]}')
        #ax.hlines(trans[1],0,1000,colors=["orange"],linestyles="dashed")
        

        AirDist = self.Approach(remainingH)

        totalDist = GroundRoll + flaredist + AirDist
        #ax.vlines(totalDist,0,50)
        airtime = numpy.linspace(GroundRoll+flaredist,totalDist,5)
        height = numpy.linspace(trans[1].to("meters"),self.obsH.to("meters"),5)
        ax.plot(airtime,height,linestyle="--",color="red",label="Approach")

        #Circle for transition
        radius = self.LoadFactor.getRadiusPullUp(self.n_flare,self.V_F).magnitude
        ellipse = Arc((GroundRoll,radius),radius*2,radius*2,theta1=270,theta2=270+self.apprachA,color='orange',linestyle="--",label="Flare")
        ax.add_patch(ellipse)
       
        #plot ground roll
        ax.hlines(0,0,GroundRoll,colors=["green"],linestyles="dashed",label = "Ground Roll")
        ax.set_xbound(0,totalDist)
        ax.set_ybound(-0.5,self.obsH.magnitude + 1)
        ax.set_xlabel("Distance from stop (meters)")
        ax.set_ylabel("Height (meters)")
        ax.set_title("Landing Diagram")

        #Adding annotations to indicate velocity 
        bbox = dict(boxstyle="round", fc="0.8")
        arrowprops = dict(arrowstyle="->",connectionstyle="angle,angleA=15,angleB=-30,rad=10")
        VelApprach = r'$V_A$ ' + str(round(self.V_A,1))
        ax.annotate(VelApprach,(numpy.average(airtime),numpy.average(height)),arrowprops=arrowprops,bbox=bbox,xytext=(-45,15), textcoords='offset points')

        VelFlare = r'$V_F$ ' + str(round(self.V_F,1))
        ax.annotate(VelFlare,(GroundRoll.magnitude + radius*numpy.sin(numpy.deg2rad(self.apprachA/2)),radius - radius*numpy.cos(numpy.deg2rad(self.apprachA/2))),arrowprops=arrowprops,bbox=bbox,xytext=(-40,20), textcoords='offset points')

        VelTouch = r'$V_{TD}$ ' + str(round(self.V_TD,1))
        ax.annotate(VelTouch,(GroundRoll.magnitude,0),arrowprops=arrowprops,bbox=bbox,xytext=(-100,20), textcoords='offset points')

        #Add info box
        textstr = f'Landing Parameters \n Altitude: {self.Alt} \n ApprAng: {self.apprachA:0.2f}° \n Obstacle: {self.obsH}'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.05, 0.7, textstr, transform=ax.transAxes, fontsize=10, bbox=props)

        return fig

if __name__ == "__main__":
    #TESTING OF CLASS, NOT ACTUALL CRAFT PROPTERIES 
    OppaStoppa = Craft("OppaStoppa")
    ur = OppaStoppa.ur
    OppaStoppa.Atmosphere = Atmosphere.Atmosphere(300,286.21, 9.77774,1.19,15,ur)
    atmo = OppaStoppa.Atmosphere
    OppaStoppa.weight_empty = 28 * 9.81 * ur.newtons
    OppaStoppa.weight_takeoff = 28 * 9.81 * ur.newtons
    OppaStoppa.CLmax = 1.45
    OppaStoppa.CLrolling = 0.43 #assuming AOA of 2.5 degrees


    """Defining the draggy components of our craft"""
    #Wing defined: NAME, AIRFOIL, SWEEP, AREA, SPAN, CHORD, angleZeroLift, AngleStall, TC, XC, AreaWIngObscured, atmosphere
    N4312Wing = Wing3d("Oppa Main Wing","NACA 4312",34.83,1.29,2.34,.532,-1,15,.12,.3,.02,atmo)
    MainWing = N4312Wing

    HorizontalTail = Wing3d("HT","NACA 0012",26.56,.206,.73,.274,0,15,0.12,0.03,0.12,atmo)
    VerticalTail = Wing3d("VT","NACA 0012",26.56,.206,.73,.274,0,15,0.12,0.03,0.12,atmo)

    #Fuselage defined: NAME, Length, AreaTop, AreaSide, maxCrossSectionArea, Interf, MainwingArea, atmosphere
    MainFuselage = Fuselage("Oppa Fuselage",1.926,.48,.39,.019,.02,MainWing.Area,atmo)

    #Gear defined: NAME, CD_component, FrontalArea, MainwingArea, Interf, atmosphere
    TailGear = FixedGear("Tail Gear",0.25,0.002,MainWing.Area,1.2,atmo)

    """Defining our engines"""

    

    #Define Prop
    prop = Propeller("16X8",0.4064,8,ur)
    def thrustproppoly(advr):
        try:
            advr  = advr.magnitude
        except:
            advr = advr
        thr = 0.122 - (0.0138 * advr) + (0.0709* advr**2) - (0.287 * advr**3) + (0.137 * advr**4)
        return thr
    prop.thrust_polynomal_func = thrustproppoly

    def powerproppoly(advr):
        try:
            advr  = advr.magnitude
        except:
            advr = advr
        pwr = 0.0583 + (0.0127 * advr) + (0.23* advr**2) - (0.416 * advr**3) + (0.151 * advr**4)
        return pwr
    prop.power_polynomal_func = powerproppoly

    Motor = ElectricMotor("V804 KV170",5200,0.90,8000,prop,ur)

    Battery = Battery(139.76,4.32 * 2,45,ur)
    #Battery.print_state()

    OppaStoppa.powertrain = [Motor,Motor]
    OppaStoppa.mainwing = MainWing

    OppaStoppa.dragcomponents = [MainWing,MainFuselage,HorizontalTail,VerticalTail]

    OppaStoppa.compute_components()
    OppaStoppa.PosStruturalNlimit = 5
    OppaStoppa.NegStruturalNlimit = -3


    Takeoff = Takeoff(OppaStoppa, 300,1.5,5)
    print("Total Takeoff: " + str(Takeoff.DoTakeoff().to("meter")))
    print()
    Landing = Landing(OppaStoppa,300,1.1,3,OppaStoppa.weight_takeoff)
    print("Landing " + str(Landing.DoLanding().to("meter")))

    TOG = Takeoff.Graph()
    TOG.legend()
    #plt.show()

    LOG = Landing.Graph()
    #LOG.legend()
    plt.show()