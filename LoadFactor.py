from craft import Craft
import CraftStatistics
import Atmosphere
from Components import *
import math
import numpy as np
import matplotlib.pyplot as plt


class LoadFactor():
    LoadCraft: Craft
    ActiveAtmos: Atmosphere

    def __init__(self,craft:Craft) -> None:
        """LoadFactor takes in a craft and performs the analysis requred to generate a V-n Diagram as well as find various other performace charecteristics 

        Args:
            craft (Craft): Craft to be used 
        """        
        self.LoadCraft: Craft= craft
        self.ActiveAtmos = craft.Atmosphere
        self.ur = self.ActiveAtmos.ur
        self.g = 9.81 * self.ur.meters / self.ur.seconds ** 2

    def weight_from_str(self, WEIGHT: str)-> float:
        """Gets the crafts weight from a string.

        Args:
            WEIGHT (str): Either "TAKEOFF" for the provided crafts takeoff weight, "EMPTY" for empty weight , "AVE" average the TOW and EW, or a float (ie: 10.1). When a float is provided it will be cast from string to the float.

        Returns:
            float: weight in Newtons
        """
        if WEIGHT.upper() == "TAKEOFF":
            weight = self.LoadCraft.weight_takeoff
        elif WEIGHT.upper() == "EMPTY":
            weight = self.LoadCraft.weight_empty
        elif WEIGHT.upper() == "AVE":
            weight = (self.LoadCraft.weight_empty + self.StatsCraft.weight_takeoff)/2
        else:
            try:
                weight = float(WEIGHT)
            except:
                weight = 0
                print("Conversion to float failed")
        return weight

    def getRadiusPullUp(self,LoadFactor,velocity):
        """Finds the radius of a pull up for the given velocity and loadfactor

        Args:
            LoadFactor (_type_): LoadFactor
            velocity (_type_): Velocity m/s

        Returns:
            _type_: Radius of pullup (or pushup) (m)
        """        
        RPU = numpy.power(velocity,2) / (self.g * (LoadFactor - 1))
        return RPU

    def getRadiusPullDown(self,LoadFactor,velocity):
        """Finds the radius of a pull down 

        Args:
            LoadFactor (_type_): _description_
            velocity (_type_): Velocity m/s

        Returns:
            _type_: Radius of pull down (or pushdown) (m)
        """        
        RPD = numpy.power(velocity,2) / (self.g * (LoadFactor + 1))
        return RPD
    
    def getRadiusLevelTurn(self,LoadFactor,velocity):
        """gets radius of a level turn

        Args:
            LoadFactor (_type_): Load Factor
            velocity (_type_): Velocity (m/s)

        Returns:
            _type_: Radius of turn (m)
        """        
        RLT = numpy.power(velocity,2) / (self.g * numpy.sqrt(numpy.power(LoadFactor,2) - 1))
        return RLT
    
    def getn(self, Vinf, CL, alt, Weight):
        """Finds the aerodynamicly possible n value

        Args:
            Vinf (_type_): Velocity
            CL (_type_): Coeff of lift
            alt (_type_): Altitude
            Weight (_type_): Weight of craft (Newton)
        """        
        dens = self.ActiveAtmos.dens_trop_alt(alt)
        weight = self.weight_from_str(Weight)
        n = (0.5 * dens * (Vinf ** 2) * self.LoadCraft.mainwing.Area * CL) / weight
        #print("0.5 * " + str(dens) + )
        #print(f'Nvalue: {n}\nDens {dens}\nWeight {weight}\nVel^2 {numpy.power(vinf,2)}\nWing Area {self.LoadCraft.mainwing.Area}')
        return(n)
    
    def getAeroStructMeetV(self, NstructPOS, alt, CLMAX, WEIGHT):
        """Find the Velocity where the Aero limits of craft meet the structural limit

        Args:
            NstructPOS (_type_): Structural n limit
            alt (_type_): Altidude
            CLMAX (_type_): Max CL value
            WEIGHT (_type_): Weight of craft (newton)

        Returns:
            _type_: Velocity where n values meet
        """        
        weight = self.weight_from_str(WEIGHT)
        dens = self.ActiveAtmos.dens_trop_alt(alt)
        vinf = numpy.sqrt(((2 * NstructPOS)/(dens * CLMAX)) * (weight/self.LoadCraft.mainwing.Area))
        return vinf

    def genVnDiagram(self,Vmin,Vmax,numpoints,alt,Weight,clmin,clmax)-> plt:
        """Returns a matplotlib figure of the Vn diagram for the craft

        Args:
            Vmin (_type_): Lower velocity limit of the graph, should be below stall speed of craft    
            Vmax (_type_): Upper velocity limit of graph
            numpoints (_type_): number of points to plot 
            alt (_type_): Altitude to use
            Weight (_type_): Weight of craft (newton)
            clmin (_type_): Max CL value of craft
            clmax (_type_): Max CL of craft in a push down (negative AoA)

        Returns:
            _type_: Figure of V-n Chart
        """        
        weight = self.weight_from_str(Weight)
        atmos = self.LoadCraft.Atmosphere
        nplus = self.LoadCraft.PosStruturalNlimit
        nminus = self.LoadCraft.NegStruturalNlimit
        #Clmax = self.LoadCraft.MainWing.CLmax
        #Clmin = self.LoadCraft.MainWing.CLmin

        vstar = self.getAeroStructMeetV(nplus,alt,clmax,Weight)

        vArr = np.linspace(Vmin,Vmax,numpoints)
        nplusarr = np.zeros(vArr.shape)
        nminusarr = np.zeros(vArr.shape)

        for vel in enumerate(vArr):
            #print("HELLO " + str(self.getn(vel[1],clmax,alt,Weight)))
            nplusarr[vel[0]] = self.getn(vel[1],clmax,alt,Weight).magnitude
            nminusarr[vel[0]] = (-1) * self.getn(vel[1],clmin,alt,Weight).magnitude

        locofn1p = np.abs(nplusarr - 1)
        LocPn1 = locofn1p.argmin()
        locofn1n = np.abs(nminusarr + 1)
        LocNn1 = locofn1n.argmin()

        fig, ax = plt.subplots()
        maskedNplus = np.ma.masked_where((nplusarr < 1) | (nplusarr > self.LoadCraft.PosStruturalNlimit),nplusarr) #Mask to just data below nmax
        ABOVENplus = np.ma.masked_where((nplusarr < self.LoadCraft.PosStruturalNlimit),nplusarr) #Mask to just data below nmax

        ax.plot(vArr,maskedNplus,linewidth=1,label="Positive Aero Limit",color='blue')
        ax.plot(vArr,ABOVENplus,linewidth=1,color='blue',linestyle="--")

        maskedNminus = np.ma.masked_where((nminusarr > -1) | (nminusarr < (-1 *self.LoadCraft.NegStruturalNlimit)),nminusarr)
        BELOWNminus = np.ma.masked_where((nminusarr > (-1 *self.LoadCraft.NegStruturalNlimit)),nminusarr)

        ax.plot(vArr,maskedNminus,linewidth=1,label="Negative Aero Limit",color='orange')
        ax.plot(vArr,BELOWNminus,linewidth=1,color='orange',linestyle='--')


        ax.hlines(self.LoadCraft.PosStruturalNlimit,Vmin,Vmax,colors="red",linestyles="dotted",label="+/- Structual Limit")
        ax.vlines(vArr[LocPn1],0,1,colors="blue")

        ax.hlines((-1) * self.LoadCraft.NegStruturalNlimit,Vmin,Vmax,colors="red",linestyles="dotted")
        ax.vlines(vArr[LocNn1],-1,0,colors="orange")


        ax.vlines(vstar,0,nplus,colors="grey",linestyles="dotted",label="V*")

        ax.set_ybound((-1) * (self.LoadCraft.NegStruturalNlimit+1), self.LoadCraft.PosStruturalNlimit + 1 )

        
        strTitle = "V-n Diagram"
        ax.set(xlabel='Velocity (m/s)', ylabel='LoadFactor',title=strTitle)
        ax.spines.bottom.set_position('zero')
        
        ax.spines.top.set_color('none')
        ax.legend(loc='center right')

        return fig



 

if __name__ == "__main__":
    OppaStoppa = Craft("OppaStoppa")
    OppaStoppa.Atmosphere = Atmosphere(300,286.21, 9.77774,1.19,170,OppaStoppa.ur)
    ur = OppaStoppa.ur
    atmo = OppaStoppa.Atmosphere
    OppaStoppa.weight_empty = 4450 * 9.81 * ur.newton
    OppaStoppa.weight_takeoff = 5225 * 9.81 * ur.newton


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

    #Structural Limits of craft at n 5,-2
    LFA = LoadFactor(OppaStoppa)
    vinf = LFA.ActiveAtmos.Vinfinity
    alt = LFA.ActiveAtmos.Altitude

    SradPD = LFA.getRadiusPullDown(LFA.LoadCraft.NegStruturalNlimit,vinf)
    SradPU = LFA.getRadiusPullUp(LFA.LoadCraft.PosStruturalNlimit,vinf)
    SradLT = LFA.getRadiusLevelTurn(LFA.LoadCraft.PosStruturalNlimit,vinf)

    CLMAX = 1.45
    CLMIN = 0.86 #Actually -0.86

    wstr = "TAKEOFF"

    nCLMAX = LFA.getn(vinf,CLMAX,alt,wstr)
    nCLMIN = LFA.getn(vinf,CLMIN,alt,wstr)

    AradPD = LFA.getRadiusPullDown(nCLMIN,vinf)
    AradPU = LFA.getRadiusPullUp(nCLMAX,vinf)
    AradLT = LFA.getRadiusLevelTurn(nCLMAX,vinf)

    print("Load at:" )
    print(" Craft weight used: " + wstr +"->" + str(LFA.weight_from_str(wstr)) )
    LFA.ActiveAtmos.printAtmos()

    print("\nCrafts Structural n+, n-: " + str(OppaStoppa.PosStruturalNlimit) +", " + str(OppaStoppa.NegStruturalNlimit))
    print(" Radius for PullUp: " + str(SradPU) )
    print(" Radius for PullDown: " + str(SradPD) )
    print(" Radius for LevelTurn: " + str(SradLT) )

    print("\nCrafts Aero n+, n-: " + str(nCLMAX.magnitude) +", " + str(nCLMIN.magnitude))
    print(" Radius for PullUp: " + str(AradPU) )
    print(" Radius for PullDown: " + str(AradPD) )
    print(" Radius for LevelTurn: " + str(AradLT) )

    Vmeeting = LFA.getAeroStructMeetV(OppaStoppa.PosStruturalNlimit,alt,CLMAX,wstr)
    print("\nV* = " + str(Vmeeting) )

    fig = LFA.genVnDiagram(40,200,400,300,"TAKEOFF",CLMIN,CLMAX)

    plt.show()



    