import math
import matplotlib.pyplot as plt
import numpy as np
from Craft import *
from Atmosphere import Atmosphere

class CraftStatistics():
    StatsCraft: Craft
    Active_Atmosphere: Atmosphere

    def __init__(self, Craft: Craft) -> None:
        self.StatsCraft = Craft
        self.Active_Atmosphere = Craft.Atmosphere

    def get_ThrustAvailable_jet(self,alt):
        """Taking into consideration all sources of thrust on the craft, finds the thrust available at the given altitude.
        Assumes m = 1

        Args:
            alt (_type_): Altitude above sea level in meters
        """        
        maxthurust = self.StatsCraft.get_max_thrust()
        aval_thrust = self.Active_Atmosphere.jet_thrust_available_alt(maxthurust,alt,1)
        return aval_thrust
    
    def graph_ThrustAvailable(self,Alt_lower_lim, Alt_upper_lim,numPoints)-> plt:
        """Returns a Matplotlib plot of thrust available vs alt

        Args:
            Alt_lower_lim (_type_): Lower alt limit in meters above sea
            Alt_upper_lim (_type_): Upper alt limit
            numPoints (_type_): number of points

        Returns:
            plt: Plot containing Thrust available vs altitude
        """        
        alt_array = np.linspace(Alt_lower_lim,Alt_upper_lim,num=numPoints) #Array of numbers between lower and upper (inclusive)
        dens_array = np.zeros(alt_array.shape)
        thrust_array = np.zeros(alt_array.shape) #Make new array for thrust
        temp_array = np.zeros(alt_array.shape) #Make new array for thrust


        for i in enumerate(alt_array): #Iterate over altitude and calculate Thrust
            thrust_array[i[0]] = self.get_ThrustAvailable_jet(i[1])
            dens_array[i[0]] = self.Active_Atmosphere.dens_trop_alt(i[1])
            temp_array[i[0]] = self.Active_Atmosphere.temp_trop(i[1])



        fig, ax = plt.subplots()
        #ax.plot(alt_array,dens_array,linewidth=1,label="Dens")
        #ax.plot(alt_array,temp_array,linewidth=1,label="temp")

        ax.plot(alt_array,thrust_array,linewidth=2,label="Thrust Available (N)")
        ax.set(xlabel='Altitude (meters)', ylabel='Thrust available(N)',title='Altitude vs Thrust Available')

        return fig
    

if __name__ == "__main__":
    OppaStoppa = Craft("OpptaStoppa")
    OppaStoppa.Atmosphere = Atmosphere(300,286.21, 9.77774,1.19,76)
    atmo = OppaStoppa.Atmosphere
    OppaStoppa.weight_empty = 4450 * 9.81
    OppaStoppa.weight_takeoff = 5225 * 9.81


    """Defining the draggy components of our craft"""
    #Wing defined: NAME, AIRFOIL, SWEEP, AREA, SPAN, CHORD, angleZeroLift, AngleStall, TC, XC, AreaWIngObscured, atmosphere
    MainWing = Wing3d("Oppa Main Wing","NACA 4312",34.87,25.26,9.14,1.76,-4,17,0.12,0.3,6.56,atmo)
    HorizontalTail = Wing3d("HT","NACA 0012",26.57,4.58,2.44,0.915,0,15,0.12,0.3,0.12,atmo)
    VerticalTail = Wing3d("VT","NACA 0012",26.57,4.58/2,2.44/2,0.915,0,15,0.12,0.3,0,atmo)

    #Fuselage defined: NAME, Length, AreaTop, AreaSide, maxCrossSectionArea, Interf, MainwingArea, atmosphere
    MainFuselage = Fuselage("Oppa Fuselage",7.51,8.11,5.24,1.00,1.0,MainWing.Area,atmo)

    #Gear defined: NAME, CD_component, FrontalArea, MainwingArea, Interf, atmosphere
    TailGear = FixedGear("Tail Gear",0.25,0.196129,MainWing.Area,1.2,atmo)

    """Defining our engines"""
    #Engine defined: NAME, TSFC, BSFC, MaxThrust, MaxPower, efficency.
    #Note that for a turbojet we dont really need BSFC or power
    WilliamsFJ33 = Engine("Willams FJ33",13.77,0,8210,0,0.9)
    OppaStoppa.powertrain = [WilliamsFJ33,WilliamsFJ33]

    OppaStoppa.dragcomponents = [MainWing,MainFuselage,HorizontalTail,VerticalTail,TailGear]

    OppaStoppa.compute_components()


    MyCraftStats = CraftStatistics(OppaStoppa)
    print(OppaStoppa.get_max_thrust())
    print(MyCraftStats.get_ThrustAvailable_jet(2))
    ThrustAvailableCurve = MyCraftStats.graph_ThrustAvailable(0,10000,1000)
    ThrustAvailableCurve.legend()
    plt.show()
    