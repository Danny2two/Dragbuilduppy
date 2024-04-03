import math
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from craft import *
from Atmosphere import Atmosphere

class CraftStatistics():
    StatsCraft: Craft
    Active_Atmosphere: Atmosphere

    def __init__(self, Craft: Craft) -> None:
        self.StatsCraft = Craft
        self.Active_Atmosphere = Craft.Atmosphere

        print("Initalizing Statistics for Craft: " + self.StatsCraft.name)
        print("CD0: " + str(self.StatsCraft.Cd0))
        print("Ozwald: " + str(self.StatsCraft.mainwing.OswaldE))

    #PART F
    def weight_from_str(self, WEIGHT: str)-> float:
        """Gets the crafts weight from a string.

        Args:
            WEIGHT (str): Either "TAKEOFF" for the provided crafts takeoff weight, "EMPTY" for empty weight , "AVE" average the TOW and EW, or a float (ie: 10.1). When a float is provided it will be cast from string to the float.

        Returns:
            float: weight in Newtons
        """
        if WEIGHT.upper() == "TAKEOFF":
            weight = self.StatsCraft.weight_takeoff
        elif WEIGHT.upper() == "EMPTY":
            weight = self.StatsCraft.weight_empty
        elif WEIGHT.upper() == "AVE":
            weight = (self.StatsCraft.weight_empty + self.StatsCraft.weight_takeoff)/2
        else:
            try:
                weight = float(WEIGHT)
            except:
                weight = 0
                print("Conversion to float failed")
        return weight

    #PART F
    def get_ThrustAvailable_jet(self,alt):
        """Taking into consideration all sources of thrust on the craft, finds the thrust available at the given altitude.
        Assumes m = 1

        Args:
            alt (_type_): Altitude above sea level in meters
        """        
        maxthurust = self.StatsCraft.get_max_thrust()
        aval_thrust = self.Active_Atmosphere.jet_thrust_available_alt(maxthurust,alt,1)
        return aval_thrust
    
    #PART F
    def get_PowerAvailable_jet(self,alt,velocity):
        return (self.get_ThrustAvailable_jet(alt) * velocity)
    
    #PART F
    def get_PowerRequired_alt_jet(self,alt,velocity,WEIGHT: str):
        """Finds power requred for SLF at given conditions

        Args:
            alt (_type_): Altitude above sealevel
            velocity (_type_): Velocity in m/s
            WEIGHT (str): Either "TAKEOFF" for the provided crafts takeoff weight, "EMPTY" for empty weight , "AVE" average the TOW and EW, or a float (ie: 10.1). When a float is provided it will be cast from string to the float.

        Returns:
            Pr: Power requred 
        """        
        weight = self.weight_from_str(WEIGHT)

        dens = self.Active_Atmosphere.dens_trop_alt(alt)
        k = calc_K_value(self.StatsCraft.mainwing.OswaldE,self.StatsCraft.mainwing.AR)
        Pr = calc_PowerReq(dens,velocity,self.StatsCraft.mainwing.Area,self.StatsCraft.Cd0,k,weight)
        return Pr
    
    #PART F
    def graph_PowerAval_vs_PowerReq(self,Alt_lower, Alt_upper,numPoints, Velocity, WEIGHT: str,GRAPH_EXCESS: bool = False,SENDRAW: bool = False):
        """Returns a MPL figure of power requred and available vs alt

        Args:
            Alt_lower (_type_): Altitude above sealevel lower lim
            Alt_upper : Alt upper lim
            velocity (_type_): Velocity in m/s
            WEIGHT (str): Either "TAKEOFF" for the provided crafts takeoff weight, "EMPTY" for empty weight , "AVE" average the TOW and EW, or a float (ie: 10.1). When a float is provided it will be cast from string to the float.

        Returns:
            plt: Plot
        """     
        weight = self.weight_from_str(WEIGHT)
        
        alt_array = np.linspace(Alt_lower,Alt_upper,num=numPoints) #Array of numbers between lower and upper (inclusive)
        prA_array = np.zeros(alt_array.shape)
        prR_array = np.zeros(alt_array.shape)
        for i in enumerate(alt_array): #Iterate over altitude and calculate Thrust
            prA_array[i[0]]= self.get_PowerAvailable_jet(i[1],Velocity) /1000
            prR_array[i[0]]= self.get_PowerRequired_alt_jet(i[1],Velocity,str(weight)) /1000

        if GRAPH_EXCESS:
            prEx_arr = prA_array - prR_array
            if SENDRAW:
                return [alt_array,prEx_arr]
            fig, ax = plt.subplots()
            ax.plot(alt_array,prEx_arr,linewidth=2,label="Power Excess (kW)")
            ax.plot(alt_array,prA_array,linewidth=1,label="Power Available (kW)",linestyle="--")
            ax.plot(alt_array,prR_array,linewidth=1,label="Power Required (kW)",linestyle="--")

            strTitle = "Altitude vs Excess power"
            ax.hlines(0,Alt_lower,Alt_upper,colors="red",linestyles="dotted",label="Zero Excess")
        else:
            fig, ax = plt.subplots()
            strTitle = "Altitude vs Power Available & Required"
            ax.plot(alt_array,prA_array,linewidth=2,label="Power Available (kW)")
            ax.plot(alt_array,prR_array,linewidth=2,label="Power Required (kW)")


        ax.set(xlabel='Altitude (meters)', ylabel='Power (kW)',title=strTitle)
        textstr = "$V_\infty = $" + str(Velocity) + "$\dfrac{m}{s}$" +"\nweight =" + str(round(weight / 1000,2)) + "kN"
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.85, 0.85, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top',horizontalalignment='center', bbox=props)
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)

        return fig

    #PART F
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
        ax.plot(alt_array,thrust_array,linewidth=2,label="Thrust Available (N)")
        ax.set(xlabel='Altitude (meters)', ylabel='Thrust available(N)',title='Altitude vs Thrust Available')
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)

        return fig
    
    def get_ROC_vel_alt(self,alt,vel,WEIGHT):
        powA = self.get_PowerAvailable_jet(alt,vel)
        powR = self.get_PowerRequired_alt_jet(alt,vel,WEIGHT)
        weight = self.weight_from_str(WEIGHT)
        powEx = powA - powR
        return (powEx/weight)
    
    def get_angle_climb_jet(self,Alt, WEIGHT):
        '''MAX ANGLE OF CLIMB
        '''
        weight = self.weight_from_str(WEIGHT)
        thA = self.get_ThrustAvailable_jet(Alt)
        k = calc_K_value(self.StatsCraft.mainwing.OswaldE,self.StatsCraft.mainwing.AR)
        sintheta = (thA / weight) - (1/calc_CL_CDmax(k,self.StatsCraft.Cd0))
        theta = np.arcsin(sintheta)
        return theta
    
    def graph_angle_max_ANGLE_OF_CLIMB(self, alt_lower,alt_upper,numPoints, WEIGHT):
        weight = self.weight_from_str(WEIGHT)
        alt_array = np.linspace(alt_lower,alt_upper,num=numPoints) #Array of numbers between lower and upper (inclusive)
        aoa_array = np.vectorize(self.get_angle_climb_jet)(alt_array,WEIGHT)
        aoad = np.rad2deg(aoa_array)
        fig, ax = plt.subplots()
        ax.plot(alt_array,aoad, label = "AOA")
        ax.set(xlabel='Altitude (meters)', ylabel='Angle of attack',title='Altitude vs maxAOA')
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)

        return fig

    def graph_ROC(self,alt_Lower,alt_Upper,numPoints,Velocity,WEIGHT,INFEETMIN: bool = False,SENDRAW: bool = False)-> plt:
        """Graphs rate of climb vs Altitude for the given velocity and weight.

        Args:
            alt_Lower (_type_): Lower limit for altitude
            Velocity (_type_): Velocity to be used in excess energy calculation (m/s)
            WEIGHT (_type_): Either "TAKEOFF" for the provided crafts takeoff weight, "EMPTY" for empty weight , "AVE" average the TOW and EW, or a float (newtons) (ie: 10.1). When a float is provided it will be cast from string to the float.
            INFEETMIN (bool, optional): Graph in ROC in feet/min rather than m/s Defaults to False.

        Returns:
            _type_: A plot of ROC vs Altitide
        """        
        weight = self.weight_from_str(WEIGHT)
    
        PowerCurve = self.graph_PowerAval_vs_PowerReq(alt_Lower,alt_Upper,numPoints,Velocity,str(weight),GRAPH_EXCESS=True,SENDRAW=True)
        """About  PowerCurve = [PowerCurve[0], (PowerCurve[1] * 1000)/weight] 
        this is deviding every number in powercurve[1] (our excess power) by the weight. 
        """
        if INFEETMIN:
            PowerCurve = [PowerCurve[0], ((PowerCurve[1] * 1000)/weight) * 3.28084 * 60] # Convert from m/s to ft/min
            if SENDRAW: #If we want, we can send back the data without graphing it.
                return PowerCurve
            
            fig, ax = plt.subplots()
            ax.plot(PowerCurve[0],PowerCurve[1],linewidth=2,label="Rate of Climb (ft / m)")
            ax.set(xlabel='Altitude (meters)', ylabel='Rate of Climb (ft/min)',title='Altitude vs RoC')
            ax.hlines(100,PowerCurve[0][0],PowerCurve[0][len(PowerCurve[0]) - 1],colors="red",linestyles="dotted",label="100 ft/min")
        else:
            PowerCurve = [PowerCurve[0], (PowerCurve[1] * 1000)/weight] #Convert to m/s
            if SENDRAW:
                return PowerCurve
            
            fig, ax = plt.subplots()
            ax.plot(PowerCurve[0],PowerCurve[1],linewidth=2,label="Rate of Climb (m/s)")
            ax.set(xlabel='Altitude (meters)', ylabel='Rate of Climb (m/s)',title='Altitude vs RoC')
            ax.hlines(0.508,PowerCurve[0][0],PowerCurve[0][len(PowerCurve[0]) - 1],colors="red",linestyles="dotted",label="0.508 m/s")
        textstr = "$V_\infty = $" + str(Velocity) + "$\dfrac{m}{s}$" +"\nweight =" + str(round(weight / 1000,2)) + "kN"
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.85, 0.85, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top',horizontalalignment='center', bbox=props)
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)

        return fig
    
    def graph_ROC_3d(self,alt_Lower,alt_Upper,numPoints,Velocity_min,Velocity_max,num_vel_points,WEIGHT,INFEETMIN: bool = False,SENDRAW: bool = False)-> plt:
        """Returns a 3d plot of RATE OF CLIMB over ALTUDUDE and VELOCITY

        Args:
            alt_Lower (_type_): lower alt limit
            alt_Upper (_type_): upper alt limit
            numPoints (_type_): number of altitude points
            Velocity_min (_type_): lower velocity limit
            Velocity_max (_type_): upper velocity limit
            num_vel_points (_type_): number of velocity points
            WEIGHT (_type_): WEIGHT (_type_): Either "TAKEOFF" for the provided crafts takeoff weight, "EMPTY" for empty weight , "AVE" average the TOW and EW, or a float (newtons) (ie: 10.1). When a float is provided it will be cast from string to the float.
            INFEETMIN (bool, optional): _description_. Defaults to False.
            SENDRAW (bool, optional): _description_. Defaults to False.

        Returns:
            plt: 3d plot
        """
        alt_array = np.linspace(alt_Lower,alt_Upper,num=numPoints) #y
        vel_array = np.linspace(Velocity_min,Velocity_max,num=num_vel_points)#x

        vectorized = np.vectorize(self.get_ROC_vel_alt) #Numpy magic

        X, Y = np.meshgrid(vel_array, alt_array)
        Z = vectorized(Y,X,WEIGHT)

        fig = plt.figure()
        ax = plt.axes(projection='3d',computed_zorder=False)
        cs = ax.contour(X,Y,Z, levels=[0.508],cmap=cm.summer,offset=+0,zorder = 5)
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,linewidth=0, antialiased=False,alpha = 1,zorder = 1)
        ax.set_xlabel("Velocity (m/s)")
        ax.set_ylabel("Altitude (m)")
        ax.set_zlabel("ROC (m/s)")
        ax.view_init(elev=30, azim=120, roll=0)

        textstr = "Green line: ROC = 100 ft/min"
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0, 0,350, textstr, transform=ax.transAxes, fontsize=10,verticalalignment='top',horizontalalignment='center', bbox=props)
        ax.set_box_aspect(aspect=None, zoom=0.8)
        fig.add_axes(ax)
        return fig

    def get_MAX_ROC_jet(self,alt,WEIGHT,RETURN_VEL: bool = False):
        '''Using max ROC for jets equation from lecture
        '''
        weight = self.weight_from_str(WEIGHT)
        wingarea = self.StatsCraft.mainwing.Area
        cd0 = self.StatsCraft.Cd0
        thrust = self.get_ThrustAvailable_jet(alt)
        k = calc_K_value(self.StatsCraft.mainwing.OswaldE,self.StatsCraft.mainwing.AR)
        ldmax = calc_CL_CDmax(k,cd0)
        dens = self.Active_Atmosphere.dens_trop_alt(alt)
        Z = 1 + np.sqrt(1 + (3/(np.power(ldmax,2) * np.power(thrust/weight,2))))

        if RETURN_VEL:
            V_rc = np.power((((thrust/weight) * (weight/wingarea))/(3 * dens * cd0)) * Z,0.5)
            return V_rc
        
        ST1 = np.power(((weight/wingarea)*Z)/(3 * dens * cd0 ),0.5)
        ST2 = np.power(thrust/weight,3/2)
        ST3 = 1 - (Z/6)- (3/(2 * np.power(thrust/weight,2) * np.power(ldmax,2) * Z))
        rcmax = ST1 * ST2 * ST3

        return rcmax

    def graph_MAX_ROC_JET(self,alt_lower,alt_upper,numPoints,WEIGHT,PLOT_CEILING: bool = False):
        weight = self.weight_from_str(WEIGHT)
        alt_array = np.linspace(alt_lower,alt_upper,num=numPoints) #Array of numbers between lower and upper (inclusive)
        ROC_array = np.vectorize(self.get_MAX_ROC_jet)(alt_array,WEIGHT)
        Vel_array = np.vectorize(self.get_MAX_ROC_jet)(alt_array,WEIGHT,RETURN_VEL = True)

        fig, ax = plt.subplots()
        ax.plot(alt_array,ROC_array, label = " Max Rate of climb")
        ax1 = ax.twinx()
        ax1.plot(alt_array,Vel_array, label = "Velocity",color="orange")
        ax1.set(xlabel='Altitude (meters)', ylabel='Velocity m/s')

        ax.set(xlabel='Altitude (meters)', ylabel='ROC (m/s)',title='Altitude vs Max RoC')
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)
        ax.hlines(0.508,alt_array[0],alt_array[len(alt_array) -1 ],colors="red",linestyles="dotted",label="0.508 m/s")
        
        if PLOT_CEILING:
            minarr = np.absolute(ROC_array - 0.508)
            locOfCelling = minarr.argmin()
            ax.vlines(alt_array[locOfCelling],ROC_array.min(),ROC_array.max(),colors="Purple",linestyles="dotted",label="Service Ceiling")
            textstr = "weight =" + str(round(weight / 1000,2)) + "kN\n" + "Service Ceiling: " + str(round(alt_array[locOfCelling],1)) + "m"
        else:
            textstr = "weight =" + str(round(weight / 1000,2)) + "kN"
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.5, 0.85, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top',horizontalalignment='center', bbox=props)
        return fig

       

if __name__ == "__main__":
    OppaStoppa = Craft("OppaStoppa")
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
    OppaStoppa.mainwing = MainWing

    OppaStoppa.compute_components()


    MyCraftStats = CraftStatistics(OppaStoppa)

    MAXRoc = MyCraftStats.graph_MAX_ROC_JET(0,12000,1000,"AVE",PLOT_CEILING=True)
    MAXRoc.legend()

    plt.show()
    