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
        """CraftStatistics serves as a layer to make graphing properties of the craft simple.

        Args:
            Craft (Craft): _description_
        """        
        self.StatsCraft = Craft
        self.Active_Atmosphere = Craft.Atmosphere
        self.ur = self.StatsCraft.ur

        #print("Initalizing Statistics for Craft: " + self.StatsCraft.name)
        #print("CD0: " + str(self.StatsCraft.Cd0))
        #print("Ozwald: " + str(self.StatsCraft.mainwing.OswaldE))

    #PART F
    def weight_from_str(self, WEIGHT:str)-> float:
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
                print("Conversion to float failed: " + WEIGHT)
        return weight

    #PART F
    def get_ThrustAvailable_jet(self,alt,vel=0):
        """Taking into consideration all sources of thrust on the craft, finds the thrust available at the given altitude.
        Assumes m = 1

        Args:
            alt (_type_): Altitude above sea level in meters
        """        
        maxthurust = self.StatsCraft.get_max_thrust()
        aval_thrust = self.Active_Atmosphere.jet_thrust_available_alt(maxthurust,alt,1)
        return aval_thrust
    
    def get_ThrustAvailable_Elec(self,alt,vel,Mrpm =0):
        motor:ElectricMotor = self.StatsCraft.powertrain[0] #gets one of our motors
        prop:Propeller = motor.prop # get the prop on the motor 
        propDia = prop.diameter
        if Mrpm != 0:
            rpm = Mrpm * self.ur.rpm
        else:
            rpm = motor.maxRPM #we want max power avail so max RPM
        RPS = rpm.to("rad/sec")
        advr = prop.calc_advance_ratio(vel,RPS) #determine advance ratio
        dens = self.Active_Atmosphere.dens_trop_alt(alt)
        CT = prop.thrust_polynomal_func(advr)
        Thrust = CT * dens * numpy.power(RPS /(2 * np.pi),2) * numpy.power(propDia, 4)

        Thrust = Thrust * len(self.StatsCraft.powertrain)
        return Thrust

    
    def get_thrustAvailable(self,alt,vel,Mrpm = 0):
        """Gets thrust available reguardless of type

        Args:
            alt (_type_): alt
            vel (_type_): velocity
            Mrpm (int, optional): Electric motor RPM. Defaults to 0, meaning max rpm

        Returns:
            _type_: _description_
        """        
        #print(type(self.StatsCraft.powertrain[0]).__name__)
        if type(self.StatsCraft.powertrain[0]).__name__ == "ElectricMotor":
            thr =  self.get_ThrustAvailable_Elec(alt,vel,Mrpm)
        else:
            thr =  self.get_ThrustAvailable_jet(alt,vel)


        return (thr)

    def get_Power(self,alt,vel,Mrpm = 0):
        """Gets power available for a craft reguardless of type

        Args:
            alt (_type_): altitude
            vel (_type_): velocity
            Mrpm (int, optional): Electric motor RPM. Defaults to 0, meaning max RPM
        """        
        #print(type(self.StatsCraft.powertrain[0]).__name__)
        if type(self.StatsCraft.powertrain[0]).__name__ == "ElectricMotor":
            pow = self.get_PowerAvailable_Elec(alt,vel,Mrpm)
        else:
            pow = self.get_PowerAvailable_jet(alt,vel)

        return(pow)

    #PART F
    def get_PowerAvailable_jet(self,alt,velocity):
        """Returns the power available from a jet at given alt

        Args:
            alt (_type_): Altiude 
            velocity (_type_): Velocity of incoming air (m/s)

        Returns:
            _type_: Power available (watt)
        """        
        return ((self.get_thrustAvailable(alt,velocity) * velocity).magnitude * self.ur.watt)
    
    def get_PowerAvailable_Elec(self,alt,vel,Mrpm = 0):
        """Gets the power avalible from Electric motors

        Returns:
            _type_: Power avalible in watts
        """   
        motor:ElectricMotor = self.StatsCraft.powertrain[0] #gets one of our motors
        prop:Propeller = motor.prop # get the prop on the motor
        if Mrpm != 0:
            rpm = Mrpm * self.ur.rpm
        else: 
            rpm = motor.maxRPM #we want max power avail so max RPM
        RPS = rpm.to("rad/sec")
        advr = prop.calc_advance_ratio(vel,RPS) #determine advance ratio
        CT = prop.thrust_polynomal_func(advr)
        CP = prop.power_polynomal_func(advr) #CHANGE THIS MAYBE
        dens = self.Active_Atmosphere.dens_trop_alt(alt)
        Eff = (CT*advr)/CP
        #print(f"Effic: {Eff}")

        PA = (CP * dens * numpy.power(RPS / (2 * numpy.pi),3) * numpy.power(prop.diameter, 5)) * len(self.StatsCraft.powertrain) #calc power from coeff of power
        #PA = self.get_ThrustAvailable_Elec(alt,vel) * vel # calc power from Thrust * vel

        #for motor in self.StatsCraft.powertrain: #total up all power availible  # calc from effic
        #    PA += motor.MaxPower
        #
        #PA = PA * motor.effic # take into account motor effic
        #PA = PA * Eff # take into account Prop effic
        return PA
    
    #PART F
    def get_PowerRequired_alt_jet(self,alt,velocity,WEIGHT):
        """Finds power requred for SLF at given conditions

        Args:
            alt (_type_): Altitude above sealevel
            velocity (_type_): Velocity in m/s
            WEIGHT (str): Either "TAKEOFF" for the provided crafts takeoff weight, "EMPTY" for empty weight , "AVE" average the TOW and EW, or a float (ie: 10.1). When a float is provided it will be cast from string to the float.

        Returns:
            Pr: Power requred 
        """        
        if type(WEIGHT) == str:
            weight = self.weight_from_str(WEIGHT)
        else:
            weight = WEIGHT
    

        if(self.ur.get_dimensionality(alt) != self.ur.get_dimensionality( 1 * self.ur.meters)):
            alt = alt * self.ur.meters
        if self.ur.get_dimensionality(velocity) != self.ur.get_dimensionality(1 * self.ur.m / self.ur.s):
            velocity = velocity * self.ur.m / self.ur.seconds
        

        dens = self.Active_Atmosphere.dens_trop_alt(alt)
        k = calc_K_value(self.StatsCraft.mainwing.OswaldE,self.StatsCraft.mainwing.AR)
        #print(f'w: {weight}, dens: {dens}, K: {k}, alt: {alt}, vel: {velocity}')
        Pr = calc_PowerReq(dens,velocity,self.StatsCraft.mainwing.Area,self.StatsCraft.Cd0,k,weight)
        return Pr
    
    #PART F
    def graph_PowerAval_vs_PowerReq(self,Alt_lower, Alt_upper,numPoints, Velocity, WEIGHT,GRAPH_EXCESS: bool = False,SENDRAW: bool = False):
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

        try:  #catch if velocity is provided without units
            print(self.ur.get_dimensionality(Velocity))
        except:
            Velocity = Velocity * self.ur.m / self.ur.s
        
        alt_array = np.linspace(Alt_lower,Alt_upper,num=numPoints) #Array of numbers between lower and upper (inclusive)
        prA_array = np.zeros(alt_array.shape)
        prR_array = np.zeros(alt_array.shape)

        for i in enumerate(alt_array): #Iterate over altitude and calculate Thrust
            #print(self.get_PowerAvailable_jet(i[1],Velocity*self.ur.m/self.ur.s))
            prA_array[i[0]]= self.get_Power(i[1],Velocity).magnitude /1000 #need to strip units :(
            prR_array[i[0]]= self.get_PowerRequired_alt_jet(i[1],Velocity,WEIGHT).magnitude /1000

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
        textstr = "$V_\infty = $" + str(Velocity) + "$\dfrac{m}{s}$" +"\nweight =" + str(round(weight.to("kilonewton"),2))
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
            thrust_array[i[0]] = self.get_thrustAvailable(i[1],self.Active_Atmosphere.Vinfinity).magnitude
            dens_array[i[0]] = self.Active_Atmosphere.dens_trop_alt(i[1]).magnitude
            temp_array[i[0]] = self.Active_Atmosphere.temp_trop(i[1]).magnitude

        fig, ax = plt.subplots()
        ax.plot(alt_array,thrust_array,linewidth=2,label="Thrust Available (N)")
        ax.set(xlabel='Altitude (meters)', ylabel='Thrust available(N)',title='Altitude vs Thrust Available')
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)

        return fig
    
    #PART F
    def get_ROC_vel_alt(self,alt,vel,WEIGHT):
        """Gets the achivable ROC for a jet

        Args:
            alt (_type_): Altitude
            vel (_type_): Velocity m/s
            WEIGHT (_type_): Weight newtons

        Returns:
            _type_: Rate of climb (m/s)
        """        
        weight = self.weight_from_str(WEIGHT)
        powA = self.get_Power(alt,vel)
        powR = self.get_PowerRequired_alt_jet(alt,vel,weight)
        #print(f'Poweravil: {powA}')
        #print(f'power R: {powR}')
        powEx = powA - powR
        #print(powEx)
        return ((powEx/weight).magnitude) #had to strip units for vecorization

    #PART F
    def get_angle_climb_jet(self,Alt, WEIGHT):
        """Returns the maximum angle of climb for a jet powered craft

        Args:
            Alt (_type_): Altitude
            WEIGHT (_type_): Weight (newton)

        Returns:
            _type_: max angle of climb (rad)
        """        
        weight = self.weight_from_str(WEIGHT)
        thA = self.get_thrustAvailable(Alt,self.Active_Atmosphere.Vinfinity)
        k = calc_K_value(self.StatsCraft.mainwing.OswaldE,self.StatsCraft.mainwing.AR)
        sintheta = (thA / weight) - (1/calc_CL_CDmax(k,self.StatsCraft.Cd0))
        theta = np.arcsin(sintheta)
        return theta

    #PART F
    def graph_angle_max_ANGLE_OF_CLIMB(self, alt_lower,alt_upper,numPoints, WEIGHT):
        """Gets the max angle of constant climb

        Args:
            alt_lower (_type_): _description_
            alt_upper (_type_): _description_
            numPoints (_type_): _description_
            WEIGHT (_type_): _description_

        Returns:
            _type_: _description_
        """        
        weight = self.weight_from_str(WEIGHT)
        alt_array = np.linspace(alt_lower,alt_upper,num=numPoints) #Array of numbers between lower and upper (inclusive)
        aoa_array = np.vectorize(self.get_angle_climb_jet)(alt_array,WEIGHT)
        aoad = (np.rad2deg(aoa_array))*self.ur.degrees
        fig, ax = plt.subplots()
        ax.plot(alt_array,aoad, label = "AOA")
        ax.set(xlabel='Altitude (meters)', ylabel='Angle of attack',title='Altitude vs maxAOA')
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)

        return fig

    #PART F
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
        #Velocity = Velocity * ur.m /ur.s
    
        PowerCurve = self.graph_PowerAval_vs_PowerReq(alt_Lower,alt_Upper,numPoints,Velocity,WEIGHT,GRAPH_EXCESS=True,SENDRAW=True)
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
        textstr = "$V_\infty = $" + str(Velocity) + "$\dfrac{m}{s}$" +"\nweight =" + str(round(weight.to("kilonewton"),2))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.85, 0.85, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top',horizontalalignment='center', bbox=props)
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)

        return fig

    #PART F
    def graph_ROC_3d(self,alt_Lower,alt_Upper,numPoints,Velocity_min,Velocity_max,num_vel_points,WEIGHT,INFEETMIN: bool = False,SENDRAW: bool = False)-> plt: # working post units
        """Returns a 3d plot of RATE OF CLIMB over ALTUDUDE and VELOCITY

        Args:
            alt_Lower (_type_): lower alt limit
            alt_Upper (_type_): upper alt limit
            numPoints (_type_): number of altitude points
            Velocity_min (_type_): lower velocity limit
            Velocity_max (_type_): upper velocity limit
            num_vel_points (_type_): number of velocity points
            WEIGHT (_type_): Either "TAKEOFF" for the provided crafts takeoff weight, "EMPTY" for empty weight , "AVE" average the TOW and EW, or a float (newtons) (ie: 10.1). When a float is provided it will be cast from string to the float.
            INFEETMIN (bool, optional): _description_. Defaults to False.
            SENDRAW (bool, optional): _description_. Defaults to False.

        Returns:
            plt: 3d plot
        """
        alt_array = np.linspace(alt_Lower,alt_Upper,num=numPoints) #y
        vel_array = np.linspace(Velocity_min,Velocity_max,num=num_vel_points)#x

        vectorized = np.vectorize(self.get_ROC_vel_alt) #Numpy magic

        X, Y = np.meshgrid(vel_array, alt_array)
        #print(f'X: {X}, Y: {Y}')
        Z = (vectorized(Y,X,WEIGHT))

        fig = plt.figure()
        ax = plt.axes(projection='3d',computed_zorder=False)
        cs = ax.contour(X,Y,Z, levels=[0.508],cmap=cm.summer,offset=+0,zorder = 5)
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,linewidth=0, antialiased=False,alpha = 1,zorder = 1)
        ax.set_xlabel("Velocity (m/s)")
        ax.set_ylabel("Altitude (m)")
        ax.set_zlabel("ROC (m/s)")
        ax.view_init(elev=30, azim=120)

        textstr = "Green line: ROC = 100 ft/min"
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0, 0,350, textstr, transform=ax.transAxes, fontsize=10,verticalalignment='top',horizontalalignment='center', bbox=props)
        ax.set_box_aspect(aspect=None, zoom=0.8)
        fig.add_axes(ax)
        return fig

    #PART F
    def get_MAX_ROC_jet(self,alt,WEIGHT,RETURN_VEL: bool = False): #working post units
        """Returns Max ROC for a jet

        Args:
            alt (_type_): Altitude (meters)
            WEIGHT (_type_): weight (newtons)
            RETURN_VEL (bool, optional): Return velcity of Max ROC rather than ROC. Defaults to False.

        Returns:
            _type_: Rate of climb (m/s)
        """        

        if(self.ur.get_dimensionality(alt) != self.ur.get_dimensionality( 1 * self.ur.meters)):
            #print("oof")
            alt = alt * self.ur.meters        
        weight = self.weight_from_str(WEIGHT)
        wingarea = self.StatsCraft.mainwing.Area
        cd0 = self.StatsCraft.Cd0
        thrust = self.get_thrustAvailable(alt,self.Active_Atmosphere.Vinfinity)
        k = calc_K_value(self.StatsCraft.mainwing.OswaldE,self.StatsCraft.mainwing.AR)
        ldmax = calc_CL_CDmax(k,cd0)
        dens = self.Active_Atmosphere.dens_trop_alt(alt)
        #print(f'GetMaxRocJet: \nweight: {weight} \nwingarea: {wingarea}\n thrust: {thrust}\n dens: {dens}\n Alt: {alt}')
        Z = 1 + np.sqrt(1 + (3/(np.power(ldmax,2) * np.power(thrust/weight,2))))

        if RETURN_VEL:
            V_rc = np.power((((thrust/weight) * (weight/wingarea))/(3 * dens * cd0)) * Z,0.5)
            return V_rc.magnitude
        
        ST1 = np.power(((weight/wingarea)*Z)/(3 * dens * cd0 ),0.5)
        ST2 = np.power(thrust/weight,3/2)
        ST3 = 1 - (Z/6)- (3/(2 * np.power(thrust/weight,2) * np.power(ldmax,2) * Z))
        rcmax = ST1 * ST2 * ST3
        #print("RCMAX: " + str(rcmax))
        return rcmax.magnitude

    #PART F
    def graph_MAX_ROC_JET(self,alt_lower,alt_upper,numPoints,WEIGHT,PLOT_CEILING: bool = False):
        """Returns a plot of absolute Max rate of climb over the given altitude range, Also plots the Velocity requred to achive that RoC

        Args:
            alt_lower (_type_): Lower alt limit (meters)
            alt_upper (_type_): Upper alt limit (meters)
            numPoints (_type_): number of graph points
            WEIGHT (_type_): Either "TAKEOFF" for the provided crafts takeoff weight, "EMPTY" for empty weight , "AVE" average the TOW and EW, or a float (newtons) (ie: 10.1). When a float is provided it will be cast from string to the float.
            PLOT_CEILING (bool, optional): Plot where the ROC equals 100 ft/min. Defaults to False.

        Returns:
            _type_: _description_
        """        
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
    
    def graph_prop_effic_advr(self,vmin,vmax,Motorrpm = 0):
        """Graphs prop effic over velocity

        Args:
            vmin (_type_): v min
            vmax (_type_): v max
            Motorrpm (int, optional): RPM to use in math. Defaults to 0, meaning Motors max rpm

        Returns:
            _type_: _description_
        """        
        motor:ElectricMotor = self.StatsCraft.powertrain[0]
        prop= motor.prop
        velarr = np.linspace(vmin,vmax,(vmax - vmin))
        advrarr = np.zeros(velarr.shape)
        effarr = np.zeros(velarr.shape)
        powerarr = np.zeros(velarr.shape)
        Thrustarr = np.zeros(velarr.shape)
        if Motorrpm != 0:
            rpm = Motorrpm * self.ur.rpm
        else:
            rpm = motor.maxRPM #we want max power avail so max RPM
        RPS = rpm.to("rad/sec")
        alt = self.Active_Atmosphere.Altitude

        idx = 0
        for i in velarr:
            print(f'calling adv with {RPS} and {i}')
            advrarr[idx] = prop.calc_advance_ratio(i*self.ur.m / self.ur.s,RPS)
            powerarr[idx] = self.get_Power(alt,(i*self.ur.m / self.ur.s),Motorrpm).magnitude
            Thrustarr[idx] = self.get_thrustAvailable(alt,i*self.ur.m / self.ur.s,Motorrpm).magnitude
            effarr[idx] = prop.get_effic_from_advr(advrarr[idx])
            print(f'Got {effarr[idx]}')
            idx += 1

        fig, ax = plt.subplots()
        ax.plot(velarr,powerarr / 100,label = "Power hW")
        ax.plot(velarr,Thrustarr,label = "Thrust N")
        
        ax2 = ax.twinx()
        ax2.plot(velarr,advrarr, label = "adv ratio", color="grey", linestyle= "--")
        ax2.plot(velarr,effarr, label = "Effic", color = "green", linestyle = "--")

        ax2.axes.set_ylabel("Ratio")
        ax2.legend(loc="upper right")
        ax.legend(loc="upper left")
        ax.axes.set_xlabel("Velocity m/s")
        ax.set_title("Prop Efficiency Model")
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)

        textstr = f'Prop RPM: {rpm}'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        fig.text(0.5, 0.2,textstr, transform=ax.transAxes, fontsize=10,verticalalignment='top',horizontalalignment='center', bbox=props)
        #print(effarr)
        return fig
    
    def graph_prop_effic_advr_alt(self,altmin,altmax, Motorrpm = 0):
        """Graphs prop efficency over altitude

        Args:
            altmin (_type_): min alt
            altmax (_type_): max alt
            Motorrpm (int, optional): RPM to use in calculation. Defaults to 0, meaning the motors max RPM

        Returns:
            _type_: _description_
        """        
        motor:ElectricMotor = self.StatsCraft.powertrain[0]
        prop= motor.prop
        altarr = np.linspace(altmin,altmax)
        advrarr = np.zeros(altarr.shape)
        effarr = np.zeros(altarr.shape)
        powerarr = np.zeros(altarr.shape)
        Thrustarr = np.zeros(altarr.shape)
        if Motorrpm != 0:
            rpm = Motorrpm * self.ur.rpm
        else:
            rpm = motor.maxRPM #we want max power avail so max RPM
        RPS = rpm.to("rad/sec")
        vel = self.Active_Atmosphere.Vinfinity

        idx = 0
        for i in altarr:
            advrarr[idx] = prop.calc_advance_ratio(vel,RPS)
            effarr[idx] = prop.get_effic_from_advr(advrarr[idx])
            powerarr[idx] = self.get_Power(i,vel).magnitude
            Thrustarr[idx] = self.get_thrustAvailable(i,vel).magnitude
            idx += 1

        fig, ax = plt.subplots()
        ax.plot(altarr,powerarr / 100,label = "Power hW")
        ax.plot(altarr,Thrustarr,label = "Thrust N")
        
        ax2 = ax.twinx()
        ax2.plot(altarr,advrarr, label = "adv ratio", color="grey", linestyle= "--")
        ax2.plot(altarr,effarr, label = "Effic", color = "green", linestyle = "--")

        ax2.axes.set_ylabel("Ratio")
        ax2.legend(loc="upper right")
        ax.legend(loc="upper left")
        ax.axes.set_xlabel("Alt")
        ax.set_title("Prop Efficiency Model")
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)
        textstr = f'Prop RPM: {rpm}'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        fig.text(0.5, 0.2,textstr, transform=ax.transAxes, fontsize=10,verticalalignment='top',horizontalalignment='center', bbox=props)
        #print(effarr)
        return fig
   
    def graph_MAX_ROC_PROP(self,alt_lower,alt_upper,numpoints,WEIGHT):
        weight = self.weight_from_str(WEIGHT)
        motor:ElectricMotor = self.StatsCraft.powertrain[0]
        prop:Propeller = motor.prop 
        wing = self.StatsCraft.mainwing
        k = calc_K_value(wing.OswaldE,wing.AR)
        LDmax = calc_CL_CDmax(k,self.StatsCraft.Cd0)
        alt_lower = alt_lower * self.ur.m 
        alt_upper = alt_upper * self.ur.m

        altarr = np.linspace(alt_lower,alt_upper,numpoints)
        velarr = np.zeros(altarr.shape)
        rcarr = np.zeros(altarr.shape)

        idx = 0
        for i in altarr:
            dens = self.Active_Atmosphere.dens_trop_alt(i)
            vel = np.sqrt((2/dens) * np.sqrt(k/(3*self.StatsCraft.Cd0)) * (weight/wing.Area))
            velarr[idx] = vel.magnitude
            adv = prop.calc_advance_ratio(vel,motor.maxRPM.to("rad/sec"))
            propeff = prop.get_effic_from_advr(adv)
            #print(f'Adv: {adv}, eff: {prop.get_effic_from_advr(adv)}, Rad/s: {motor.maxRPM.to("rad/sec")}')
            #rc = (( motor.MaxPower * motor.effic * propeff)/weight) - vel * (1.155/LDmax)
            rc = (( self.get_PowerAvailable_Elec(i,vel))/weight) - vel * (1.155/LDmax)

            rcarr[idx] = rc.magnitude 
            idx +=1

        fig, ax = plt.subplots()
        fig.text(0.5, 0.95, self.StatsCraft.name, horizontalalignment="center",fontsize = 10)
        ax2 = ax.twinx()
        ax2.plot(altarr,velarr, label = "Velocity",linestyle = "--",color = "grey")
        ax2.axes.yaxis.set_label("Velocity (m/s)")
        ax.plot(altarr,rcarr,label="Rate of climb")
        ax.hlines(0.508,alt_lower,alt_upper,colors=["red"],linestyles=['dotted'],label="0.508 m/s (100 ft/min)")
        minarr = np.absolute(rcarr - 0.508)
        locOfCelling = minarr.argmin()
        ax.vlines(altarr[locOfCelling],rcarr.min(),rcarr.max(),colors="Purple",linestyles="dotted",label="Service Ceiling")
        textstr = "Service Ceiling: " + str(round(altarr[locOfCelling],1))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.5, 0.85, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top',horizontalalignment='center', bbox=props)
        ax.axes.set_ylabel("Meters per Second")
        ax.axes.set_xlabel("Altitude (Meters)")
        ax.set_title("MAX Rate of Climb vs Altitude")


        return fig




       

if __name__ == "__main__":
        
        #TESTING OF CLASS, NOT ACTUALL CRAFT PROPTERIES 
        OppaStoppa = Craft("OppaStoppa")
        ur = OppaStoppa.ur
        OppaStoppa.Atmosphere = Atmosphere(300,286.21, 9.77774,1.19,15,ur)
        atmo = OppaStoppa.Atmosphere
        OppaStoppa.weight_empty = 24 * 9.81 * ur.newtons
        OppaStoppa.weight_takeoff = 24 * 9.81 * ur.newtons
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

        MyCraftStats = CraftStatistics(OppaStoppa)

        MAXRoc = MyCraftStats.graph_MAX_ROC_PROP(0,18000,300,"AVE")
        MAXRoc.legend()



        #powavr = MyCraftStats.graph_PowerAval_vs_PowerReq(0,12000,300,25,"AVE",GRAPH_EXCESS=True)
        #powavr.legend()

        #ROC3d = MyCraftStats.graph_ROC_3d(0,12000,50,15,50,50,"AVE")

        #ROC = MyCraftStats.graph_ROC(0,18000,100,33,"AVE")
        #ROC.legend()

        #AOC = MyCraftStats.graph_angle_max_ANGLE_OF_CLIMB(0,12000,200,"AVE")

        #THRA = MyCraftStats.graph_ThrustAvailableS(0,12000,100)

        #test = MyCraftStats.graph_prop_effic_advr(0,200)
        #test2 = MyCraftStats.graph_prop_effic_advr_alt(0,18000,6321)




        #plt.legend()
        plt.show()
        
        
