import numpy as np
import matplotlib.pyplot as plt
import pint 
import TakeOffLanding
import Atmosphere
from LoadFactor import LoadFactor
from craft import Craft
from StFl import *
from Components import *
from CraftStatistics import CraftStatistics


#Kinda assembling the Avengers with these import statements
#This var will control if you want to display all the graphs or not.
SHOWGRAPHS = True #Warning, LOTS of graphs. Like damn near every one we made.
ROUGH_PROP_EFFIC = 0.8 #VERY rough estimate of prop effic, Code includes a more complicatesd Prop effic model, 0.8 is about where our props max out.


#First we need to define our craft.
OppaStoppa = Craft("OppaStoppa") #Empty Craft, with just a name.
ur = OppaStoppa.ur # Getting the unit reg thats created with the craft to use.

#Giving the craft an atmosphere to use for calculations involving atmosphereic quantites
#In the future, when a function does not ask for atmosphereic conditions, it will use the atmospheres defaults assingned here.
OppaStoppa.Atmosphere = Atmosphere(600,286.21, 9.77774,1.19,25,ur) 
atmo = OppaStoppa.Atmosphere
#Defining the crafts weights, becuse we are electric, they are the same
OppaStoppa.weight_empty = 15 * 9.81 * ur.newtons
OppaStoppa.weight_takeoff = 15 * 9.81 * ur.newtons
#Data from the wings CL vs alpha graph
OppaStoppa.CLmax = 1.45
OppaStoppa.Clmin = 0.9
OppaStoppa.CLrolling = 0.43 #assuming AOA of 2.5 degrees

#Setting load maxes
OppaStoppa.PosStruturalNlimit = 5
OppaStoppa.NegStruturalNlimit = 3


"""Defining the draggy components of our craft, each component has methods to calculate its own contribution to the crafts CD0"""
#Wing defined: NAME, AIRFOIL, SWEEP, AREA, SPAN, CHORD, angleZeroLift, AngleStall, TC, XC, AreaWIngObscured, atmosphere
N4312Wing = Wing3d("Oppa Main Wing","NACA 4312",      34.83,1.29,    2.34,   .532,-1,15,0.12,0.3,0.02,atmo)
MainWing = N4312Wing
HorizontalTail = Wing3d("Horizontal Tail","NACA 0012",26.56,.206 * 2,.73 * 2,.274,0 ,15,0.12,0.3,0.008,atmo)
VerticalTail = Wing3d("Vertical Tail","NACA 0012",    26.56,.206,    .73,    .274,0 ,15,0.12,0.3,0.004,atmo)

#Fuselage defined: NAME, Length, AreaTop, AreaSide, maxCrossSectionArea, Interf, MainwingArea, atmosphere
MainFuselage = Fuselage("Oppa Fuselage",1.926,.48,.39,.019,1,MainWing.Area,atmo)

#Gear defined: NAME, CD_component, FrontalArea, MainwingArea, Interf, atmosphere
Gear = FixedGear("Tail Gear",0.25,0.008,MainWing.Area,1.2,atmo) #defined here but not included in crafts dragcompoennts 

#MotorNacelle = Nacelle("MotorNacelle", length, maxCrossAra,interf, AreaWet, Mainwing.Area,atmo)

"""Defining our engines, For the majority of the project we where using Turbojet engines,
Only in the last few weeks did we decide to change to Electric prop, so the implementation of props is somewhat rushed."""

prop = Propeller("16X8",0.4064,8,ur) #Define a basic propeller, 0.4m in dia, pitch of 8 
def thrustproppoly(advr): # To calcualate thrust available, we are using a best fit polynomial from test data of the props
    try:
        advr  = advr.magnitude
    except:
        advr = advr
    thr = 0.122 - (0.0138 * advr) + (0.0709* advr**2) - (0.287 * advr**3) + (0.137 * advr**4)
    return thr
prop.thrust_polynomal_func = thrustproppoly

def powerproppoly(advr): # To calcualate power available, we are using a best fit polynomial from test data of the props
    try:
        advr  = advr.magnitude
    except:
        advr = advr
    pwr = 0.0583 + (0.0127 * advr) + (0.23* advr**2) - (0.416 * advr**3) + (0.151 * advr**4)
    return pwr
prop.power_polynomal_func = powerproppoly

Motor = ElectricMotor("V804 KV170",5200,0.90,8000,prop,ur) # Defining our electric motors. Maxpower, Effic, MAXRPM, prop

Battery = Battery(139.76,4.32 * 2,45,ur) #Defining our battery pack. Density, Mass, voltage
#Battery.print_state() #Print our batteries state, to validate the parameters are correct.

OppaStoppa.powertrain = [Motor,Motor] #Adding the motors to the craft, when the craft then neets to calcualte power or thrust available, it will add both motors together 
OppaStoppa.mainwing = MainWing #Give the aircraft a refrence to its mainwing, so it can access it area other properties 

OppaStoppa.dragcomponents = [MainWing,MainFuselage,HorizontalTail,VerticalTail] #List of all of the compoenents that will add to CD0

OppaStoppa.compute_components() #Compute the CD0

OppaStoppa.Atmosphere.printAtmos() #Prints the details of the atmosphere
print('##Begin C_D0 Calculations')
OppaStoppa.print_stats_components() #Prints Details About overall C_d0 Calculation
print("Total CD0: " + str(OppaStoppa.Cd0))#Print TOtal CD0
print('##End C_D0 Calculations\n')


k = calc_K_value(MainWing.OswaldE,MainWing.AR)# K calc for mainwing
q = calc_dynpressure(OppaStoppa.Atmosphere.Density,OppaStoppa.Atmosphere.Vinfinity) #Dynamic pressure at cruse.
print("dynamic pressure" + str(q) + " K: " + str(k))
#print(f'Vstall calc inputs: Density: {OppaStoppa.Atmosphere.Density} Weight: {OppaStoppa.weight_takeoff} WingArea: {MainWing.Area}')
vstall = calc_Vstall(OppaStoppa.Atmosphere.Density,OppaStoppa.weight_takeoff,MainWing.Area,OppaStoppa.CLmax) #Calcualte the stall speed

#Find various cl/cd ratios
clcdmax = calc_CL_CDmax(k,OppaStoppa.Cd0)
clcd2max = calc_CL_CDmaxRangeJet(k,OppaStoppa.Cd0)
cl23cd = calc_CL_CDmaxenduranceProp(k,OppaStoppa.Cd0)
print("CL/CD values:")
print(" CL/CD max: " + str(clcdmax))
print(" CL^0.5/CD max :" + str(clcd2max))
print(" CL^3/2/CD max :" + str(cl23cd))

#Calcualte endurance assumnng a prop effic of 0.8
maxrendurance = calc_endurance_electric_prop(Battery.MaxEnergy,Motor.effic,ROUGH_PROP_EFFIC,atmo.Density,MainWing.Area,OppaStoppa.weight_takeoff,cl23cd)
#Calcualte range assiuming a prop effic of 0.8
maxrange = calc_range_electric_prop(Battery.MaxEnergy,Motor.effic,ROUGH_PROP_EFFIC,OppaStoppa.weight_takeoff,clcdmax)
#calculate Velocity of Max range
vinfMR = calc_Vinf_MaxJetEndurance(atmo.Density,OppaStoppa.weight_takeoff,MainWing.Area,OppaStoppa.Cd0,k)
#calcualte Vlocity of max endurance
vinfME = calc_Vinf_MaxPropEndurance(atmo.Density,OppaStoppa.weight_takeoff,MainWing.Area,OppaStoppa.Cd0,k)

#Output various stats about Level flight at cruise conditions
Powerreq = calc_PowerReq(OppaStoppa.Atmosphere.Density,OppaStoppa.Atmosphere.Vinfinity,MainWing.Area,OppaStoppa.Cd0,k,OppaStoppa.weight_takeoff).to("watt")
print("\n##Cruise Condition Statistics")
print("| Power Req: " + str(Powerreq))
print("| Thrust req: " + str(calc_ThrustReq(q,MainWing.Area,OppaStoppa.Cd0,k,OppaStoppa.weight_takeoff).to("newton")))
print("| V Stall: " + str(vstall.to("meter/second")))
estCL = calc_req_CL(OppaStoppa.weight_takeoff,q,MainWing.Area)
print("| Estimated CL for routine flight: " + str(estCL.magnitude))
print("| Max Endurance:" + str(maxrendurance.to("hours")) )
print("|  -at Vinf: " + str(vinfME.to("meters/second")))
print("| Max Range:" + str(maxrange.to("kilometers")))
print("|  -at Vinf: " + str(vinfMR.to("meters/second")))
print("##End Cruise Condition Statisitcs \n")

print("#BEGINNING FLIGHT")
Battery.print_state()
MyStats = CraftStatistics(OppaStoppa)

#Begin takeoff section
print(f'\n##Begin Takeoff Section')
TakeoffClimbAngle = 6
MyTakeoff = TakeOffLanding.Takeoff(OppaStoppa,300,1.15,TakeoffClimbAngle) #Create a takeoff
print ("Total Takeoff Distance: " + str(MyTakeoff.DoTakeoff())) #Perform a takeoff and output stats
PlanesRealDist = MyTakeoff.get_DistTraveled() # [Roll Dist, airDIst]
#print(f'Plane traveled a total distance of {PlanesRealDist}')

#Full throttle for GroundRoll section of take off, assuming a linear acceleration the average speed would be 1.2 vstall
powerCons = Motor.MaxPower *2 #both motors at full power (This is consupmtion not output)
avevel = MyTakeoff.V_LO / 2
timeRoll = PlanesRealDist[0] / avevel
EnergyConsumedTORoll = (powerCons * timeRoll).to("J") #done calculating Energy used in roll

#Power requred for SLF
PowerReqTOclimb = calc_PowerReq(atmo.dens_trop_alt(MyTakeoff.Alt),atmo.Vinfinity,MainWing.Area,OppaStoppa.Cd0,k,OppaStoppa.weight_takeoff).to("watt")
thrReqTOclimb = (PowerReqTOclimb / atmo.Vinfinity).to('newton')
threxcessReqTOclimb = OppaStoppa.weight_takeoff * numpy.sin(numpy.deg2rad(MyTakeoff.climbAngle))
#print(f'Pow req: {PowerReqTOclimb}, Thrust for SLF after TO: {thrReqTOclimb}, Excess to climb: {threxcessReqTOclimb}')
#print(f"Thrust Requred for post liftoff climb: {threxcessReqTOclimb + thrReqTOclimb}, @({MyTakeoff.climbAngle * ur.deg}, {atmo.Vinfinity})")
PowerClimb = ((threxcessReqTOclimb + thrReqTOclimb) * atmo.Vinfinity) / (Motor.effic * ROUGH_PROP_EFFIC) #Rough power to climb assiming props are 80% efficient 
climbTime = PlanesRealDist[1] / atmo.Vinfinity #Rough time of climb
EnergyConsumedTOClimb = (PowerClimb * climbTime).to("J") #Done calculating energy to climb
Battery.discharge((EnergyConsumedTOClimb + EnergyConsumedTORoll)) #Consuming the power for takeoff

print(f'\n #Energy To Takeoff#')#Print out the power stats
print(f'  | Power for Takeoff roll: {powerCons}, Rolling for {timeRoll}. Roll will use {EnergyConsumedTORoll}.')
print(f'  | Power to maintain climb: {PowerClimb.to("watt")}, Climbing for {climbTime}. Climb will use {EnergyConsumedTOClimb}. ')
Battery.print_state() #Checkbattery stats
if SHOWGRAPHS:
    to = MyTakeoff.Graph() #potentally show graph
print(f'##End Takeoff Section')

print(f'\n##Begin climb to cruse altitude.')
heighttoclimb = atmo.Altitude - (MyTakeoff.height_Obs + MyTakeoff.Alt)
aveAlt = (atmo.Altitude + (MyTakeoff.height_Obs + MyTakeoff.Alt))/2
aveROC = MyStats.get_ROC_vel_alt(aveAlt,atmo.Vinfinity,"AVE") * ur.m / ur.s #Unrestricted climb at full power (Vel pinned to 25 m/s, all excess power going into climb)
TsecondClimb = heighttoclimb / aveROC
fullpowerclimbEnergy = (powerCons * TsecondClimb).to("J") #Full motor power for the duration of climb, Not taking into account effic as thats done in the "get_roc" func
Battery.discharge(fullpowerclimbEnergy)# Discharge that energy 
print(f" | Remaining altitude to gain: {heighttoclimb} with an average ROC of {aveROC}. Climb for {TsecondClimb}")
print(f' | Full power climb will consume {fullpowerclimbEnergy}')
Battery.print_state()
print(f'##End Climb to cruse altitude.') 

print(f"\n##Begin cruse")#Start cruise at 600 m alt
cruseTimeSec = 3600 * ur.sec
atmo.printAtmos()    
powerReqCruse = calc_PowerReq(atmo.Density,atmo.Vinfinity,MainWing.Area, OppaStoppa.Cd0,k,OppaStoppa.weight_takeoff).to("watt") #Get power for SLF
energyCruse = ((powerReqCruse / (Motor.effic * ROUGH_PROP_EFFIC)) * cruseTimeSec).to("J") #Joules for entire flight
print(f"Steady level cruse at {atmo.Altitude} altitude and {atmo.Vinfinity} requires {powerReqCruse} of power.")
print(f" -Power Draw from battery will be {powerReqCruse / (Motor.effic * ROUGH_PROP_EFFIC)} due to inefficency.")
print(f"Cruse for {cruseTimeSec}, using {energyCruse} of energy.")

Battery.discharge(energyCruse) #End 
Battery.print_state()
print(f"##End cruse")

print(f'\n##Begin Decent')
print(f'At this point we have very low energy remaining. Its ideal to glide to a landing.')
print(f'Optimal Glide is at CL/CD max -> Max range stats for our prop plane.')
print(f' | Cl/CD max: {clcdmax:.4f} at {(vinfMR.magnitude * ur.m / ur.s):0.2f}')
OptimGlideTheta = numpy.arctan(1/clcdmax)
RateOfDecent = vinfMR * numpy.sin(OptimGlideTheta)
print(f' | Resulting Glide ang: {numpy.rad2deg(OptimGlideTheta.magnitude):0.2f}°, Rate of Decent: {(RateOfDecent.magnitude * ur.m / ur.s):0.2f}.')

MyLanding = TakeOffLanding.Landing(OppaStoppa,300,1.2,numpy.rad2deg(OptimGlideTheta.magnitude),OppaStoppa.weight_takeoff)
AlttoDrop = (atmo.Altitude -MyLanding.obsH) # we can idealy drop right into a landing 
RangeGlide = (AlttoDrop / (numpy.tan(OptimGlideTheta)))
print(f' | Craft has {AlttoDrop} of elevation, it can glide for {RangeGlide}!')
print(f' | By glinding we preserve our power for any adjustments needed for landing.')
Battery.print_state()
print(f"##End Decent")

print(f"\n##Begin Landing")
print(f" |->If we apprach landing at our glide slope angle we can save energy, but its going to be a hot landing.")
MyLanding.apprachA = numpy.rad2deg(OptimGlideTheta.magnitude) #IF YOU WANT TO CHANGE APPROACH ANGLE DO IT HERE
print(f' |-->Landing Approach Angle: {MyLanding.apprachA}')
PowerReqTOland = calc_PowerReq(atmo.dens_trop_alt(MyLanding.Alt),MyLanding.V_A,MainWing.Area,OppaStoppa.Cd0,k,OppaStoppa.weight_takeoff).to("watt") #power req for SLF at landing conditions
thrReqTOland = (PowerReqTOland / MyLanding.V_A).to('newton') #Thrust for slf
threxcessReqTOApprach = OppaStoppa.weight_takeoff * numpy.sin(numpy.deg2rad(-MyLanding.apprachA)) #Weight is accually helping us here
reqThrustFromprop = (threxcessReqTOApprach + thrReqTOland) #how much thust do we really need?
powerReqApproach = (reqThrustFromprop * MyLanding.V_A)
print(" | Total Landing distance " + str(MyLanding.DoLanding()))
if SHOWGRAPHS:
    MyLanding.Graph()
print(f' | We apprached at {MyLanding.apprachA} degress at {MyLanding.V_A}.')
ApprTime = MyLanding.Approach(MyLanding.obsH - MyLanding.Flare()[1]) / MyLanding.V_A
ApprEnergy = (powerReqApproach * ApprTime).to("J")
print(f' | This required {ApprEnergy} of energy')
Battery.discharge(ApprEnergy)
Battery.print_state()
print("##End Landing")
print("\n##END Flight")
if SHOWGRAPHS:
    PRvsPA = MyStats.graph_PowerAval_vs_PowerReq(0,20000,300,25,"AVE",GRAPH_EXCESS=True)
    PRvsPA.legend()
    ROC = MyStats.graph_MAX_ROC_PROP(0,20000,300,"AVE")
    ROC.legend()
    PrVsVel = MyStats.graph_powerReq_vs_Vinf(600,10,70)
    PrVsVel.legend()
    MyLoadFactor = LoadFactor(OppaStoppa)
    VnDiag = MyLoadFactor.genVnDiagram(10,65,200,atmo.Altitude,"TAKEOFF",OppaStoppa.Clmin,OppaStoppa.CLmax)
    plt.show()