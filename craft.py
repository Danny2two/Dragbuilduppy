import math
import os
os.environ['PINT_ARRAY_PROTOCOL_FALLBACK'] = "0"
import pint
import Atmosphere
from Components import *
from Drag import Drag
from StFl import *

class Craft():
    ur = pint.UnitRegistry()
    ur.setup_matplotlib()
    Atmosphere: Atmosphere
    name = "" #Name for craft 
    mass = 0.0 * ur.m
    weight_empty = 0 * ur.newton #4350 * 9.81 #kg * g
    weight_takeoff = 0 * ur.newton #5085 * 9.81 #kg * g
    Cd0 = 0
    Cl = 0
    PosStruturalNlimit = 1
    NegStruturalNlimit = -1
    dragcomponents: DragyComponent = []#List of of components of drag
    powertrain = [] #list of engines
    mainwing: Wing3d

    def __init__(self,name) -> None:
        self.name = name

    def list_components(self):
        for i in self.dragcomponents:
            print(i.Name)

    def compute_components(self):
        self.Cd0 = 0
        for i in self.dragcomponents:
            i.compute()
            #print(f'computing CD0 for {i.Name}')
            #print(f'  CD0: {i.getCD0()}')
            self.Cd0 += i.getCD0()

    def print_stats_components(self):
        for i in self.dragcomponents:
            i.printStats()

    def get_max_thrust(self):
        thr = 0
        for i in self.powertrain:
            thr += i.MaxThrust
        return thr

if __name__ == "__main__":
    OppaStoppa = Craft("OpptaStoppa")
    OppaStoppa.Atmosphere = Atmosphere(300,286.21, 9.77774/4,1.19,76,OppaStoppa.ur)
    ur = OppaStoppa.ur
    atmo = OppaStoppa.Atmosphere
    OppaStoppa.weight_empty = 4246.00052 * 9.81 * ur.newton
    OppaStoppa.weight_takeoff = 5021.00052 * 9.81 * ur.newton

    MainWing = Wing3d("Oppa Main Wing","NACA 4312",34,25.26,10,1.76,-4,17,0.12,0.3,6.56,atmo)
    MainFuselage = Fuselage("Oppa Fuselage",7.51,8.11,5.24,1.00,1.0,MainWing.Area,atmo)
    HorizontalTail = Wing3d("HT","NACA 0012",26.57,4.58,2.44,0.915,0,15,0.12,0.3,0.12,atmo)
    VerticalTail = Wing3d("VT","NACA 0012",26.57,4.58/2,2.44/2,0.915,0,15,0.12,0.3,0,atmo)
    TailGear = FixedGear("Tail Gear",0.25,0.196129,MainWing.Area,1.2,atmo) 

    OppaStoppa.dragcomponents = [MainWing,MainFuselage,HorizontalTail,VerticalTail,TailGear]

    OppaStoppa.compute_components()
    OppaStoppa.print_stats_components()
    print("Total CD0: " + str(OppaStoppa.Cd0))
    OppaStoppa.Atmosphere.printAtmos()
    k = calc_K_value(MainWing.OswaldE,MainWing.AR)
    q = calc_dynpressure(OppaStoppa.Atmosphere.Density,OppaStoppa.Atmosphere.Vinfinity)
    print("Power Req: " + str(calc_PowerReq(OppaStoppa.Atmosphere.Density,OppaStoppa.Atmosphere.Vinfinity,MainWing.Area,OppaStoppa.Cd0,k,OppaStoppa.weight_takeoff).to("watt")))
    print("Thrust req: " + str(calc_ThrustReq(q,MainWing.Area,OppaStoppa.Cd0,k,OppaStoppa.weight_takeoff).to("newton")))
    print("V Stall: " + str(calc_Vstall(OppaStoppa.Atmosphere.Density,OppaStoppa.weight_takeoff,MainWing.Area,1.3858).to('meters/second')))
    estCL = calc_req_CL(OppaStoppa.weight_takeoff,q,MainWing.Area)
    print("Estimated CL for routine flight: " + str(estCL))

