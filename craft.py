import math
from Components import *
from Drag import Drag
from StFl import *

class craft():
    name = "" #Name for craft
    mass = 0.0
    weight_empty = 4.44822 * 9414.356 #9414.356 lbs
    weight_takeoff = 0
    Cd0 = 0
    Cl = 0
    dragcomponents = list[DragyComponent] #List of of components of drag

    def __init__(self,name) -> None:
        self.name = name

    def list_components(self):
        for i in self.dragcomponents:
            print(i.Name)

    def compute_components(self):
        for i in self.dragcomponents:
            i.compute()

    def print_stats_components(self):
        for i in self.dragcomponents:
            i.printStats()

if __name__ == "__main__":
    OppaStoppa = craft("OpptaStoppa")

    MainWing = Wing3d()

    '''
    print(calc_PowerReq(1.19,90,25.26,0.00698,0.093,4.44822 * 9414.356))
    q = calc_dynpressure(1.19,90)
    print(calc_ThrustReq(q,25.23,0.0069849887954259,0.093,(4.44822 * 9414.356)))
    '''