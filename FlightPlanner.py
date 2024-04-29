import numpy as np
import matplotlib.pyplot as plt
import pint 
from craft import Craft
import TakeOffLanding
import CraftStatistics
import Atmosphere
import StFl

#Kinda assembling the Avengers with these import statements

class FlightPlanner:
    def __init__(self,Craft: Craft) -> None:
        self.FlightCraft = Craft
        self.ur = Craft.ur
        self.atmosphere: Atmosphere = Craft.Atmosphere
        self.Takeoff = TakeOffLanding.Takeoff
        self.Landing = TakeOffLanding.Landing
        pass