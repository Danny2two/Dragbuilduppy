import math
from Components import *
from Drag import Drag
from StFl import *

class craft():
    name = "" #Name for craft
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