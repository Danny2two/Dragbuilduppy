import math
import pint
import numpy as np
"""ATMOSPHERE CONDITIONS

Alt = 300 #meters
Temp = 286.21 #Kelvin
Press = 9.77774/4 #N/m^2
Dens = 1.19 #kg/m^3  𝜌
Default_Vinf = 76 #m/s
"""
class Atmosphere:
    Altitude = 0 #*ur.meter
    Temp = 0 #* ur.kelvin
    Pressure = 0 #* ur.newton / (ur.meter ** 2)
    Density = 0 #* ur.kg / (ur.meter ** 3)#𝜌
    Vinfinity = 0 #* ur.m / ur.s
    '''
    #Refrence Conditions
    ThrustA_SeaLevel = 0 #N
    Temp_SeaLevel = 294 #K
    Density_SeaLevel = 1.225 #kg/m^3
    g = 8.91 #* ur.m / ur.s /ur.s #m/s/s
    a = -6.5*(math.pow(10,-3)) #* ur.kelvin / ur.m #K/m
    R = 287 #J/kg K'''

    
    def __init__(self,Altitude, Temp, Pressure, Density, Vinf,UNITREG:pint.UnitRegistry) -> None:
        """Variables for the atmosphere
            GIVE UNITLESS QUANITIES, UNITS ASIGNED BY ATMOSPHERE CLASS

        Args:
            Altitude (_type_): Altitude (meters)
            Temp (_type_): Temperature (kelvin)
            Pressure (_type_): Pressure (N/m^2)
            Density (_type_): Desnity of atmosphere (𝜌, kg/m^3)
            Vinf (_type_): Velocity of incoming airstream (m/s)
        """        
        self.ur = UNITREG
        ur = self.ur

        self.Temp_SeaLevel = 294 * ur.kelvin
        self.Density_SeaLevel = 1.225 * ur.kg / ur.meter ** 3
        self.g = 8.91 * ur.m / ur.s /ur.s #m/s/s
        self.a = -6.5*(math.pow(10,-3)) * ur.kelvin / ur.m #K/m
        self.R = 287 * ur.joules / (ur.kg * ur.kelvin)#J/kg 

        self.Altitude = Altitude * ur.meters
        self.Temp = Temp * ur.kelvin
        self.Pressure = Pressure * ur.newton / (ur.meter ** 2)
        self.Density = Density * ur.kg / (ur.meter ** 3)
        self.Vinfinity = Vinf * ur.m / ur.s
    
    def printAtmos(self):
        print(f"Atmosphere State:\n ALT: {self.Altitude:P}")
        print(f"Temp: {self.Temp:P}") 
        print(f"Pressure: {self.Pressure:P}") 
        print(f"Density: {self.Density:P}")
        print(f"V_inf: {self.Vinfinity:P}")        

    def temp_trop(self, altitude):
        """Calculates temperature from altitude in troposphere
        Args:
        altitude: altitude in meters
        Returns:
        temperature: temp in Kelvin
        """

        if(self.ur.get_dimensionality(altitude) != self.ur.get_dimensionality( 1 * self.ur.meters)):
            altitude = altitude * self.ur.meters
        T = self.Temp_SeaLevel + (self.a*(altitude))
        return T

    def dens_trop(self, temperature):
        """Calculates density from temperature in troposphere
        
        Args:
        temperature: temp in kelvin
        Returns:
        density: density kg/m^3
        """
        dens = self.Density_SeaLevel * np.power((temperature/self.Temp_SeaLevel), (((-1 * self.g)/(self.a * self.R))-1))
        return dens
    
    def dens_trop_alt(self, alt):
        """Calculates density from temperature in troposphere
        Args:
        alt: Altitude in meters above sea level
        Returns:
        density: density kg/m^3
        """
        temperature = self.temp_trop(alt)
        dens = self.dens_trop(temperature)
        return dens

    def jet_thrust_available(self,refThrust, density, m):
        """Calculates thrust available
        Args:
        refThrust: refrence thrust at sea level
        density: density at alt
        m:
        Returns:
        Thrust: Thrust available
        """
        thrust = refThrust * np.power((density/self.Density_SeaLevel),m)
        return thrust
    
    def jet_thrust_available_alt(self,refThrust, alt, m):
        """Calculates thrust available
        Args:
        refThrust: refrence thrust at sea level
        alt: Altitude in meters above sea level
        m:
        Returns:
        Thrust: Thrust available
        """
        density = self.dens_trop_alt(alt)
        thrust = refThrust * np.power((density/self.Density_SeaLevel),m)
        return thrust