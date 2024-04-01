import math
"""ATMOSPHERE CONDITIONS

Alt = 300 #meters
Temp = 286.21 #Kelvin
Press = 9.77774/4 #N/m^2
Dens = 1.19 #kg/m^3  𝜌
Default_Vinf = 76 #m/s
"""
class Atmosphere:
    Altitude = 0
    Temp = 0
    Pressure = 0
    Density = 0 #𝜌
    Vinfinity = 0

    #Refrence Conditions
    ThrustA_SeaLevel = 0 #N
    Temp_SeaLevel = 294 #K
    Density_SeaLevel = 1.225 #kg/m^3
    g = 8.91 #m/s/s
    a = -6.5*(math.pow(10,-3)) #K/m
    R = 287 #J/kg K

    
    def __init__(self,Altitude, Temp, Pressure, Density, Vinf) -> None:
        """Variables for the atmosphere

        Args:
            Altitude (_type_): Altitude (meters)
            Temp (_type_): Temperature (kelvin)
            Pressure (_type_): Pressure (N/m^2)
            Density (_type_): Desnity of atmosphere (𝜌, kg/m^3)
            Vinf (_type_): Velocity of incoming airstream (m/s)
        """        
        self.Altitude = Altitude
        self.Temp = Temp
        self.Pressure = Pressure
        self.Density = Density
        self.Vinfinity = Vinf
    
    def printAtmos(self):
        print("Atmosphere State:\n ALT: " + str(self.Altitude) +
               " m\n Temp: " + str(self.Temp) + 
               " k\n Pressure: " + str(self.Pressure)+ 
               " N/m^2\n Density: " + str(self.Density)+ 
               " kg/m^3\n V_inf: " + str(self.Vinfinity) + "m/s")        

    def temp_trop(self, altitude):
        """Calculates temperature from altitude in troposphere
        Args:
        altitude: altitude in meters
        Returns:
        temperature: temp in Kelvin
        """
        T = self.Temp_SeaLevel + (self.a*(altitude))
        return T

    def dens_trop(self, temperature):
        """Calculates density from temperature in troposphere
        
        Args:
        temperature: temp in kelvin
        Returns:
        density: density kg/m^3
        """
        dens = self.Density_SeaLevel * math.pow((temperature/self.Temp_SeaLevel), (((-1 * self.g)/(self.a *-1))))
        return dens
    
    def dens_trop_alt(self, alt):
        """Calculates density from temperature in troposphere
        Args:
        alt: Altitude in meters above sea level
        Returns:
        density: density kg/m^3
        """
        temperature = self.temp_trop(alt)
        dens = self.Density_SeaLevel * math.pow((temperature/self.Temp_SeaLevel), (((-1 * self.g)/(self.a *-1))))
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
        thrust = refThrust * math.pow((density/self.Density_SeaLevel),m)
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
        thrust = refThrust * math.pow((density/self.Density_SeaLevel),m)
        return thrust