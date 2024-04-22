import math
import numpy
import pint

def calc_req_CL(weight, Q_inf, wingarea):
    """Calculates requred CL for steady level flight

    Args:
        weight (_type_): weight of craft in newton
        Q_inf (_type_): Dynamic pressure
        wingarea (_type_): Main wing area

    Returns:
        _type_: Estimated CL
    """    
    CL = (weight/(Q_inf * wingarea))
    return CL

def calc_ThrustReq(Q_inf, wingarea,CD_0,K,weight):
    """Calculates requred thrust for SLF

    Args:
        Q_inf (_type_): Dynamic pressure
        wingarea (_type_): mainwing area
        CD_0 (_type_): CD_0 for craft
        K (_type_): K coef
        weight (_type_): Craft weight in newtons

    Returns:
        _type_: thrust requred in newtons
    """    
    Thrust = Q_inf * wingarea * CD_0 + ((K*numpy.power(weight,2))/(Q_inf*wingarea))
    return Thrust

def calc_Vstall(P_inf, weight,wingarea,CL_max):
    """Calculates the stall velocity of the craft

    Args:
        P_inf (_type_): Density of air
        weight (_type_): weight of craft in newton
        wingarea (_type_): wing area in square meters
        CL_max (_type_): Max CL

    Returns:
        _type_: Stall Velocity
    """    
    #print(f'pinf: {P_inf}, w: {weight}, Wingarea: {wingarea}')
    vstall = numpy.sqrt((2/P_inf) * (weight/wingarea) * (1/CL_max))
    return vstall

def calc_PowerReq(P_inf,V_inf,wingarea,CD_0,K,weight:float):
    """Calculates power requred for flight 

    Args:
        P_inf (_type_): Atmospheric density (kg/m^3)
        V_inf (_type_): Velocity @ infinity
        wingarea (_type_): Wing planform area
        CD_0 (_type_): Crafts Cd_0
        K (_type_): k value
        weight (_type_): craft weight (N)

    Returns:
        _type_: _description_
    """    
    #print(f'pinf: {P_inf}, vinf: {V_inf}, warea: {wingarea}, weight: {weight}')
    pr = 0.5*P_inf*numpy.power(V_inf,3)*wingarea*CD_0 + ((2*K*numpy.power(weight,2))/(P_inf*V_inf*wingarea))
    return pr

def calc_dynpressure(P_inf,V_inf):
    """Calculates dynamic pressure

    Args:
        P_inf (_type_): Density (kg/m^3)
        V_inf (_type_): Velocity m/s

    Returns:
        _type_: Dynamic pressure
    """    
    dynpress = 0.5*P_inf*numpy.power(V_inf,2)
    return dynpress

def calc_K_value(OswaldEff,AspectRatio):
    """Calculates K value

    Args:
        OswaldEff (_type_): Oswalds efficany
        AspectRatio (_type_): _description_

    Returns:
        _type_: _description_
    """    
    k= 1/(OswaldEff * math.pi * AspectRatio)
    return k

def calc_CL_CDmax(K, CD0):
    """Finds (CL/CD)max given K and CDO usefull for:
    PROP: max range
    JET: max endurace, Velocity of climb at ROC max

    Args:
        K (_type_): K value
        CD0 (_type_): C_D0

    Returns:
        Float: Ratio = (CL/CD)max
    """    
    clcd = numpy.sqrt(1/(4*K*CD0))
    return clcd

def calc_Vinf_MaxJetEndurance(P_inf,weight,wingArea,Cd_0,K):
    """Calculates V_Inf to achieve max endurance with a jet

    Args:
        P_inf (float): Atmospheric Density (kg/m^3)
        weight (float): weight of craft (N)
        wingArea (float): Area of wing (m^2)
        Cd_0 (_type_): _description_
        K (_type_): _description_

    Returns:
        float: V infinity to get max endurance (m/s)
    """    
    vinf = numpy.sqrt((2/P_inf) * (weight/wingArea)* numpy.sqrt(K/Cd_0))
    return vinf

def calc_Vinf_MaxJetRange(P_inf,weight,wingArea,Cd_0,K):
    """Calculates V_Inf to achieve max range with a jet

    Args:
        P_inf (float): Atmospheric Density (kg/m^3)
        weight (float): weight of craft (N)
        wingArea (float): Area of wing (m^2)
        Cd_0 (_type_): _description_
        K (_type_): _description_

    Returns:
        float: V infinity to get max range (m/s)
    """    
    vinf = numpy.sqrt((2/P_inf) * (weight/wingArea)* numpy.sqrt((3*K)/Cd_0))
    return vinf

def calc_CL_CDmaxRangeJet(K,CD_0):
    """Calculates Cl^1/2 / CD max usefull for max range of jets

    Args:
        K (_type_): K value
        CD_0 (_type_): CD_0 value

    Returns:
        _type_: (Cl^1/2 / CD)max
    """
    clcd = (3/4) * numpy.power((1/(3 * K * numpy.power(CD_0,3))),(1/4))
    return clcd


def calc_endurance_turbojet(TSFC,Cl_Cd,weightTO,weightEm):
    """Calculates Enduance for a turbojet with the provided characteristitcs

    Args:
        TSFC (g/N/s): Thrust Specific fuel Consumption *IN N fuel / sec / N thrust*
        Cl_Cd (Coeff): Cl / Cd max
        weightTO (N): Takeoff weight
        weightEm (N): Empty weight

    Returns:
        float: Endurance (hours)
    """    
    endur =(1/TSFC) * Cl_Cd * (math.log((weightTO/weightEm),math.e))
    return endur

def calc_max_range_jet(p_inf,wingarea,TSFC,clcd,weighttakeoff,weightempty):
    """Calculates range for a turbojet at a given parameters

    Args:
        p_inf (_type_): Atmospheric Density (kg/m^3)
        wingarea (_type_): _description_
        TSFC (_type_): Thrust specific fuel consumption *IN N fuel / sec / N thrust*
        clcd (_type_): Cl/Cd ratio (Cl^1/2 / Cd for max range)
        weighttakeoff (_type_): Craft weight at takeoff (N)
        weightempty (_type_): Craft empty weight (N)

    Returns:
        _type_: Range
    """    
    range = (2/1) * numpy.sqrt(2/(p_inf*wingarea)) * (1/TSFC) * clcd * (numpy.sqrt(weighttakeoff) - numpy.sqrt(weightempty))
    return range

def calc_endurance_electric_prop(Batt_Energy,Motor_eff,Prop_eff,p_inf,wingarea, weight,CLCD):
    """Calcualtes Endurance for a Electric Prop power train

    Args:
        Batt_Energy (_type_): Total energy in battery (joules)
        Motor_eff (_type_): Motor Effic
        Prop_eff (_type_): Prop Effic
        p_inf (_type_): Atmosphereic density
        wingarea (_type_): Mainwing area (m^2)
        weight (_type_): Craft Weight (N)
        CLCD (_type_): CL^(3/2) / CD max

    Returns:
        _type_: Endurace in seconds 
    """    
    Endurance = ((Batt_Energy * Motor_eff * Prop_eff * numpy.sqrt(p_inf * wingarea))/(numpy.sqrt(2)* numpy.power(weight,3/2))) * CLCD
    return Endurance

def calc_range_electric_prop(Batt_Energy,Motor_eff,Prop_eff,weight,CLCD):
    """Calcualtes Range of an electric prop power trian

    Args:
        Batt_Energy (_type_): Total Battery energy (joules)
        Motor_eff (_type_): Motor eff
        Prop_eff (_type_): Prop eff
        weight (_type_): Craft weight (N)
        CLCD (_type_): CL/CD max

    Returns:
        _type_: _description_
    """    
    Range = ((Batt_Energy * Motor_eff * Prop_eff)/weight) * CLCD
    return Range
    

def calc_Lift(CL,V_inf,Wingarea,dens):
    """Calculates Lift at given conditions

    Args:
        CL (_type_): Coefficent of lift
        V_inf (_type_): Incoming airstream velocity
        Wingarea (_type_): Manwing area
        dens (_type_): Air density 

    Returns:
        _type_: Lift in newtons
    """    
    Lift = 0.5 * dens * numpy.power(V_inf,2) * Wingarea * CL
    return Lift

def calc_Drag(CD0,V_inf,Wingarea,dens,K,CL):
    """Calculates Drag for given conditions

    Args:
        CD0 (_type_): Crafts CD_0
        V_inf (_type_): Incoming air velocity
        Wingarea (_type_): Manwing area
        dens (_type_): Density of air
        K (_type_): K coeff
        CL (_type_): CL at curent conditions

    Returns:
        _type_: _description_
    """    
    Drag = 0.5 * dens * numpy.power(V_inf,2) * Wingarea * (CD0 + (K * numpy.power(CL,2)))
    return Drag