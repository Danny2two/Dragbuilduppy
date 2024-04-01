import math

def calc_req_CL(weight, Q_inf, wingarea):
    CL = (weight/(Q_inf * wingarea))
    return CL

def calc_ThrustReq(Q_inf, wingarea,CD_0,K,weight):
    Thrust = Q_inf * wingarea * CD_0 + ((K*math.pow(weight,2))/(Q_inf*wingarea))
    return Thrust

def calc_Vstall(P_inf, weight,wingarea,CL_max):
    vstall = math.sqrt((2/P_inf) * (weight/wingarea) * (1/CL_max))
    return vstall

def calc_PowerReq(P_inf,V_inf,wingarea,CD_0,K,weight):
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
    pr = 0.5*P_inf*math.pow(V_inf,3)*wingarea*CD_0 + ((2*K*math.pow(weight,2))/(P_inf*V_inf*wingarea))
    return pr

def calc_dynpressure(P_inf,V_inf):
    """Calculates dynamic pressure

    Args:
        P_inf (_type_): Density (kg/m^3)
        V_inf (_type_): Velocity m/s

    Returns:
        _type_: Dynamic pressure
    """    
    dynpress = 0.5*P_inf*math.pow(V_inf,2)
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
    JET: max endurace 

    Args:
        K (_type_): K value
        CD0 (_type_): C_D0

    Returns:
        Float: Ratio = (CL/CD)max
    """    
    clcd = math.sqrt(1/(4*K*CD0))
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
    vinf = math.sqrt((2/P_inf) * (weight/wingArea)* math.sqrt(K/Cd_0))
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
    vinf = math.sqrt((2/P_inf) * (weight/wingArea)* math.sqrt((3*K)/Cd_0))
    return vinf

def calc_CL_CDmaxRangeJet(K,CD_0):
    """Calculates Cl^1/2 / CD max usefull for max range of jets

    Args:
        K (_type_): K value
        CD_0 (_type_): CD_0 value

    Returns:
        _type_: (Cl^1/2 / CD)max
    """
    clcd = (3/4) * math.pow((1/(3 * K * math.pow(CD_0,3))),(1/4))
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
    range = (2/1) * math.sqrt(2/(p_inf*wingarea)) * (1/TSFC) * clcd * (math.sqrt(weighttakeoff) - math.sqrt(weightempty))
    return range