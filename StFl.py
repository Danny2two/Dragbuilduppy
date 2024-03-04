import math

def calc_requred_CL(weight, Q_inf, wingarea):
    CL = (weight/(Q_inf * wingarea))
    return CL

def calc_ThrustReq(Q_inf, wingarea,CD_0,K,weight):
    Thrust = Q_inf * wingarea * CD_0 + ((K*math.pow(weight,2))/(Q_inf*wingarea))
    return Thrust

def calc_Vstall(P_inf, weight,wingarea,CL_max):
    vstall = math.sqrt((2/P_inf) * (weight/wingarea) * (1/CL_max))
    return vstall

def calc_PowerReq(P_inf,V_inf,wingarea,CD_0,K,weight):
    pr = 0.5*P_inf*math.pow(V_inf,3)*wingarea*CD_0 + ((2*K*math.pow(weight,2))/(P_inf*V_inf*wingarea))
    return pr

def calc_dynpressure(P_inf,V_inf):
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
    """Finds (CL/CD)max given K and CDO

    Args:
        K (_type_): K value
        CD0 (_type_): C_D0

    Returns:
        Float: Ratio = (CL/CD)max
    """    
    clcd = math.sqrt(1/(4*K*CD0))
    return clcd