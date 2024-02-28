import math

def calc_CL(weight, Q_inf, wingarea):
    CL = (weight/(Q_inf * wingarea))
    return CL

def calc_ThrustReq(Q_inf, wingarea,CD_0,K,weight):
    Thrust = Q_inf * wingarea * CD_0 + ((K*math.pow(weight,2))/(Q_inf*wingarea))
    return Thrust

def calc_Vstall(P_inf, weight,wingarea,CL_max):
    vstall = math.sqrtI((2/P_inf) * (weight/wingarea) * (1/CL_max))
    return vstall

def calc_PowerReq(P_inf,V_inf,wingarea,CD_0,K,weight):
    pr = 0.5*P_inf*math.pow(V_inf,3)*wingarea*CD_0 + ((2*K*math.pow(weight,2))/(P_inf*V_inf*wingarea))
    return pr

def calc_dynpressure(P_inf,V_inf):
    dynpress = 0.5*P_inf*math.pow(V_inf,2)
    return dynpress