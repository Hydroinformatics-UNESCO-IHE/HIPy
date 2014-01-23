# -*- coding: utf-8 -*-
"""
Created on Mon Jun 03 14:17:49 2013

@author: chaco3

HBV-96
This is the HBV-96 implementation by Juan Chacon at UNESCO-IHE, Delft, NL.
This is meant to be part of the wesenseit project,of citizens observatories
"""

def Prec(double T,double LTT,double UTT,double P,double RFCF,double SFCF, double TFAC):
    '''
    Precipitation routine
    ---------------------
    
    Below the lower treshold level, all precipitation is snowfall.
    Similarly, all the precipitation above upper temperature treshold is
    rainfall. In between, there is a linear mixture between raifall and
    snowfall.
    
    Parameters
    ----------
        **T -- Measured temperature [C]
        **LTT -- Lower temperature treshold [C]
        **UTT -- Upper temperature treshold [C]
        **P -- Precipitation [mm]
        **RFCF -- Rainfall corrector factor
        **SFCF -- Snowfall corrector factor
    Returns
    -------
        **RF - Rainfall [mm]
        **SF - Snowfall [mm]
    '''
    
    cdef double RF,SF
    if T <= LTT:
        RF = 0.0
        SF = P*SFCF

    elif T >= UTT: 
        RF = P*RFCF
        SF = 0.0
        
    else:
        RF = ((T-LTT)/(UTT-LTT)) * P * RFCF
        SF = (1.0-((T-LTT)/(UTT-LTT))) * P * SFCF
    
    return RF/TFAC, SF/TFAC

# Snow Routine
def Snow(double CFMAX,double TFAC,double T,double TTM,double CFR,double CWH,
         double RF,double SF,double WCOld,double SPOld):
    '''
    Snow routine
    ------------
        At first, comparison of temperature is made. if temperature is below
        treshold, melting is happening, otherwise, refreezing. if the water
        content in the snow pack is bigger than water holding capacity, excess
        infiltrates soil.
        
    Parameters
    ----------
        **CFMAX -- Day degree factor
        **TFAC -- Temperature correction factor
        **T -- Temperature [C]
        **TTM -- Temperature treshold for Melting [C]
        **CFR -- Refreezing factor
        **CWH -- Capacity for water holding in snow pack
        **RF -- Rainfall [mm]
        **SF -- Snowfall [mm]
        **WCOld -- Water content in previous state [mm]
        **SPOld -- Snow pack in previous state [mm]
        
    Returns
    -------
        **IN -- Infiltration [mm]
        **WCNew -- Water content in posterior state [mm]
        **SPNew -- Snowpack in posterior state [mm]   
    '''
             
    cdef double MELT, SPNew, WCInt, REFR, IN, WCNEW
    if T > TTM:
    
        if CFMAX*(T-TTM)<SPOld+SF:
            MELT = CFMAX*(T-TTM)
        else:
            MELT = SPOld+SF
    
        SPNew = SPOld+SF-MELT
        WCInt = WCOld+MELT+RF
    
    else: 
        if CFR*CFMAX*(TTM-T)<WCOld+RF: 
            REFR = CFR*CFMAX*(TTM-T) 
        else:
            REFR = WCOld+RF 
          
        SPNew = SPOld+SF+REFR 
        WCInt = WCOld-REFR+RF 

    if WCInt>CWH*SPNew: 
        IN = WCInt-CWH*SPNew 
        WCNew = CWH*SPNew 
    else:
        IN = 0.0 
        WCNew = WCInt 
    
    return IN, WCNew, SPNew

## Soil routine
def Soil(double FC,double BETA,double ETF,double T,double TM,double ECORR,
         double LP,double TFAC,double CFLUX,double IN,double EP,double SMOld,
         double UZOld):
    '''
    Soil routine
    ------------
        At first, comparison of temperature is made. if temperature is below
        treshold, melting is happening, otherwise, refreezing. if the water
        content in the snow pack is bigger than water holding capacity, excess
        infiltrates soil.
        
    Parameters
    ----------
        **FC -- Filed capacity
        **BETA -- Shape coefficient for effective precipitation separation
        **ETF -- Total potential evapotranspiration
        **T -- Temperature
        **TM -- Average long term temperature
        **ECORR -- Evapotranspiration corrector factor
        **LP -- Soil wilting point
        **TFAC -- Time conversion factor
        **CFLUX -- Capilar flux in the root zone
        **IN -- actual infiltration
        **EP -- actual evapotranspiration
        **SMOld -- Previous soil moisture value
        **UZOld -- Previous Upper zone value
        
    Returns
    -------
        **SMNew -- New value of soil moisture
        **UZInt1 -- New value of direct runoff into upper zone
    '''
    
    cdef double R, EPInt, EA, CF, SMNew, UZInt1
    EP = EP/TFAC
    R = ((SMOld/FC)** BETA) * IN
    EPInt = (1.0+ETF*(T-TM))*ECORR*EP
    
    if SMOld/(LP*FC) < 1.0: 
        EA = (SMOld/(LP*FC))*EPInt 
    else: 
        EA = EPInt 
    
    if CFLUX*(1.0-(SMOld/FC)) < UZOld:
        CF = CFLUX*(1.0-(SMOld/FC))
    else:
        CF = UZOld 
    
    SMNew = SMOld+(IN-R)+CF-EA
    if SMNew < 0.0:
        SMNew = 0.0
        
    UZInt1 = UZOld+R-CF 
    
    return SMNew,UZInt1

def Soil_2(double FC,double BETA,double ETF,double T,double TM,double ECORR,
         double LP,double TFAC,double CFLUX,double IN,double EP,double SMOld,
         double UZOld):
    '''
    Soil routine
    ------------
        At first, comparison of temperature is made. if temperature is below
        treshold, melting is happening, otherwise, refreezing. if the water
        content in the snow pack is bigger than water holding capacity, excess
        infiltrates soil.
        
    Parameters
    ----------
        **FC -- Filed capacity
        **BETA -- Shape coefficient for effective precipitation separation
        **ETF -- Total potential evapotranspiration
        **T -- Temperature
        **TM -- Average long term temperature
        **ECORR -- Evapotranspiration corrector factor
        **LP -- Soil wilting point
        **TFAC -- Time conversion factor
        **CFLUX -- Capilar flux in the root zone
        **IN -- actual infiltration
        **EP -- actual evapotranspiration
        **SMOld -- Previous soil moisture value
        **UZOld -- Previous Upper zone value
        
    Returns
    -------
        **SMNew -- New value of soil moisture
        **UZInt1 -- New value of direct runoff into upper zone
    '''
    
    cdef double R, EPInt, EA, CF, SMNew, UZInt1
    EP = EP/TFAC
    R = ((SMOld/FC)** BETA) * IN
    EPInt = ECORR*EP
    
    if SMOld/(LP*FC) < 1.0: 
        EA = (SMOld/(LP*FC))*EPInt 
    else: 
        EA = EPInt 
    
    if CFLUX*(1.0-(SMOld/FC)) < UZOld:
        CF = CFLUX*(1.0-(SMOld/FC))
    else:
        CF = UZOld 
    
    SMNew = SMOld+(IN-R)+CF-EA
    if SMNew < 0.0:
        SMNew = 0.0
        
    UZInt1 = UZOld+R-CF 
    
    return SMNew,UZInt1

## Response routine
def Resp(double TFAC,double PERC,double ALPHA,double K,double K1,double AREA,
         double LZOld,double UZInt1):
             
    cdef double LZInt1, UZInt2, Q0, Q1, UZNew, LZNew, QNew
    if PERC < UZInt1: 
        LZInt1 = LZOld+PERC
    else:
        LZInt1 = LZOld+UZInt1 


    if UZInt1 > PERC: 
        UZInt2 = UZInt1-PERC
    else:
        UZInt2 = 0.0

    Q0 = K*(UZInt2**(1.0+ALPHA))
    Q1 = K1*LZInt1
    
    UZNew = UZInt2-(Q0) 
    LZNew = LZInt1-(Q1)
    
    QNew = AREA*(Q0+Q1)/(86.4)
    return QNew,UZNew,LZNew

"""
## Routing routine

#  Data load
#  v = Input vector
#  p = Parameter vector
#  St = State vector
#  x = Aditional parameters (values that multiply old states in the model)

# Input variables 
# P = Total precipitation v(0)
# T = actual temperature v(1)
# ETF = Total potential evapotranspiration v(2)   input
# TM = daily long term mean temperature v(3)  input

# Parameter Set
# TT = Limit temperature for rain/snow precipitation p(0)
# TTI = temperature treshold for linear mix of snow/rain precipitation p(1)
# TTM = Limit temperature for melting p(2)
# CFMAX = Degree day factor (measures the temperature variation along the day) p(3)
# FC = Field Capacity p(4)
# ECORR = Evapotranspiration corrector factor p(5)
# EP = Long term mean potential evapotranspiration p(6)
# LP = Soil moisture value where soil moisture reaches maximum potential
# evapotranspiration p(7)
# K = Upper zone response coefficient p(8)
# K1 = Lowe zone response coefficient p(9)
# ALPHA = upper zone runoff coefficient p(10)
# BETA = Controls the contribution of the increase in the soil moisture or
# to the response function p(11)
# CWH = Maximum amount of water that can be stored in snow pack p(12)
# CFR = Refreezing factor p(13)
# CFLUX = Capilar flux p(14)
# PERC = Percolation p(15)
# RFCF = Rainfal correction factor p(16)
# SFCF = Snowfall correction factor p(17)

# Non optimised parameters
# TFAC = Time factor p(18) = dt/86400
# AREA = Catchment area p(19)
"""

def Run(p,p2,v,St):
    '''
    This is the main module for the HBV run, as described in Lindstrom, 1997
    
    This script receives
        p = Parameter vector
        p2 = non optimisable parameter vector
        v = inputs
        St = Old states of the model
    
    This script returns
        QNew = Outflow
        St = Posterior states of the model
    '''
    
#    SPNew,SMNew,UZNew,LZNew,WCNew = 0.0,0.0,0.0,0.0,0.0   
#    QNew = 0.0
    
    ## Parse of parameters from input vector to model
    LTT = p[0]
    UTT = p[1]
    TTM = p[2]
    CFMAX = p[3]
    FC = p[4]
    ECORR = p[5]
    ETF = p[6]
    LP = p[7]
    K = p[8]
    K1 = p[9]
    ALPHA = p[10]
    BETA = p[11]
    CWH = p[12]
    CFR = p[13]
    CFLUX = p[14]
    PERC = p[15]
    RFCF = p[16]
    SFCF = p[17]
    
    ## Non optimisable parameters
    TFAC = p2[0]
    AREA = p2[1]
    
    ## Parse of Inputs
    Pr = v[0] # Precipitation [mm]
    T = v[1] # Temperature [C]
    EP = v[2] # Long terms (monthly) Evapotranspiration [mm]
    TM = v[3] #Long term (monthly) average temperature [C]

    ## Parse of states
    SPOld = St[0]
    SMOld = St[1]
    UZOld = St[2]
    LZOld = St[3]
    WCOld = St[4]
    
    RF,SF = Prec(T,LTT,UTT,Pr,RFCF,SFCF,TFAC)
    IN,WCNew,SPNew = Snow(CFMAX,TFAC,T,TTM,CFR,CWH,RF,SF,WCOld,SPOld)
    SMNew,UZInt1 = Soil(FC,BETA,ETF,T,TM,ECORR,LP,TFAC,CFLUX,IN,EP,SMOld,UZOld)
    QNew,UZNew,LZNew = Resp(TFAC,PERC,ALPHA,K,K1,AREA,LZOld,UZInt1)
    
    return QNew,[SPNew,SMNew,UZNew,LZNew,WCNew]

def Run_2(p,p2,v,St):
    '''
    This is the main module for the HBV run, as described in Lindstrom, 1997
    
    This script receives
        p = Parameter vector
        p2 = non optimisable parameter vector
        v = inputs
        St = Old states of the model
    
    This script returns
        QNew = Outflow
        St = Posterior states of the model
    '''
    
#    SPNew,SMNew,UZNew,LZNew,WCNew = 0.0,0.0,0.0,0.0,0.0   
#    QNew = 0.0
    
    ## Parse of parameters from input vector to model
    LTT = p[0]
    UTT = p[1]
    TTM = p[2]
    CFMAX = p[3]
    FC = p[4]
    ECORR = p[5]
    ETF = p[6]
    LP = p[7]
    K = p[8]
    K1 = p[9]
    ALPHA = p[10]
    BETA = p[11]
    CWH = p[12]
    CFR = p[13]
    CFLUX = p[14]
    PERC = p[15]
    RFCF = p[16]
    SFCF = p[17]
    
    ## Non optimisable parameters
    TFAC = p2[0]
    AREA = p2[1]
    
    ## Parse of Inputs
    Pr = v[0] # Precipitation [mm]
    T = v[1] # Temperature [C]
    EP = v[2] # Long terms (monthly) Evapotranspiration [mm]
    TM = v[3] #Long term (monthly) average temperature [C]

    ## Parse of states
    SPOld = St[0]
    SMOld = St[1]
    UZOld = St[2]
    LZOld = St[3]
    WCOld = St[4]
    
    RF,SF = Prec(T,LTT,UTT,Pr,RFCF,SFCF,TFAC)
    IN,WCNew,SPNew = Snow(CFMAX,TFAC,T,TTM,CFR,CWH,RF,SF,WCOld,SPOld)
    SMNew,UZInt1 = Soil_2(FC,BETA,ETF,T,TM,ECORR,LP,TFAC,CFLUX,IN,EP,SMOld,UZOld)
    QNew,UZNew,LZNew = Resp(TFAC,PERC,ALPHA,K,K1,AREA,LZOld,UZInt1)
    
    return QNew,[SPNew,SMNew,UZNew,LZNew,WCNew]


def HBVCalibPre(P, Flow, AvgPrec, Temper, Et,LLTemp, P2, St, WU):
    St = [[50.0,50.0,50.0,50.0,50.0]]
    cdef int fail = 0    
    QCal = []
    St2 = St[0][:]
    g = []
    for i in xrange(len(Flow)):
        v = [AvgPrec[i],Temper[i],Et[i],LLTemp[i]]
        HBVOut = Run(P,P2,v,St2)
        if HBVOut[0] < 0:
            fail = 1
            return 10000, g, fail
        if min(HBVOut[1])<0:    
            fail = 1
            return 10000, g, fail
        QCal.append(HBVOut[0])
        St2 = HBVOut[1]
    F = NSE(Flow[WU:len(Flow)],QCal[WU:len(Flow)])

    return -F, g, fail

import numpy
def NSE(x,y):
    '''
    x = measured
    y - Simulated
    '''
#    import numpy
    Erro = numpy.square(numpy.subtract(x,y))
    Erro2 = numpy.square(numpy.subtract(x,numpy.average(x)))
    if Erro.any < 0:
        return(-99999999)
    cdef double F = 1-(1.*numpy.sum(Erro)/numpy.sum(1.*Erro2))
    return(F)

def RMSE(x,y):
    '''
    x = measured
    y - Simulated
    '''
    Erro = numpy.array([],dtype='int64')
    Erro = numpy.square(numpy.subtract(x,y))
    if Erro.any < 0:
        return(99999999)
    F = numpy.sqrt(1.*sum(Erro)/len(x))
    return(F)