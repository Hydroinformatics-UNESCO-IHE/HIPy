# -*- coding: utf-8 -*-
"""
Performance Functions

Reference to the methods:

P Krause, DP Boyle, F Bäse. Comparison of different efficiency criteria for 
hydrological model assessment. Advances in Geosciences.5, 89-97, EGU 2005.

x - calculated value
y - recorded value
q - Quality tag (0-1)
"""
import numpy
import scipy.stats
# constants import 

# Or something similar,.... like importing a module of contants 
#or something
ERROR_CODE = -9999

def DataValid(x,y,q):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    
    if len(x) != len(y):
        print 'Calculated and recorded series length do not match'
        return ERROR_CODE
        
    if len(q) != len(y):
        print 'Quality tags vector length do not match with measurements'
        return ERROR_CODE
    
    if max(q) == 0:
        print 'Totally unreliable data, check quality tags'
        return ERROR_CODE
    
    if numpy.amin(q) < 0:
        print 'Quality tags cannot be negative'
        return ERROR_CODE
    
    if numpy.amax(q) > 1:
        print 'Quality tags cannot be greater than 1'        
        return ERROR_CODE
        
    try:
        numpy.sum(x)
    except ValueError:
        print 'Calculated data might contain non-numerical values'
        return ERROR_CODE
        
    try:
        numpy.sum(y)
    except ValueError:
        print 'Recorded data might contain non-numerical values'
        return ERROR_CODE
    
    try:   
        numpy.sum(q)
    except ValueError:
        print 'Quality tags might contain non-numerical values'
        return ERROR_CODE
        
    x = numpy.array(x)
    y = numpy.array(y)    
    q = numpy.array(q)
    return x,y,q

def RMSE(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))
    
    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE    
        
    Erro = numpy.square(x-y)*q    
    cdef float F = numpy.sqrt(1.*numpy.sum(Erro)/(numpy.sum(q)))
    return F
    
def NRMSE(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))
    
    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE    
        
    Erro = numpy.square(x-y)*q    
    cdef float F = (numpy.sqrt(1.*numpy.sum(Erro)/(numpy.sum(q))))/(numpy.amax(y)-numpy.amin(y))
    return F

def RSR(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))
    
    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE    
        
    Erro = numpy.square(x-y)*q    
    Erro2 = numpy.square(y-numpy.average(y))*q
    cdef float F = numpy.sqrt(1.*numpy.sum(Erro)/(numpy.sum(q)))/numpy.sqrt(1.*numpy.sum(Erro2))
    return F
    
def RSD(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))
    
    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
        

    cdef float F = numpy.std(x)/numpy.std(y)
    return F

def Bias(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))
    
    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
    
    cdef float F = numpy.sum(numpy.subtract(x,y)*q)/numpy.sum(q)
    return F

def PBias(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))
    
    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
        
    cdef float F = numpy.sum(numpy.subtract(x,y)*q)/numpy.sum(y*q)
    return F

def MAE(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))    
    
    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
    cdef float F = (1./sum(q))*numpy.abs(numpy.sum(x-y))
    return F
    
def MSE(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))    
    
    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE    
    Erro = numpy.square(x-y)*q    
    cdef float F = numpy.sum(Erro)/numpy.sum(q)
    return F   
    
def PercVol(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))   
    
    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
    cdef float sy = numpy.sum(y)
    cdef float F = (numpy.sum(x)-sy)/sy
    return F

def NSE(x,y,q='def',j=2.0):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    j - exponent to modify the inflation of the variance (standard NSE j=2)
    """
    if q is 'def':
        q = numpy.ones(len(y))   

    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
        
    a = numpy.sum(numpy.power(x-y,j)*q)
    b = numpy.sum(numpy.power(y-numpy.average(y),j)*q)
    cdef float F = 1.0 - a/b
    return F

def LNSE(x,y,q='def',j=2.0):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    j - exponent to modify the inflation of the variance (standard NSE j=2)
    """
    if q is 'def':
        q = numpy.ones(len(y))   

    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
    
    x = numpy.log(x)
    y = numpy.log(y)
    
    a = numpy.sum(numpy.power(x-y,j)*q)
    b = numpy.sum(numpy.power(y-numpy.average(y),j)*q)
    cdef float F = 1.0 - a/b
    return F

def RNSE(x,y,q='def',j=2.0):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    j - exponent to modify the inflation of the variance (standard NSE j=2)
    """
    if q is 'def':
        q = numpy.ones(len(y))   

    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
        
    a = numpy.sum(numpy.power(numpy.subtract(x,y)/numpy.average(y),j)*q)
    b = numpy.sum(numpy.power(numpy.subtract(y,numpy.average(y))/numpy.average(y),j)*q)
    cdef float F = 1.0 - a/b
    return F

def IOA(x,y,q='def',j=2.0):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))   

    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
        
    a = numpy.sum(numpy.power(numpy.subtract(x,y)*q,j))
    b = numpy.sum(numpy.power(numpy.abs(numpy.subtract(x,y)*q)+numpy.abs(numpy.subtract(y,numpy.average(y))*q),j))
    cdef float F = 1.0 - a/b
    return F

def RIOA(x,y,q='def',j=2.0):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))   

    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
        
    a = numpy.sum(numpy.power(numpy.subtract(x,y)*q/y,j))
    b = numpy.sum(numpy.power(numpy.abs(numpy.subtract(x,y)*q/y)+numpy.abs(numpy.subtract(y,numpy.average(y))*q/y),j))
    cdef float F = 1.0 - a/b
    return F
    
def PearsonR(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))   

    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
        
    return scipy.stats.pearsonr(x,y)[0]

def SpearmanR(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))   

    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
        
    return scipy.stats.spearmanr(x,y)[0]

def DetCoef(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))   

    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
        
    return numpy.power(scipy.stats.pearsonr(x,y)[0],2)

def COP(x,y,q='def'):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    """
    if q is 'def':
        q = numpy.ones(len(y))   

    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
    
    a = numpy.sum(numpy.power(numpy.subtract(y[1:]-x[1:])*q,2))
    b = numpy.sum(numpy.power(numpy.subtract(y[1:]-y[:-1])*q,2))
    cdef float F = 1.0 - a/b
    return F

def KGE(x,y,q='def',s=[1,1,1]):
    """
    Performance Functions
    x - calculated value
    y - recorded value
    q - Quality tag (0-1)
    s - inflation of different parameters (see:Gupta, 2009)
    """
    if q is 'def':
        q = numpy.ones(len(y))   

    x,y,q = DataValid(x,y,q)
    if min(x) == ERROR_CODE:
        return ERROR_CODE
    
    if min(q) == 0:
        for i in xrange(0,len(x)):
            CleanDum = []
            if q[i] == 0:
                CleanDum.append(i)
        
        x = numpy.delete(x,CleanDum)
        y = numpy.delete(y,CleanDum)
            
    r = scipy.stats.pearsonr(x,y)[0]
    alp = numpy.std(x)/numpy.std(y)
    bet = (numpy.mean(x)-numpy.mean(y))/numpy.std(y)
    ED = numpy.sqrt((s[0]*(r-1))**2+(s[1]*(alp-1))**2+(s[2]*(bet-1))**2)
    KGE = 1-ED
    return KGE