# -*- coding: utf-8 -*-
"""
============================
Variogram libraries function
============================
Implemented by Juan Chacon @ UNESCO-IHE
Integrated Water Systems and Governance Department
Hydroinformatics Laboratory

This library contains some semivariogram functions. At this stage the use of \
Matern semivariogram function is still presenting instabilities, and will be \
adjusted in a posterior stage.

* Pre requisites
    you will need the following libraries, not coming alongside with the\ 
    Anaconda ditribution (recommended)

* Functions
    * SVExponential: Expnential semivariogram computation.
    * SVGaussian: Gaussian semivariogram computation.
    * SVPower: Power semivariogram computation
    * SVSpherical: Spherical semivariogram computation
    * SVCubic: Cubic semivariogram computation
    * SVPentaspherical: Pentaspherical semivariogram computation
    * SVSinehole: Sinehole semivariogram computation
    * SVMatern: Matérn semivariogram computation
    * OptFunMaster: Function to calculate adjustment between observed and \
    simulated semivariogram

* Use policy
    * You should include the respective citation to the authors
    * If you find this tool usefull, you will give the main author a beer next\
    time you see him/her :)
    
* References
    * http://people.ku.edu/~gbohling/cpe940/Variograms.pdf
    * http://bit.ly/17z0aDw
"""

import numpy
import scipy.special
from numpy import linalg

## Performance Metrics
def RMSE(x,y):
    '''
    Calculates Root Mean Squared Error between two data series. \n
    
    Parameters
    ----------
        **x and y -- Data series which are intercambiable \n
    
    Returns
    -------
        **RMSE -- Value of the root mean squared error between data series
    '''
    Erro = numpy.square(numpy.subtract(x,y))
    if Erro.any < 0:
        return 9999
    cdef float F = numpy.sqrt(1.*sum(Erro)/len(x))
    return F

## Semivariogram lists
    # h = Lag dostance
    # S = Partial sill
    # R = Range
    # N = Nugget Effect
    # a = Dimensionless exponent
    # v = Matern parameter

def SVExponential(h,x):
    '''
    Calculate the value of the theoretical semivariogram 
    
    Parameters
    ----------
        **x -- vector of model parameters
        **h -- distance (lag) for computation of semivariogram
    
    Returns
    -------
        **SV -- value of semivariogram at lag x
    '''    
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]   
    cdef float SV
    if h == 0:
        return S
    if h/R > 3:
        return N+S
    try:
        SV = S - (-N + S * (1-numpy.exp(-1.*h/R)))
    except OverflowError:
        SV = N+S
    return SV
    
def SVGaussian(h,x):
    '''
    Calculate the value of the theoretical semivariogram 
    
    Parameters
    ----------
        **x -- vector of model parameters
        **h -- distance (lag) for computation of semivariogram
    
    Returns
    -------
        **SV -- value of semivariogram at lag x
    '''
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float SV
    if h == 0:
        return S
    if h/R > 1:
        return N+S
    try:
        SV = S - (-N + S * (1-numpy.exp(-1.*numpy.square(h)/R)))
    except OverflowError:
        SV = N+S
    return SV

def SVPower(h,x):
    '''
    Calculate the value of the theoretical semivariogram 
    
    Parameters
    ----------
        **x -- vector of model parameters
        **h -- distance (lag) for computation of semivariogram
    
    Returns
    -------
        **SV -- value of semivariogram at lag x
    '''
    cdef float S = x[0]
    cdef float N = x[2]
    cdef float a = x[3]
    cdef float SV
    if h == 0:
        return S
    if a > 2:
        a = 2
    try:
        SV = S - (-N + S*numpy.power(h,a))
    except OverflowError:
        SV = N+S
    return SV
    
def SVSpherical(h,x):
    '''
    Calculate the value of the theoretical semivariogram 
    
    Parameters
    ----------
        **x -- vector of model parameters
        **h -- distance (lag) for computation of semivariogram
    
    Returns
    -------
        **SV -- value of semivariogram at lag x
    '''
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float SV, hr
    if h == 0:
        return S+N
    if h > R:
        return 0
    hr = 1.*h/R
    SV = S - (-N + (S * ((3./2)*hr - ((1./2)*(hr**3.)))))
#    if SV < 0:
#        SV = 0
        
    return SV
 
def SVCubic(h,x):
    '''
    Calculate the value of the theoretical semivariogram 
    
    Parameters
    ----------
        **x -- vector of model parameters
        **h -- distance (lag) for computation of semivariogram
    
    Returns
    -------
        **SV -- value of semivariogram at lag x
    '''
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float SV
    if h == 0:
        return S
    if h > R:
        return N+S
    SV = S - (-N + S * (7.0*numpy.power((1.0*h/R),2.0) - 
        (35.0/4.0)*numpy.power(1.0*h/R,3.0) + (7.0/2.0)*numpy.power(h/R,5.0) - 
        (3.0/4.0)*numpy.power(h/R,7.0)))
    return SV 

def SVPentaspherical(h,x):
    '''
    Calculate the value of the theoretical semivariogram 
    
    Parameters
    ----------
        **x -- vector of model parameters
        **h -- distance (lag) for computation of semivariogram
    
    Returns
    -------
        **SV -- value of semivariogram at lag x
    '''
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float SV
    if h == 0:
        return S
    if h > R:
        return N+S
    SV = S - (-N + S * ((15./8)*(1.*h/R)-(5./4)*numpy.power(1.*h/R,3)+
        (3./8)*numpy.power(1.*h/R,5)))
    return SV
    
def SVSinehole(h,x):
    '''
    Calculate the value of the theoretical semivariogram 
    
    Parameters
    ----------
        **x -- vector of model parameters
        **h -- distance (lag) for computation of semivariogram
    
    Returns
    -------
        **SV -- value of semivariogram at lag x
    '''
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float SV
    if h == 0:
        return S
    if h > R:
        h = R
    SV = S - (-N + S * (1.-numpy.sin(numpy.pi*1.*h/R)/(numpy.pi*1.*h/R)))
    return SV

def SVMatern(h,x):
    '''
    Calculate the value of the theoretical semivariogram 
    
    Parameters
    ----------
        **x -- vector of model parameters
        **h -- distance (lag) for computation of semivariogram
    
    Returns
    -------
        **SV -- value of semivariogram at lag x
    '''
    cdef float S = x[0]
    cdef float R = x[1]
    cdef float N = x[2]
    cdef float v = x[4]
    cdef float SV
    
    if h == 0:
        return S
        
    if v < 0:
        return 9999
    
    a = 1/(scipy.special.gamma(v)*(2**(v-1)))
    b = numpy.sqrt(2*v)*h/R
    c = scipy.special.kv(v,b)
    SV = (S-N)*a*(b**v)*c
    return SV
    
def optFunMaster(x,SVExp,j,VarFunArr):
    '''
    Calculate the value of the RMSE between theoretical and experimental \
    semivariograms
    
    Parameters
    ----------
        **x -- Vector of model parameters
        **SVExp -- Experimental semivariogram vector
        **j -- Variogram model in the array
        **VarFunArr -- Array of semivariogram models
    
    Returns
    -------
        **RMSE -- Root mean squared error between observed and theoretical \
        semivariogram
    '''
    temp = []
    temp2 = []
    cdef float F
    cdef int fail = 0
    for i in xrange (0,len(SVExp)):
        temp.append(VarFunArr[j](SVExp[i][0],x))
        temp2.append(SVExp[i][1])
#    print min(temp)
#    if min(temp) < 0:
#        F = 9999
#        fail = 1
#    else:
    F = RMSE(temp,temp2)
    return F, [], fail