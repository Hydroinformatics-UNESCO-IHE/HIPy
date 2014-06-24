# -*- coding: utf-8 -*-
"""
Created on Wed May 14 08:54:12 2014

@author: chaco3

Noise Generator

Contains White Noise
Pink Noise
and Brown Noise
"""

import numpy
import matplotlib.pyplot as plt
#White Noise
def White(x,a):
    '''
    A data for which the noise wants to be added has to be introduced
    x - vector for which noise has to be added
    a - Std Deviation of generated noise
    '''
    y = x + a * numpy.random.randn(len(x))
    return y

def Pink(x,a):
    '''
    A data for which the noise wants to be added has to be introduced
    
    composed of sum of white noises for each step    
    
    x - vector for which noise has to be added
    a - Std Deviation of max generated noise
    '''
    y = numpy.zeros(len(x))
    for i in xrange(1,len(x)):
        WN = a*numpy.random.randn(i)
#        print WN
        y += numpy.interp(range(len(x)),numpy.linspace(0,len(x),i),WN)
    return y/a

def Brown(x,a):
    '''
    A data for which the noise wants to be added has to be introduced
    x - vector for which noise has to be added
    a - maximum amplitude of the noise
    '''
    return

def Test():
    Data = numpy.zeros(1000)
#    WNoisy = White(Data,10)
#    plt.plot(WNoisy)
#    plt.show()
    
    PNoisy = Pink(Data,10)
    plt.plot(PNoisy)
    plt.show()
    
    return

