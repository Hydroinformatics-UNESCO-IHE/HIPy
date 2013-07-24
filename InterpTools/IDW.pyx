# -*- coding: utf-8 -*-
"""
IDW
"""
'''
dist - vector with distances to interpolation target
val - Measured values at given data points
power - power of interpolation weighting
'''
import numpy

def Interp(dist, val, power):
    
    if numpy.min(dist) == 0:
        Res = val[numpy.argmin(dist)]
        return Res
        
    elif numpy.min(dist) < 0:
        print 'Negative distances cannot be processed'
        return -9999
    
    else:
        d = numpy.array(dist)
        v = numpy.array(val)
        p = power
        
        ds = numpy.power(1./d,p)
        dt = numpy.sum(ds)
        W = ds/dt
        Res = numpy.sum(W*v)
        return Res