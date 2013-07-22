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

def main(ndarray[numpy.float, ndim=1] dist,
         ndarray[numpy.float, ndim=1] val, cdef float power):
    
    if min(dist) == 0:
        Res = val[argmin(dis)]
        return Res
        
    if min(dist) < 0:
        print 'Negative distances cannot be processed'
        return -9999
    
    d = numpy.array(dist)
    v = numpy.array(val)
    p = power
    
    ds = numpy.power(1./d,p)
    dt = numpy.sum(ds)
    W = ds/dt
    Res = cdef float numpy.sum(W*P)
    return Res
    