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
import tDist

def Interp(x, y, val, xt, yt, power=2.0, cutoff=0.00001):
    
    vect = numpy.transpose(numpy.vstack((x,y)))
    vecttar = numpy.transpose(numpy.vstack((xt,yt)))
    dis = tDist.target(vect, vecttar)
    
    R = []
    for i in xrange(0,len(xt)):

        if numpy.min(dis[i]) < cutoff:
            R.append(val[numpy.argmin(dis)])
            
        else:
            d = numpy.array(dis)
            v = numpy.array(val)
            
            ds = numpy.power(1./d,power)
            dt = numpy.sum(ds)
            W = ds/dt
            R.append(numpy.sum(W*v))
    
    return R

#def test():
#    x = [1,2,3]
#    y = [1,2,3]
#    val = [5,3,1]
#    xt = [1,4,5]
#    yt = [3,4,6]
#    
#    k = Interp(x,y,val,xt,yt)
#    return k