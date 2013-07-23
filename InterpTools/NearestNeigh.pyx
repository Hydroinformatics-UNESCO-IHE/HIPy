# -*- coding: utf-8 -*-
"""
Nearest Neighborhood
"""

from numpy import min, argmin

def main(ndarray[numpy.float, ndim=1] dist,
         ndarray[numpy.float, ndim=1] val):
    
    if min(dist) < 0:
        print 'Negative distances cannot be processed'
        return -9999
        
    Res = val[argmin(dis)]
    return Res
        

