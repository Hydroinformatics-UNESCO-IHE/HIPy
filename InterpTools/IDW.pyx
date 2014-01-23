# -*- coding: utf-8 -*-
'''
================================
Inverse distance weighting (IDW)
================================

Implemented by Juan Chacon @ UNESCO-IHE
Integrated Water Systems and Governance Department
Hydroinformatics Laboratory

This library contains tools for interpolation using inverse distance weighting.

* Pre requisites
    you will need the following libraries, not coming alongside with the \ 
    Anaconda ditribution (recommended)

* Functions
    *Interp
    *Interp_bat

* Use policy
    * You should include the respective citation to the authors
    * If you find this tool usefull, you will give the main author a beer next\
    time you see him/her :)
    
* References
    * http://people.ku.edu/~gbohling/cpe940/Kriging.pdf

dist - vector with distances to interpolation target
val - Measured values at given data points
power - power of interpolation weighting
'''
import numpy
import Dist
import csv
import DataSave
import DataLoad

def Interp(x, y, val, xt, yt, power=2.0, cutoff=0.00001):
    '''
    Algorithm for IDW interpolation. This contains the main build to perform
    interpolation using any type of IDW procedure.
    
    Parameters
    ----------
        **x --
        **y --
        **val --
        **xt --
        **yt --
        **power --
        **cutoff --
    '''
    vect = numpy.transpose(numpy.vstack((x,y)))
    vecttar = numpy.transpose(numpy.vstack((xt,yt)))
    dis = Dist.target(vect, vecttar)

    R = []
    for i in xrange(0,len(xt)):

        if numpy.min(dis[i]) < cutoff:
            R.append(val[numpy.argmin(dis)])
            
        else:
            d = numpy.array(dis[i])
            v = numpy.array(val)
            
            ds = numpy.power(1./d,power)
            dt = numpy.sum(ds)
            W = ds/dt
            R.append(numpy.sum(W*v))
    
    return R

def Interp_bat(Loc, POI, Prec, power=2.0, cutoff=0.00001, tmin=0, tmax='def',
               valmin=-1e99, valmax=1e99):
    '''
    Parser for interpolation of several scenarios.
    '''
    Z = []
    ZAvg = []    
    if tmax == 'def':
        tmax = len(Prec)
        
    Loc = numpy.array(Loc)
    POI = numpy.array(POI)    
    
    x = Loc[:,0]
    y = Loc[:,1]
    
    xt = POI[:,0]
    yt = POI[:,1]
    
    for i in xrange(tmin,tmax):
        temp = Interp(x,y,Prec[i],xt,yt,power,cutoff)
        numpy.clip(temp,valmin,valmax)        
        Z.append(temp)
        ZAvg.append(numpy.average(temp))
        
    return Z, ZAvg