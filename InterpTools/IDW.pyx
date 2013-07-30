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
import Dist
import csv
import DataSave
import DataLoad

def Interp(x, y, val, xt, yt, power=2.0, cutoff=0.00001):
    '''
    Algorithm for IDW interpolation. This contains the main build to perform
    interpolation using any type of IDW procedure.
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

def Interp_bat(Loc, POI, Prec, power=2.0, cutoff=0.00001, tmin=0, tmax='def'):
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
        Z.append(temp)
        ZAvg.append(numpy.average(temp))
        
    return Z, ZAvg

def simple_IDW(SiteInfo, XYTargets, DataRecord):
    Loc, POI, Prec = DataLoad.lcsv(SiteInfo, XYTargets, DataRecord)
    Z, ZAvg = Interp_bat(Loc, POI, Prec)
    DataSave.spkl(Z, 'PrecField')
    DataSave.spkl(ZAvg, 'AvgPrec')
    return 'IDW done'