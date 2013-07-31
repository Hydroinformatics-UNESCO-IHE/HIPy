# -*- coding: utf-8 -*-
"""
Linear Interpolation

Wrapper of scipy function
"""

from scipy.interpolate import griddata
import numpy

import DataLoad
import DataSave

def Interp(x, y, val, xt, yt):
    R = griddata((x,y) ,val , (xt,yt), method='linear')
    return R
    
def Interp_bat(Loc, POI, Prec, tmin=0, tmax='def', valmin=-1e99, valmax=1e99):
    if tmax == 'def':
        tmax = len(Prec)
    
    Loc = numpy.array(Loc)
    POI = numpy.array(POI)
    
    xt = POI[:,0]
    yt = POI[:,1]
    x = Loc[:,0]
    y = Loc[:,1]
    Prec = numpy.array(Prec)
    
    Z = []
    ZAvg = []
    for i in xrange(tmin,tmax):
        temp = Interp(x, y, Prec[i], xt, yt)
        numpy.clip(temp,valmin,valmax)
        Z.append(temp)
        ZAvg.append(numpy.average(temp))
        
    return Z, ZAvg
    
def simple_Linear(SiteInfo, XYTargets, DataRecord):
    Loc, POI, Prec = DataLoad.lcsv(SiteInfo, XYTargets, DataRecord)
    Z, ZAvg = Interp_bat(Loc, POI, Prec)
    DataSave.spkl(Z, 'PrecField')
    DataSave.spkl(ZAvg, 'AvgPrec')
    return 'Linear interpolation, done!'