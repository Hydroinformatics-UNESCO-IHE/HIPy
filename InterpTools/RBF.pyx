# -*- coding: utf-8 -*-
"""
Radial basis funcitions
Wrapper for scipy function

function can set to: 'multiquadric','inverse','linear','cubic','gaussian', or
'thin_plate'
"""

from scipy.interpolate import Rbf
import numpy

import DataLoad
import DataSave

def Interp(x, y, val, xt, yt, function='multiquadric', epsilon=3.0):
    cls = Rbf(x, y ,val , function=function, epsilon=epsilon)
    R = cls(xt,yt)
    return R
    
def Interp_bat(Loc, POI, Prec, tmin=0, tmax='def', valmin=-1e99, valmax=1e99,
               function='multiquadric', epsilon=3.0):
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
        temp = Interp(x, y, Prec[i], xt, yt, function, epsilon)
        numpy.clip(temp,valmin,valmax)
        Z.append(temp)
        ZAvg.append(numpy.average(temp))

    return Z, ZAvg
    
def simple_MRBF(SiteInfo, XYTargets, DataRecord, function='multiquadric',
                epsilon=3.0):
    Loc, POI, Prec = DataLoad.lcsv(SiteInfo, XYTargets, DataRecord)
    Z, ZAvg = Interp_bat(Loc, POI, Prec, function=function, epsilon=epsilon)
    DataSave.spkl(Z, 'PrecField')
    DataSave.spkl(ZAvg, 'AvgPrec')
    return 'Radial basis function, done!'
