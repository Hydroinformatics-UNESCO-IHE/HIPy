# -*- coding: utf-8 -*-
"""
Nearest Neighborhood
"""

import numpy
import Dist
import DataLoad
import DataSave


def Interp(x,y,val,xt,yt):  
    
    vect = numpy.transpose(numpy.vstack((x,y)))
    vecttar = numpy.transpose(numpy.vstack((xt,yt)))
    dis = Dist.target(vect, vecttar)
    cst = numpy.argmin(dis,1)
    
    R = []
    for i in xrange(0,len(xt)):
        R.append(val[cst[i]])
    
    return R

def Interp_bat(Loc, POI, Prec, tmin=0, tmax='def'):
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
        temp = Interp(x,y,Prec[i],xt,yt)        
        Z.append(temp)
        ZAvg.append(numpy.average(temp))
    
    return Z, ZAvg

def simple_NN(SiteInfo, XYTargets, DataRecord):
    Loc, POI, Prec = DataLoad.lcsv(SiteInfo, XYTargets, DataRecord)
    Z, ZAvg = Interp_bat(Loc, POI, Prec)
    
    DataSave.spkl(Z, 'PrecField')
    DataSave.spkl(ZAvg, 'AvgPrec')
    return 'Nearest Neighborhood done!'