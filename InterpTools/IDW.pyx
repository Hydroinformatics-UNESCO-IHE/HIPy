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
    
def load_data(SiteInfo, XYTargets, DataRecord):
    '''
    Load data from CSV files as presented in the description into working files
    inside the module. No data processing taking place here.    
    
    Input files description
    ---------------------------------------------------------------------------
     This script uses CSV (comma separated) files for import and ouput of data
     The necessary files are:
       -> SiteInfo.csv - File containing coordinates of sensors in the same
           order as the other files (This is mandatory).
          Format:
          [NAME, X, Y]
           
       -> XYTargets.csv - File containing the location of sample points inside
           the catchment that are going to be used to calculate the average of
           the given variable.
          Format:
           [X, Y, Catchment Number]
           
       -> DataRecord.csv - File containing the registers of the variable that
           is going to be interpolated. For historical reasons will be noted as
           Prec (precipitation).
          Format:
           [St1Data, St2Data, ..., StnData]
    '''
    Loc = []
    with open(SiteInfo,'rb') as SiteInfo:
        Lines = csv.reader(SiteInfo)
        Lines.next()
        for row in Lines:
            Loc.append([float(row[1]),float(row[2])])
    print 'Gauge Location, Imported'
    print '' 
    
    POI = []
    with open(XYTargets,'rb') as POIf:
        Lines = csv.reader(POIf)
        Lines.next()

        for row in Lines:
            POI.append([float(row[0]),float(row[1])])

    print 'Sampling Points, Imported'
    print ''    
            
    Prec = []
    with open(DataRecord,'rb') as Data:
        Lines = csv.reader(Data)
        Lines.next()
        for row in Lines:
            Prec.append([float(x) for x in row[:]])
    print 'Data, Imported'
    print ''
    
    if len(Prec[0]) != len(Loc):
        print 'Measurements and stations are not compatible. Please check your\
                data'
    return Loc, POI, Prec