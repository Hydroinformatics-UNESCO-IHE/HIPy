# -*- coding: utf-8 -*-
"""
Data Load Module

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
"""
import csv
import cPickle

def lcsv(SiteInfo, XYTargets, DataRecord):
    '''
    Load data from CSV files as presented in the description into working files
    inside the module. No data processing taking place here.
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

def scsv(SiteInfo):
    '''
    Load data from CSV files as presented in the description into working files
    inside the module. No data processing taking place here.
    '''
    Loc = []
    with open(SiteInfo,'rb') as SiteInfo:
        Lines = csv.reader(SiteInfo)
        Lines.next()
        for row in Lines:
            Loc.append([float(row[0]),float(row[1])])
    print 'Variable Imported'
    print ''
    return Loc
    
def lpkl(SiteInfo, XYTargets, DataRecord):
    '''
    Load data from PKL files as presented in the description into working files
    inside the module. No data processing taking place here. PKL is faster to
    load than CSV in python.
    '''
    with open(SiteInfo,'rb') as SiteInfo:
        Loc = cPickle.load(SiteInfo)
    print 'Gauge Location, Imported'
    print '' 
    
    with open(XYTargets,'rb') as POIf:
        POI = cPickle.load(POIf)
    print 'Sampling Points, Imported'
    print ''    
            
    with open(DataRecord,'rb') as Data:
        Prec = cPickle.load(Data)
    print 'Data, Imported'
    print ''
    
    if len(Prec[0]) != len(Loc):
        print 'Measurements and stations are not compatible. Please check your\
                data'
    return Loc, POI, Prec
    
def slpkl(Datafile):
    '''
    Load data from PKL files as presented in the description into working files
    inside the module. No data processing taking place here. PKL is faster to
    load than CSV in python.
    '''
    with open(Datafile,'r') as foo:
        VarName = cPickle.load(foo)
    print 'Variable Imported'
    print ''
    return VarName