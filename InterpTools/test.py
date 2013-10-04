# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 17:27:42 2013

@author: chaco3
"""

import pyximport
pyximport.install()
import KrigingP
import IDW
import Kriging
import DataLoad
import NearestNeigh
import Linear
import Cubic
import RBF
import numpy


def tIDW():
    '''
    IDW Test
    '''
    Loc, POI, Prec = DataLoad.lcsv('TestData\GaugeLoc.csv',
                                       'TestData\InterpPts.csv',
                                       'TestData\Dataset.csv')
    Z, Zavg = IDW.Interp_bat(Loc, POI, Prec, 2.0, 0.00001, 0, 20)
    return 'IDW working fine!'

def tKrig():
    '''
    Kriging test    
    '''
    Loc, POIC, Prec = DataLoad.lcsv('TestData\GaugeLoc.csv',
                                       'TestData\InterpPts.csv',
                                       'TestData\Dataset.csv')
    Loc = numpy.array(Loc)/1000.0
    POIC = numpy.array(POIC)/1000.0
    SVExp, CovMea = Kriging.exp_semivariogram(Prec, Loc)
    xopt, ModOpt, VarFunArr = Kriging.theor_variogram(SVExp)
    Z, SP, ZAvg = Kriging.Krig(10.0, POIC, Loc, Prec, CovMea, ModOpt, xopt, 
                               VarFunArr,10 ,11, 'Ord')
    print Z
    print ZAvg
    return Z, ZAvg, CovMea

def tKrigP():
    '''
    Kriging test    
    '''
    Loc, POIC, Prec = DataLoad.lcsv('TestData\GaugeLoc.csv',
                                       'TestData\InterpPts.csv',
                                       'TestData\Dataset.csv')
    Loc = numpy.array(Loc)/1000.0
    POIC = numpy.array(POIC)/1000.0
    SVExp, CovMea = KrigingP.exp_semivariogram(Prec, Loc)
    xopt, ModOpt, VarFunArr = KrigingP.theor_variogram(SVExp)
    Z, SP, ZAvg = KrigingP.Krig(10.0, POIC, Loc, Prec, CovMea, ModOpt, xopt, 
                               VarFunArr,10 ,11, 'Ord')
    print Z
    print ZAvg
    return 

def tNear():
    '''
    Nearest Neighborhood test    
    '''
    Loc, POI, Prec = DataLoad.lcsv('TestData\GaugeLoc.csv',
                                       'TestData\InterpPts.csv',
                                       'TestData\Dataset.csv')
    Z, Zavg = NearestNeigh.Interp_bat(Loc, POI, Prec, 0, 20)
    return 'Nearest neighborhood working fine!'
    
def tLinear():
    '''
    Nearest Neighborhood test    
    '''
    Loc, POI, Prec = DataLoad.lcsv('TestData\GaugeLoc.csv',
                                       'TestData\InterpPts.csv',
                                       'TestData\Dataset.csv')
    Z, Zavg = Linear.Interp_bat(Loc, POI, Prec, 0, 20)
    return 'Linear Interpolation working fine!'

def tCubic():
    '''
    Cubic interpolator test
    '''
    Loc, POI, Prec = DataLoad.lcsv('TestData\GaugeLoc.csv',
                                       'TestData\InterpPts.csv',
                                       'TestData\Dataset.csv')
    Z, Zavg = Cubic.Interp_bat(Loc, POI, Prec, 0, 20)
    return 'Cubic Interpolation working fine!'

def tRBF():
    '''
    Radial basis function test    
    '''
    Loc, POI, Prec = DataLoad.lcsv('TestData\GaugeLoc.csv',
                                       'TestData\InterpPts.csv',
                                       'TestData\Dataset.csv')
    Z, Zavg = RBF.Interp_bat(Loc, POI, Prec, 0, 20)
    return 'Radial Basis Function Interpolation working fine!'