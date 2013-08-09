# -*- coding: utf-8 -*-
"""
TEST2
"""

import pyximport
pyximport.install()
import Kriging
#import cPickle
import VariogramFit
import numpy
import DataLoad

#Loc = numpy.array([[1,1],[2,2],[4,2]])
#POI = numpy.array([[3,1],[2,3],[4,5]])
#MeaRow = numpy.array([4,5,6])
#Prec = numpy.array([[4,5,6],[4,5,6]])
#Loc = numpy.transpose(Loc)
#ModOpt = 2
#CovMea = [[1.0, 0.7, 0.3],
#          [0.7, 1.0, 0.5],
#          [0.3, 0.5, 1.0]]

#VarFunArr = [VariogramFit.SVExponential, VariogramFit.SVGaussian, 
#             VariogramFit.SVSpherical, VariogramFit.SVCubic,
#             VariogramFit.SVPentaspherical, VariogramFit.SVSinehole, 
#             VariogramFit.SVPower, VariogramFit.SVMatern]
#
#xopt = [1.2, 3.0, 0.0, 1, 1]
#MaxDist = 5
#MinNumSt = 3

Loc, POIC, Prec = DataLoad.lcsv('TestData\GaugeLoc.csv',
                                   'TestData\InterpPts.csv',
                                   'TestData\Dataset.csv')
Loc = numpy.array(Loc)/1000.0
POIC = numpy.array(POIC)/1000.0
SVExp, CovMea = Kriging.exp_semivariogram(Prec, Loc)
xopt, ModOpt, VarFunArr = Kriging.theor_variogram(SVExp)
Z, SP, ZAvg = Kriging.Krig(10.0, POIC, Loc, Prec, CovMea, ModOpt, xopt, VarFunArr,10 ,11,xopt[0],'Ord')
print Z
print ZAvg
