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
from numpy import linalg
import Dist

Loc = numpy.array([[1,1],[2,2],[4,2]])
POI = numpy.array([[3,1],[2,3],[4,5]])
MeaRow = numpy.array([4,5,6])
Prec = numpy.array([[4,5,6],[4,5,6]])
#Loc = numpy.transpose(Loc)
ModOpt = 2
CovMea = [[1.0, 0.7, 0.3],
          [0.7, 1.0, 0.5],
          [0.3, 0.5, 1.0]]

VarFunArr = [VariogramFit.SVExponential, VariogramFit.SVGaussian, 
             VariogramFit.SVSpherical, VariogramFit.SVCubic,
             VariogramFit.SVPentaspherical, VariogramFit.SVSinehole, 
             VariogramFit.SVPower, VariogramFit.SVMatern]

xopt = [1.2, 3.0, 0.0, 1, 1]
MaxDist = 5
MinNumSt = 1

SVr = numpy.array(CovMea)
if linalg.det(CovMea) == 0:
    print('Non-singular covriance matrix - Sorry, cannot invert')

InvSVr = linalg.inv(SVr)

POID = Dist.target(Loc,POI)
SVm = numpy.zeros([len(POI),len(Loc)])    
for i in xrange(0,len(POI)):
    for j in xrange(0,len(Loc)):
        # Calculate variance with between observation and targets
        SVm[i][j] = VarFunArr[ModOpt](POID[i][j],xopt)

WM = []
Z = []
SP = []
for i in xrange(0,len(SVm)):
    WM.append(numpy.dot(InvSVr,numpy.transpose(numpy.array(SVm[i]))))
    Ztemp = numpy.dot(WM[i],MeaRow)
    Z.append(numpy.clip(Ztemp,0,max(MeaRow))) ## cutoff at 0 and max prec
    S = xopt[0]*numpy.ones(len(SVm[i]))-SVm[i]
    SP.append(numpy.dot(WM[i],numpy.transpose(S)))
print WM
Z2, SP2 = Kriging.Kriging_core(ModOpt,POI,Loc,VarFunArr,xopt,CovMea,MeaRow)

Z3, SP3, ZAvg3 = Kriging.Krig(MaxDist, POI, Loc, Prec, CovMea, ModOpt, xopt,
         VarFunArr, 0, 1, 2)
print Z
print Z2
print Z3