# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 18:20:07 2013

@author: chaco3
"""

import pyximport
pyximport.install()
import VariogramFit

import numpy
from numpy import linalg


def KrigInterp(ModOpt,POI,Loc,VarFunArr,xopt,k,CovMea,MeaRow):
    # Calculate distances between stations        
    SVr = numpy.array(CovMea)
    if linalg.det(CovMea) == 0:
        print('Non-singular covriance matrix - Sorry, cannot invert')
        return(9999,9999)
    InvSVr = linalg.inv(SVr)
    
    POID = numpy.zeros([len(POI),len(Loc)])
    SVm = numpy.zeros([len(POI),len(Loc)])
    
    for i in xrange(0,len(POI)):
        for j in xrange(0,len(Loc)):
            # Calculate distance from target to stations
            POID[i][j] = numpy.sqrt((POI[i][0]-Loc[j][0])**2 +
                                   (POI[i][1]-Loc[j][1])**2 )
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
    return Z,SP