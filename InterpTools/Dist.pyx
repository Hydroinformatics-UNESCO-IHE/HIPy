# -*- coding: utf-8 -*-
"""
Distance Calculator
"""
import numpy

def between(Loc):
    # Calculate distances between stations     
    POID = []    
    for i in xrange(0,len(Loc)):
        for j in xrange(0,len(Loc)):
            POID[i][j] = numpy.sqrt((Loc[i][0]-Loc[j][0])**2 +
                                   (Loc[i][1]-Loc[j][1])**2 )
    return POID

def target(Loc, Tar):
    # Calculate distances between stations and targets
    POID = []    
    for i in xrange(0,len(Tar)):
        for j in xrange(0,len(Loc)):
            # Calculate distance from target to stations
            POID[i][j] = numpy.sqrt((Tar[i][0]-Loc[j][0])**2 +
                                   (Tar[i][1]-Loc[j][1])**2 )
    return POID