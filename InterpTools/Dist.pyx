# -*- coding: utf-8 -*-
"""
Distance Calculator
"""
def DistBet(ndarray[numpy.float, ndim=2] Loc):
    # Calculate distances between stations     
    for i in xrange(0,len(Loc)):
        for j in xrange(0,len(Loc)):
            POID[i][j] = numpy.sqrt((Loc[i][0]-Loc[j][0])**2 +
                                   (Loc[i][1]-Loc[j][1])**2 )

def DistTar(ndarray[numpy.float, ndim=2] Loc, 
            ndarray[numpy.float, ndim=1] Tar):
    # Calculate distances between stations and targets
    for i in xrange(0,len(Tar)):
        for j in xrange(0,len(Loc)):
            # Calculate distance from target to stations
            POID[i][j] = numpy.sqrt((Tar[i][0]-Loc[j][0])**2 +
                                   (Tar[i][1]-Loc[j][1])**2 )