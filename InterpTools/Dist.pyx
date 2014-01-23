# -*- coding: utf-8 -*-
"""
===================
Distance Calculator
===================
Implemented by Juan Chacon @ UNESCO-IHE
Integrated Water Systems and Governance Department
Hydroinformatics Laboratory

This library contains utilities to calculate euclidean distances between \
stations, and with respect to a pre-defined target.

* Pre requisites
    you will need the following libraries, not coming alongside with the \ 
    Anaconda ditribution (recommended)

* Functions
    *between: Calculates the geometric distance between a series of points.
    *target: Calculate the geometric distance to a target, from different \ 
    points. 

* Use policy
    * You should include the respective citation to the authors
    * If you find this tool usefull, you will give the main author a beer next\
     time you see him/her :)
    
* References

"""
import numpy

def between(Loc):
    '''
    Parameters
    ----------
        **Loc -- Vector of size "[n,2]" containing "[x,y]" data pairs with \
        location of sensors
    
    Returns
    -------
        **POID -- Matrix of size "[n,n]" with distance between stations
    '''
    # Calculate distances between stations     
    POID = numpy.zeros((len(Loc),len(Loc)))     
    for i in xrange(0,len(Loc)):
        for j in xrange(0,len(Loc)):
            POID[i][j] = numpy.sqrt((Loc[i][0]-Loc[j][0])**2 +
                                   (Loc[i][1]-Loc[j][1])**2 )
    return POID

def target(Loc, Tar):
    '''
    Parameters
    ----------
        **Loc -- Vector of size "[n,2]" containing "[x,y]" the coordinates of \
        the location of the sensors
        ** Tar -- Vector of size "[m,2]" containing the targets coordinates
    
    Returns
    -------
        **POID -- Matrix of size "[n,m]" with distance between stations to \
        targets
    '''    
    # Calculate distances between stations and targets
    POID = numpy.zeros((len(Tar),len(Loc)))      
    for i in xrange(0,len(Tar)):
        for j in xrange(0,len(Loc)):
            # Calculate distance from target to stations
            POID[i][j] = numpy.sqrt((Tar[i][0]-Loc[j][0])**2 +
                                   (Tar[i][1]-Loc[j][1])**2 )
    return POID