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
     time you see him :)
    
* References

"""
f = 'initialise'
import numpy

def between(stations):
    '''
    Calculates the distance between stations.
    
    Parameters
    ----------
    stations : array_like, shape ``(n,2)``
        Vector with `x,y` coordinates of measurement stations
    
    Returns
    -------
    result_distance : array_like ``(n,n)``
        Matrix of size "[n,n]" with distance between stations
    '''
    # Calculate distances between stations     
    result_distance = numpy.zeros((len(stations),len(stations)))     
    for i in xrange(0,len(stations)):
        for j in xrange(0,len(stations)):
            result_distance[i][j] = numpy.sqrt((stations[i][0]-stations[j][0])**2 +
                                   (stations[i][1]-stations[j][1])**2 )
    return result_distance

def target(stations, targets):
    '''
    Calculates the distance between stations and targets
    
    Parameters
    ----------
    stations : array_like, shape ``(n,2)``
        Vector with `x,y` coordinates of measurement stations.
    targets : array_like, shape ``(m,2)``
        Vector with `x,y` coordinates of targets.
    
    Returns
    -------
    result_distance : ndarray, shape ``(n,m)``
        Matrix with distance between stations to targets.
    '''    
    # Calculate distances between stations and targetsgets
    result_distance = numpy.zeros((len(targets),len(stations)))      
    for i in xrange(0,len(targets)):
        for j in xrange(0,len(stations)):
            # Calculate distance from target to stations
            result_distance[i][j] = numpy.sqrt((targets[i][0]-stations[j][0])**2 +
                                   (targets[i][1]-stations[j][1])**2 )
    return result_distance