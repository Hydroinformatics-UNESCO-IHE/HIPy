# -*- coding: utf-8 -*-
"""
===========
Interpolate
===========
Implemented by Juan Chacon @ UNESCO-IHE
Integrated Water Systems and Governance Department
Hydroinformatics Laboratory

This library provides functions for interpolation using different algorithms \
based on adjustment of regression functions. Machine learning regressors are \
not implemented in this tool. Cubic, linear and radial basis functions are \
based on wrappers of scipy libraries.

* Pre requisites
    you will not need libraries, not coming alongside with the \ 
    Anaconda ditribution (recommended).

* Functions
    * interpolate - performs spatial interpolation of a time series.
    
* Use policy
    * You should include the respective citation to the authors
    * If you find this tool usefull, you will give the main author a beer next\
    time you see him :)
    
* References
    * http://docs.scipy.org/doc/scipy-0.13.0/reference/interpolate.html
"""

#Modules to import
import sys
import os
sys.path.append(os.path.abspath('..\\Utilities'))

import pyximport
pyximport.install()

from scipy.interpolate import griddata, Rbf
import numpy
import dist

ERROR_CODE = 9999

def _cubic(x,y,val,xt,yt):
    '''
    Internal function for cubic interpolation of a single instance. This is a \
    wrapper of scipy function
    
    Parameters
    ----------
    x : ndarray, shape ``(n,1)``
        Vector with `x` coordinates of predictor variables.
    y : ndarray, shape ``(n,1)``
        Vector with `y` coordinates of predictor variables.
    val : ndarray, shape ``(n,1)``
        Vector with values of the predicted variable.
    xt : ndarray, shape ``(m,1)``
        Vector with the `x` location of targets to be interpolated.
    yt : ndarray, shape ``(m,1)``
        Vector with the `y` location of targets to be interpolated.
    Returns
    -------
    r : ndarray, shape ``(m,1)``
        Vector with the interpolation results at the target location
    '''
    r = griddata((x,y) ,val , (xt,yt), method='cubic')
    return r
    
def _idw(x,y,val,xt,yt,power):
    '''
    Internal function for inverse distance weighting interpolation of a single\
     instance
    
    Parameters
    ----------
    x : ndarray, shape ``(n,1)``
        Vector with `x` coordinates of predictor variables.
    y : ndarray, shape ``(n,1)``
        Vector with `y` coordinates of predictor variables.
    val : ndarray, shape ``(n,1)``
        Vector with values of the predicted variable.
    xt : ndarray, shape ``(m,1)``
        Vector with the `x` location of targets to be interpolated.
    yt : ndarray, shape ``(m,1)``
        Vector with the `y` location of targets to be interpolated.
    power : float
        Scalar for power interpolation law.
    Returns
    -------
    r : ndarray, shape ``(m,1)``
        vector with the interpolation results at the target location.
    '''
    vect = numpy.transpose(numpy.vstack((x,y)))
    vecttar = numpy.transpose(numpy.vstack((xt,yt)))
    dis = dist.target(vect, vecttar)

    r = []
    for i in xrange(0,len(xt)):
        if numpy.min(dis[i]) < 0.0001:
            r.append(val[numpy.argmin(dis)])
            
        else:
            d = numpy.array(dis[i])
            v = numpy.array(val)
            
            ds = numpy.power(1./d,power)
            dt = numpy.sum(ds)
            W = ds/dt
            r.append(numpy.sum(W*v))
    
    return r

def _linear(x,y,val, xt, yt):
    '''
    Internal function for linear interpolation of a single instance. This is a\
     wrapper of scipy function
    
    Parameters
    ----------
    x : ndarray, shape ``(n,1)``
        Vector with `x` coordinates of predictor variables.
    y : ndarray, shape ``(n,1)``
        Vector with `y` coordinates of predictor variables.
    val : ndarray, shape ``(n,1)``
        Vector with values of the predicted variable.
    xt : ndarray, shape ``(m,1)``
        Vector with the `x` location of targets to be interpolated.
    yt : ndarray, shape ``(m,1)``
        Vector with the `y` location of targets to be interpolated.
    Returns
    -------
    r : ndarray, shape ``(m,1)``
        Vector with the interpolation results at the target location
    '''
    r = griddata((x,y) ,val , (xt,yt), method='linear')
    return r

def _nearest(x,y,val,xt,yt):  
    '''
    Internal function for nearest neighbour interpolation of a single instance.
    
    Parameters
    ----------
    x : ndarray, shape ``(n,1)``
        Vector with `x` coordinates of predictor variables.
    y : ndarray, shape ``(n,1)``
        Vector with `y` coordinates of predictor variables.
    val : ndarray, shape ``(n,1)``
        Vector with values of the predicted variable.
    xt : ndarray, shape ``(m,1)``
        Vector with the `x` location of targets to be interpolated.
    yt : ndarray, shape ``(m,1)``
        Vector with the `y` location of targets to be interpolated.
    Returns
    -------
    r : ndarray, shape ``(m,1)``
        Vector with the interpolation results at the target location
    '''    
    vect = numpy.transpose(numpy.vstack((x,y)))
    vecttar = numpy.transpose(numpy.vstack((xt,yt)))
    dis = dist.target(vect, vecttar)
    cst = numpy.argmin(dis,1)
    
    r = []
    for i in xrange(0,len(xt)):
        r.append(val[cst[i]])    
    return r

def _rbf(x, y, val, xt, yt, rbf_function, rbf_epsilon):
    '''
    Internal function for radial basis function interpolation of a single \
    instance. This is a wrapper of scipy function
    
    Parameters
    ----------
    x : ndarray, shape ``(n,1)``
        Vector with `x` coordinates of predictor variables.
    y : ndarray, shape ``(n,1)``
        Vector with `y` coordinates of predictor variables.
    val : ndarray, shape ``(n,1)``
        Vector with values of the predicted variable.
    xt : ndarray, shape ``(m,1)``
        Vector with the `x` location of targets to be interpolated.
    yt : ndarray, shape ``(m,1)``
        Vector with the `y` location of targets to be interpolated.
    rbf_function : str or callable
        string with type of RBF function to be used as
        
            'multiquadric': sqrt((r/epsilon)**2 + 1)
            'inverse': 1.0/sqrt((r/epsilon)**2 + 1)
            'gaussian': exp(-(r/epsilon)**2)
            'linear': r
            'cubic': r**3
            'quintic': r**5
            'thin_plate': r**2 * log(r)
        
    rbf_epsilon : float
        Adjustable constant for gaussian or multiquadrics functions    
    Returns
    -------
    r : ndarray, shape ``(m,1)``
        Vector with the interpolation results at the target location
    '''
    cls = Rbf(x, y ,val , function=rbf_function, epsilon=rbf_epsilon)
    r = cls(xt,yt)
    return r
    
def interpolate(stations,targets,records,method='nearest',idw_power=2.0,
           rbf_function='multiquadric', rbf_epsilon=3.0,
           tmin=0, tmax='def',val_min=-1e99,val_max=1e99):
    '''
    public function for interpolation spatially distributed data series.
    
    Parameters
    ----------
    stations : array_like, shape ``(n,2)``
        Vector with `x,y` coordinates of measurement points.
    targets : array_like, shape ``(m,2)``
        Vector with `x,y` coordinates of the target points.
    records : ndarray, shape ``(p,n)``
        Matrix with the values of recorded variables, arranged by columns.
    method : str
        method to be used in the interpolation. Available methods are::  
            'cubic' : Interpolate using a 3rd degree polynomial.
            'linear' : Make linear interpolation between neighbouring stations.
            'nearest' : Uses the value of the closest stations.
            'rbf' : Make a interpolation of the type radial basis function.
            'idw' : Interpolate usign an inverse distance weighting approach.    
    yt : ndarray, shape ``(m,1)``
        Vector with the `y` location of targets to be interpolated.
    rbf_function : str or callable
        string with type of RBF function to be used as
            'multiquadric': sqrt((r/epsilon)**2 + 1)
            'inverse': 1.0/sqrt((r/epsilon)**2 + 1)
            'gaussian': exp(-(r/epsilon)**2)
            'linear': r
            'cubic': r**3
            'quintic': r**5
            'thin_plate': r**2 * log(r)
    rbf_epsilon : float
        Adjustable constant for gaussian or multiquadrics radial basis \
        functions
    tmin : float, optional
        Time step in which the dataseries start
    tmax : float, optional
        Time step in which the dataseries finishes
    val_min : float, optional
        Low cutoff of the interpolated variable
    val_max : float, optional
        high cutoff of the interpolated variable
    Returns
    -------
    z : ndarray, shape ``((tmax-tmin),m)``
        Matrix with interpolation results for each location at each time step.
    zavg : ndarray, shape ``((tmax-tmin),1)``
        Vector with the spatial average of the interpolated variable for each \
        time step.
    '''    
    if tmax == 'def':
        tmax = len(records)
    
    stations = numpy.array(stations)
    targets = numpy.array(targets)
    
    xt = targets[:,0]
    yt = targets[:,1]
    x = stations[:,0]
    y = stations[:,1]
    records = numpy.array(records)
    
    z = []
    zavg = []
    for i in xrange(tmin,tmax):
        if method is 'cubic':
            temp = _cubic(x,y,records[i],xt,yt)
        elif method is 'linear':
            temp = _linear(x,y,records[i],xt,yt)
        elif method is 'nearest':
            temp = _nearest(x,y,records[i],xt,yt)
        elif method is 'rbf':
            temp = _rbf(x,y,records[i],xt,yt,rbf_function,rbf_epsilon)
        elif method is 'idw':
            temp = _idw(x,y,records[i],xt,yt,idw_power,val_min)
        else:
            print 'Selected method is not valid'
            return ERROR_CODE, ERROR_CODE
        
        numpy.clip(temp,val_min,val_max)
        z.append(temp)
        zavg.append(numpy.average(temp))
    
    z = numpy.array(z)
    zavg = numpy.array(zavg)
    return z, zavg