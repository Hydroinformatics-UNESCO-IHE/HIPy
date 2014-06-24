# -*- coding: utf-8 -*-
"""
=====================
Kriging Interpolation
=====================
Implemented by Juan Chacon @ UNESCO-IHE
Integrated Water Systems and Governance Department
Hydroinformatics Laboratory

This library simple and ordinary Kriging interpolation. Other Kriging \
applications such as universal kriging or universal Kriging, will be \
implemented in a posterior stage.

* Pre requisites
    you will need the following libraries, not coming alongside with the \ 
    Anaconda ditribution (recommended)

    * pyOpt - Optimisation engine, used to solve the semivariogram fitting  

* Functions
    * exp_semivariogram: Computes experimental semivariogram
    * theor_variogram: Adjust theoretical to experimental semivariogram
    * kriging_core: Performs data parsing and preprocessing for kriging \
    interpolation
    * krig: Solves the Kriging system and returns location-wise estimates of \
    points of interest
    * simple_Krig: Simple inteface for Kriging interpolation, by reading \
    standard csv files of data, and generating pkl of interpolation results
    * test: Test for Kriging module to see if everything is running as it \
    should

* Use policy
    * You should include the respective citation to the authors
    * If you find this tool usefull, you will give the main author a beer next\
    time you see him :)
    
* References
    * http://people.ku.edu/~gbohling/cpe940/Kriging.pdf
"""
#Libraries to be imported
#------------------------------------------------------------------------------
import pyximport
pyximport.install()

import sys
import os
sys.path.append(os.path.abspath('..\\Utilities'))

import variogram_fit
import dist
import numpy
from numpy import linalg
import random
from pyOpt import ALHSO, Optimization
import data_save
import data_load
import time


ERROR_CODE = 9999
#------------------------------------------------------------------------------
def _regul(lag,cov,min_bin,maxdist):
    '''
    Internal function that regularises data for semivariogram computation by \
    selecting a minimum number of bins or a maximum distance criteria for bin \
    size.
    
    Parameters
    ----------
    lag : array_like
        Variable holding the lags for regularisation 
    cov : array_like
        contains the data which is going to be averaged in the bin
    min_bin : int
        minimum number of bins of the output
    maxdist : float
        maximum lag distance of each bin
    
    Returns
    -------
    lag2 : array_like
        vector holding the regularised bin location
    cov2 : array_like
        vector holding the regularised average of the variable in the bin
    '''
    maxdif = numpy.max(lag) #always starting from 0
    num_bin_dist = (maxdif / maxdist)
    
    if num_bin_dist > min_bin:
        total_number_bins = num_bin_dist
    else:
        total_number_bins = min_bin

    dist_between_bins = maxdif / total_number_bins
    
    lag2 = []
    cov2 = []
    for i in xrange(0,int(total_number_bins)):
        #get indices of elements within the bin.
        indx = [k for k, x in enumerate(lag) if i*dist_between_bins < x <= 
                (i+1)*dist_between_bins]
        sumbin = 0        
        if indx != []:
            for j in xrange(0,len(indx)):
                sumbin += cov[indx[j]]  
            cov2.append(sumbin/len(indx))
            lag2.append(i*dist_between_bins + dist_between_bins/2)
    return lag2,cov2
    
def exp_semivariogram(records, stations):
    '''
    Public function for establishing the experimental semivariogram.
    
    Parameters
    ----------
    records : array_like, shape ``(p,n)``
        Vector for which semivariogram is going to be calculated.
    stations : array_like shape ``(n,2)``
        Vector with the `x,y` coordinates of the measurement stations.
    
    Returns
    -------
    experimental_sv : array_like, shape ``(n,n)``
        Experimental semivariogram vector composed of lag and semivariogram
    record_covariance_matrix : array_like, shape ``(n,n)``
        Covariance matrix for the recorded data. 
    '''
    ## Removal of no precipitation events   
    WetMeas = []
    for i in xrange(0,len(records)):
        if numpy.max(records[i]) > 3:
            WetMeas.append(records[i])    
    
    ## Measurement covariance
    record_covariance_matrix = numpy.cov(numpy.transpose(WetMeas))
    Dis = dist.between(stations)
    
    ## Experimental Semivariogram
    experimental_sv = []
    for i in xrange(0,len(record_covariance_matrix)-1):
        for j in xrange(i+1,len(record_covariance_matrix)):
            Cov = record_covariance_matrix[i][j]
            Lag = Dis[i][j]
            experimental_sv.append([Lag,Cov])
            
    experimental_sv = numpy.array(experimental_sv)
    Lag2, Cov2 = _regul(experimental_sv[:,0],experimental_sv[:,1],15,1)
    experimental_sv = numpy.transpose(numpy.vstack((Lag2,Cov2)))
    return experimental_sv, record_covariance_matrix

def theor_variogram(experimental_sv,Sb=(0.01,400),Rb=(2,20),Nb=(0,400),
                    ab=(0,2), vb=(0,1000),candidate_sv='def',
                    candidate_sv_tag='def'):
    '''
    Fitting of theoretical variogram
    Parameters
    ----------
        **experimental_sv** -- Experimental semivariogram ''[x,2]'', lag and \
            semivariogram \n
        **Sb** -- Boundaries on Sill of semivariogram ``(min,max)`` \n
        **Rb** -- Boundaries on Range of semivariogram ``(min,max)`` \n
        **Nb** -- Boundaries on Nugget of semivariogram ``(min,max)`` \n
        **ab** -- Boundaries on Power of semivariogram ``(min,max)`` (only \
            valid for power semivariogram) \n
        **vb** -- Boundaries on Shape parameter of semivariogram ``(min,max)``\
            (only valid for matérn type) \n
    
    Returns
    -------
        **xopt** -- Vector with optimal semivariogram parameters ``[5]`` \n
        **ModOpt** -- Pointer to optimal vector location \n
        **candidate_sv** -- Array with pointer to functions in variogram_fit module
    '''                      
    
    if candidate_sv is 'def':
        # Array with functions to be called from the Variograms library
        candidate_sv = [variogram_fit.SVExponential, variogram_fit.SVGaussian, 
                     variogram_fit.SVSpherical, variogram_fit.SVCubic,
                     variogram_fit.SVPentaspherical, variogram_fit.SVSinehole, 
                     variogram_fit.SVPower]
     
    if candidate_sv_tag is 'def':
    # Names of functions for display only
        candidate_sv_tag = ['Exponential','Gaussian','Spherical','Cubic',
                     'Pentaspherical','Sinehole','Power','Matern']
    
    # Initial seed for variogram fit
    sr = random.uniform(Sb[0],Sb[1])
    rr = random.uniform(Rb[0],Rb[1])
    nr = random.uniform(Nb[0],Nb[1])
    ar = random.uniform(ab[0],ab[1])
    vr = random.uniform(vb[0],vb[1])
    
    Var = []
    Res = []
    Mdl = [] 
    
    # Wrapper of minimisation function (RMSE) for semivariogram fitting
    def _opt_fun(x,*args):
        F, g, fail = variogram_fit.fit_function(x,experimental_sv,j,candidate_sv)
        if F == ERROR_CODE:
            fail = 1

        else:
            Var.append(x)
            Res.append(F)
            Mdl.append(j)
        return F, g, fail
    
    # Optimisation starts to minimise differences between experimental and 
    # theoretical semivariograms
    for j in xrange(0,len(candidate_sv)):   
        VarProb = Optimization('Variogram Fitting: ' + candidate_sv_tag[j], _opt_fun)
        VarProb.addObj('RMSE')
        VarProb.addVar('Sill','c',lower=Sb[0],upper=Sb[1],value=sr)
        VarProb.addVar('Range','c',lower=Rb[0],upper=Rb[1],value=rr)
        VarProb.addVar('Nugget','c',lower=Nb[0],upper=Nb[1],value=nr)
        VarProb.addVar('Exponent (a)','c',lower=ab[0],upper=ab[1],value=ar)
        VarProb.addVar('Rank (v)','c',lower=vb[0],upper=vb[1],value=vr)
        
        args = (experimental_sv, j, candidate_sv, Var, Res, Mdl)
        optmz = ALHSO()
        optmz.setOption('fileout',0)
        optmz(VarProb)

    # Get pointer to best semivariogram
    k = numpy.argmin(Res)
    xopt = Var[k]
    ModOpt = Mdl[k]
    
    return xopt, ModOpt, candidate_sv

def _kriging_core(ModOpt,single_target,stations,candidate_sv,xopt,record_covariance_matrix,records,krig_type):
    '''
    Kriging core where interpolating algorithms is taking place.
    
    Parameters
    ----------
        **ModOpt** -- Pointer to optimal semivariogram model \n
        **single_target** -- Single point of interest to interpolate (targets)\
            ``[t,2]`` \n
        **stations** -- Gauge location for interpolation ``[x,2]`` \n
        **candidate_sv** -- Array with pointer to functions in variogram_fit module
            \n
        **xopt** -- vector with optimal semivariogram parameters ``[5]`` \n
        **record_covariance_matrix** -- Gauge records covariance matrix ``[x,x]`` \n
        **records** -- Precipitation register for gauges ``[n,x]`` \n
        **krig_type** -- Type of Kriging to be used. 'Sim' for Simple and 'Ord'\
            for Ordinary Kriging
    
    Returns
    -------
        **Z** -- Interpolation for each target and time step ``[n,1]`` \n
        **SP** -- Interpolation variance field ``[1]``\n
    '''
    targetsD = dist.target(stations,[single_target])[0]
    SVm = []
    for j in xrange(0,len(stations)):
        SVm.append(candidate_sv[ModOpt](targetsD[j],xopt))

    if krig_type is 'Ord':
    #Ordinary Kriging
        record_covariance_matrix = numpy.row_stack((record_covariance_matrix,numpy.ones(len(record_covariance_matrix))))
        record_covariance_matrix = numpy.column_stack((record_covariance_matrix,numpy.ones(len(record_covariance_matrix))))
        record_covariance_matrix[-1,-1] = 0.0
        SVm.append(1.0)
        
        SVr = numpy.array(record_covariance_matrix)
        if linalg.det(record_covariance_matrix) == 0:
            print('Non-singular covriance matrix - Sorry, cannot invert \n')
            return ERROR_CODE*numpy.ones(len(records)), ERROR_CODE*numpy.ones(len(records))
        InvSVr = linalg.inv(SVr)

        Z = []        
        WM= numpy.dot(InvSVr,SVm)
        for i in xrange(0,len(records)):
            Ztemp = numpy.dot(WM[:-1],records[i])
            Ztemp = numpy.clip(Ztemp,0,max(records[i])) # cutoff at 0 and max prec
            Z.append(Ztemp)        
#        S = xopt[0]*numpy.ones(len(SVm)-1)-SVm[:-1]
        S = SVm[:-1]
        SP = xopt[0] - (numpy.dot(WM[:-1],numpy.transpose(S))) - WM[-1]
        
    elif krig_type is 'Sim':
    #Simple Kriging
        SVr = numpy.array(record_covariance_matrix)
        if linalg.det(record_covariance_matrix) == 0:
            print('Non-singular covriance matrix - Sorry, cannot invert \n')
            return ERROR_CODE*numpy.ones(len(records)), ERROR_CODE*numpy.ones(len(records))
        InvSVr = linalg.inv(SVr)
        
        Z = []        
        WM= numpy.dot(InvSVr,numpy.transpose(numpy.array(SVm)))    
        for i in xrange(0,len(records)):
            Ztemp = numpy.dot(WM,records[i])
            Ztemp = numpy.clip(Ztemp,0,max(records[i])) # cutoff at 0 and max prec
            Z.append(Ztemp)        
        S = xopt[0]*numpy.ones(len(SVm))-SVm
        SP = xopt[0] - (numpy.dot(WM,numpy.transpose(S)))
    
    else:
        print 'I pity the fool for no chosing Kriging type'
        print 'only available Ord and Sim \n'
        Z = ERROR_CODE*numpy.ones(len(records))
        SP = ERROR_CODE*numpy.ones(len(records))
    
    return Z,SP

def krig(MaxDist, targets, stations, records, record_covariance_matrix, ModOpt, xopt,
         candidate_sv, tmin=0, tmax='def', MinNumSt=3, krig_type='Sim'):
    '''
    Parsing of data for Kriging interpolation. This function contains\
    execution of interpolation as well. \n
    Parameters
    ----------
        **MaxDist** -- Initial search radius for nearby stations \n
        **targets** -- Points of interest to interpolate (targets) ``[t,2]`` \n
        **stations** -- Gauge location for interpolation ``[x,2]`` \n
        **records** -- Precipitation register for gauges ``[n,x]`` \n
        **record_covariance_matrix** -- Gauge records covariance matrix ``[x,x]`` \n
        **ModOpt** -- pointer to optimal semivariogram model \n
        **xopt** -- vector with optimal semivariogram parameters ``[5]`` \n
        **candidate_sv** -- Array with pointer to functions in variogram_fit module
            \n
        **tmin** -- Initial time step to be interpolated \n
        **tmax** -- Final time step to be interpolated \n
        **MinNumSt** -- Minimum number of stations within the search radius \n
        **krig_type** -- Type of Kriging to be used. 'Sim' for Simple and 'Ord'\
            for Ordinary Kriging
    
    Returns
    -------
        **Z** -- Interpolation for each target and time step ``[n,t]`` \n
        **SP** -- Interpolation variance field ``[t,1]``\n
        **ZAvg** -- Average of interpolated field ``[n,1]``
    '''
    if tmax == 'def':
        tmax = len(records)
    tmin = int(tmin)
    tmax = int(tmax)
    PrecSec = records[tmin:tmax]
    
    # Reduce measurements to relevant locations for the targets
    Z = []
    SP = []
    for kk in xrange(0,len(targets)):
        # Check if there are enough stations for interpolation, otherwise,
        # increase search radius
        targets_dt = dist.target(stations, [targets[kk]])[0]
        TNS = 0
        MaxDist2 = MaxDist
        while TNS <= MinNumSt:
            #Minimum distance
            selected_stations = [i for i,v in enumerate(targets_dt) if v > MaxDist2]
            TNS = len(stations)-len(selected_stations)
            MaxDist2 += 1.0
        
        # Reduction of relevant stations (reduced data and cov matrices)            
        RedLoc = numpy.delete(stations,selected_stations,0)
        reduced_records = numpy.delete(PrecSec,selected_stations,1)    
        reduced_cov_matrix = record_covariance_matrix[:]    
        reduced_cov_matrix = numpy.delete(reduced_cov_matrix,selected_stations,0)
        reduced_cov_matrix = numpy.delete(reduced_cov_matrix,selected_stations,1)
        
        # Kriging interpolation
        TempRes = _kriging_core(ModOpt,targets[kk],RedLoc,candidate_sv,xopt,reduced_cov_matrix,reduced_records,krig_type)
        if Z == []:
            Z = numpy.vstack(TempRes[0])
        else:
            temp = numpy.vstack(TempRes[0])
            Z = numpy.hstack((Z,temp))
        SP.append(TempRes[1])

    ZAvg = numpy.average(Z,1)
    return Z, SP, ZAvg

def simple_Krig(SiteInfo, XYTargets, DataRecord):
    '''
    Wrapper for Kriging interpolation of all data, and save in PKL format \n
    Parameters 
    ----------
        **SiteInfo** -- Path to file with gauge location \n
        **XYTargets** -- Path to file with interpolation target locations \n
        **DataRecord** -- Path to file with variable registries \n

    Returns
    -------
        **VarField** -- file with pickled variable field \n
        **VarUnc** -- file with pickled variance of estimated variable \n
        **AvgVar** -- file with pickled average of estimated field
    '''
    
    stations, targets, records = data_load.lcsv(SiteInfo, XYTargets, DataRecord)
    experimental_sv, record_covariance_matrix = exp_semivariogram(records, stations)
    xopt, ModOpt, candidate_sv = theor_variogram(experimental_sv)
    Z, SP, ZAvg = krig(xopt[0]/3, targets, stations, records, record_covariance_matrix, ModOpt, xopt,
         candidate_sv, tmin=0, tmax='def', MinNumSt=3, krig_type='Sim')

    return Z, SP, ZAvg

def test():
    '''
    Module testing function
    '''
    stations, targets, records = data_load.lcsv('TestData\GaugeLoc.csv',
                                       'TestData\InterpPts.csv',
                                       'TestData\Dataset.csv')
    stations = numpy.array(stations)/1000.0
    targets = numpy.array(targets)/1000.0
    experimental_sv, record_covariance_matrix = exp_semivariogram(records, stations)
    xopt, ModOpt, candidate_sv = theor_variogram(experimental_sv)
    
    Z, SP, ZAvg = krig(10.0, targets, stations, records, record_covariance_matrix, ModOpt, xopt, 
                               candidate_sv,10 ,20, 'Ord')
    print 'Ordinary Kriging working fine'
    
    Z, SP, ZAvg = krig(10.0, targets, stations, records, record_covariance_matrix, ModOpt, xopt, 
                               candidate_sv,10 ,20, 'Sim')
    print 'Simple Kriging working fine'
    return 'Ended normally, module should be working properly'