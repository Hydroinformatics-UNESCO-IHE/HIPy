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
    * Kriging_core: Performs data parsing and preprocessing for kriging \
    interpolation
    * Krig: Solves the Kriging system and returns location-wise estimates of \
    points of interest
    * simple_Krig: Simple inteface for Kriging interpolation, by reading \
    standard csv files of data, and generating pkl of interpolation results
    * test: Test for Kriging module to see if everything is running as it \
    should

* Use policy
    * You should include the respective citation to the authors
    * If you find this tool usefull, you will give the main author a beer next\
    time you see him/her :)
    
* References
    * http://people.ku.edu/~gbohling/cpe940/Kriging.pdf
"""
#Libraries to be imported
#------------------------------------------------------------------------------
import pyximport
pyximport.install()
import VariogramFit
import Dist
import scipy.special
import numpy
from numpy import linalg
import random
from pyOpt import ALHSO, Optimization
import DataSave
import DataLoad
import time

#------------------------------------------------------------------------------
def Regul(lag,cov,minbin,maxdist):
    '''
    Regularises data for semivariogram computation by selecting a minimum \
    number of bins or a maximum distance criteria for bin size.
    
    Parameters
    ----------
        **lag -- Variable holding the lags for regularisation \n
        **cov -- contains the data which is going to be averaged in the bin \n
        **minbin -- minimum number of bins of the output \n
        **maxdist -- maximum lag distance of each bin \n
    
    Returns
    -------
        **lag2 -- vector holding the regularised bin location \n
        **cov2 -- vector holding the regularised average of the variable in \
        the bin
    '''
    maxdif = numpy.max(lag) #always starting from 0
    numbind = (maxdif / maxdist)
    
    if numbind > minbin:
        numtotbin = numbind
    else:
        numtotbin = minbin

    distbetbin = maxdif / numtotbin
    lag2 = []
    cov2 = []
    for i in xrange(0,int(numtotbin)):
        #get indices of elements within the bin.
        indx = [k for k, x in enumerate(lag) if i*distbetbin < x <= 
                (i+1)*distbetbin]
        sumbin = 0        
        if indx != []:
            for j in xrange(0,len(indx)):
                sumbin += cov[indx[j]]  
            cov2.append(sumbin/len(indx))
            lag2.append(i*distbetbin + distbetbin/2)
    return lag2,cov2
    
def exp_semivariogram(Prec, Loc):
    '''
    Experimental semivariogram computation
    
    Parameters
    ----------
        **Prec** -- Value vector ``[n,x]`` for which semivariogram is going to\
         be calculated \n
        **Loc** -- Vector of size ``[x,2]`` in which are stored the *x* and \ 
        *y* location
    
    Returns
    -------
        **SVExp** -- Experimental semivariogram vector ``[x,2]`` with lag and \
        semivariogram \n
        **CovMea** -- Covariance matrix ``[x,x]`` between stations 
    '''
    ## Removal of no precipitation events   
    WetMeas = []
    for i in xrange(0,len(Prec)):
        if numpy.max(Prec[i]) > 3:
            WetMeas.append(Prec[i])    
    
    ## Measurement covariance
    CovMea = numpy.cov(numpy.transpose(WetMeas))
#    if linalg.det(CovMea) == 0:
##        print 'Warning, Covariance Matrix is non-signular. Check your data \n'
#    else:
#        print 'Determinant Covariance Matrix', str(linalg.det(CovMea)), '\n'
#    print ''

    Dis = Dist.between(Loc)
#    print 'Distance Matrix - Done \n'
    
    ## Experimental Semivariogram
    SVExp = []
    for i in xrange(0,len(CovMea)-1):
        for j in xrange(i+1,len(CovMea)):
            Cov = CovMea[i][j]
            Lag = Dis[i][j]
            SVExp.append([Lag,Cov])
            
#    print 'Experimental semivariogram - Done \n'
    SVExp = numpy.array(SVExp)
    Lag2, Cov2 = Regul(SVExp[:,0],SVExp[:,1],15,1)
    SVExp = numpy.transpose(numpy.vstack((Lag2,Cov2)))
    return SVExp, CovMea

def theor_variogram(SVExp, Sb=(0.01,400), Rb=(2,20), Nb=(0,400),
                  ab=(0,2), vb=(0,1000),VarFunArr='def',optFunNam='def'):
    '''
    Fitting of theoretical variogram
    Parameters
    ----------
        **SVExp** -- Experimental semivariogram ''[x,2]'', lag and \
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
        **VarFunArr** -- Array with pointer to functions in VariogramFit module
    '''                      
    
    if VarFunArr is 'def':
        # Array with functions to be called from the Variograms library
        VarFunArr = [VariogramFit.SVExponential, VariogramFit.SVGaussian, 
                     VariogramFit.SVSpherical, VariogramFit.SVCubic,
                     VariogramFit.SVPentaspherical, VariogramFit.SVSinehole, 
                     VariogramFit.SVPower]
     
    if optFunNam is 'def':
    # Names of functions for display only
        optFunNam = ['Exponential','Gaussian','Spherical','Cubic',
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
    def OptFun(x,*args):
        F, g, fail = VariogramFit.optFunMaster(x,SVExp,j,VarFunArr)
        if F == 9999:
            fail = 1

        else:
            Var.append(x)
            Res.append(F)
            Mdl.append(j)
        return F, g, fail
        
#    print 'Initialising Variogram fit \n'
    
    # Optimisation starts to minimise differences between experimental and 
    # theoretical semivariograms
    for j in xrange(0,len(VarFunArr)):
        
#        print 'Variogram Fitting ' + optFunNam[j] + '\n'
        
        VarProb = Optimization('Variogram Fitting: ' + optFunNam[j], OptFun)
        VarProb.addObj('RMSE')
        VarProb.addVar('Sill','c',lower=Sb[0],upper=Sb[1],value=sr)
        VarProb.addVar('Range','c',lower=Rb[0],upper=Rb[1],value=rr)
        VarProb.addVar('Nugget','c',lower=Nb[0],upper=Nb[1],value=nr)
        VarProb.addVar('Exponent (a)','c',lower=ab[0],upper=ab[1],value=ar)
        VarProb.addVar('Rank (v)','c',lower=vb[0],upper=vb[1],value=vr)
        
        args = (SVExp, j, VarFunArr, Var, Res, Mdl)
        optmz = ALHSO()
        optmz.setOption('fileout',0)
        optmz(VarProb)
    
#        print VarProb.solution(0)
#        print 'Variogram adjusted at ' + time.ctime()
    
    # Get pointer to best semivariogram
    k = numpy.argmin(Res)
    xopt = Var[k]
    ModOpt = Mdl[k]
    
    # print 'Theoretical variogram fit - Done! at ' + time.ctime()
    return xopt, ModOpt, VarFunArr

def Kriging_core(ModOpt,sPOI,Loc,VarFunArr,xopt,CovMea,Prec,KrigType):
    '''
    Kriging core where interpolating algorithms is taking place.
    
    Parameters
    ----------
        **ModOpt** -- Pointer to optimal semivariogram model \n
        **sPOI** -- Single point of interest to interpolate (targets)\
            ``[t,2]`` \n
        **Loc** -- Gauge location for interpolation ``[x,2]`` \n
        **VarFunArr** -- Array with pointer to functions in VariogramFit module
            \n
        **xopt** -- vector with optimal semivariogram parameters ``[5]`` \n
        **CovMea** -- Gauge records covariance matrix ``[x,x]`` \n
        **Prec** -- Precipitation register for gauges ``[n,x]`` \n
        **KrigType** -- Type of Kriging to be used. 'Sim' for Simple and 'Ord'\
            for Ordinary Kriging
    
    Returns
    -------
        **Z** -- Interpolation for each target and time step ``[n,1]`` \n
        **SP** -- Interpolation variance field ``[1]``\n
    '''
    POID = Dist.target(Loc,[sPOI])[0]
    SVm = []
    for j in xrange(0,len(Loc)):
        SVm.append(VarFunArr[ModOpt](POID[j],xopt))

    if KrigType is 'Ord':
    #Ordinary Kriging
        CovMea = numpy.row_stack((CovMea,numpy.ones(len(CovMea))))
        CovMea = numpy.column_stack((CovMea,numpy.ones(len(CovMea))))
        CovMea[-1,-1] = 0.0
        SVm.append(1.0)
        
        SVr = numpy.array(CovMea)
        if linalg.det(CovMea) == 0:
            print('Non-singular covriance matrix - Sorry, cannot invert \n')
            return 9999*numpy.ones(len(Prec)), 9999*numpy.ones(len(Prec))
        InvSVr = linalg.inv(SVr)

        Z = []        
        WM= numpy.dot(InvSVr,SVm)
        for i in xrange(0,len(Prec)):
            Ztemp = numpy.dot(WM[:-1],Prec[i])
            Ztemp = numpy.clip(Ztemp,0,max(Prec[i])) # cutoff at 0 and max prec
            Z.append(Ztemp)        
#        S = xopt[0]*numpy.ones(len(SVm)-1)-SVm[:-1]
        S = SVm[:-1]
        SP = xopt[0] - (numpy.dot(WM[:-1],numpy.transpose(S))) - WM[-1]
        
    elif KrigType is 'Sim':
    #Simple Kriging
        SVr = numpy.array(CovMea)
        if linalg.det(CovMea) == 0:
            print('Non-singular covriance matrix - Sorry, cannot invert \n')
            return 9999*numpy.ones(len(Prec)), 9999*numpy.ones(len(Prec))
        InvSVr = linalg.inv(SVr)
        
        Z = []        
        WM= numpy.dot(InvSVr,numpy.transpose(numpy.array(SVm)))    
        for i in xrange(0,len(Prec)):
            Ztemp = numpy.dot(WM,Prec[i])
            Ztemp = numpy.clip(Ztemp,0,max(Prec[i])) # cutoff at 0 and max prec
            Z.append(Ztemp)        
        S = xopt[0]*numpy.ones(len(SVm))-SVm
        SP = xopt[0] - (numpy.dot(WM,numpy.transpose(S)))
    else:
        print 'I pity the fool for no chosing Kriging type'
        print 'only available Ord and Sim \n'
        Z = 9999*numpy.ones(len(Prec))
        SP = 9999*numpy.ones(len(Prec))
    
    return Z,SP

def Krig(MaxDist, POI, Loc, Prec, CovMea, ModOpt, xopt,
         VarFunArr, tmin=0, tmax='def', MinNumSt=3, KrigType='Sim'):
    '''
    Parsing of data for Kriging interpolation. This function contains\
    execution of interpolation as well. \n
    Parameters
    ----------
        **MaxDist** -- Initial search radius for nearby stations \n
        **POI** -- Points of interest to interpolate (targets) ``[t,2]`` \n
        **Loc** -- Gauge location for interpolation ``[x,2]`` \n
        **Prec** -- Precipitation register for gauges ``[n,x]`` \n
        **CovMea** -- Gauge records covariance matrix ``[x,x]`` \n
        **ModOpt** -- pointer to optimal semivariogram model \n
        **xopt** -- vector with optimal semivariogram parameters ``[5]`` \n
        **VarFunArr** -- Array with pointer to functions in VariogramFit module
            \n
        **tmin** -- Initial time step to be interpolated \n
        **tmax** -- Final time step to be interpolated \n
        **MinNumSt** -- Minimum number of stations within the search radius \n
        **KrigType** -- Type of Kriging to be used. 'Sim' for Simple and 'Ord'\
            for Ordinary Kriging
    
    Returns
    -------
        **Z** -- Interpolation for each target and time step ``[n,t]`` \n
        **SP** -- Interpolation variance field ``[t,1]``\n
        **ZAvg** -- Average of interpolated field ``[n,1]``
    '''
    if tmax == 'def':
        tmax = len(Prec)
    tmin = int(tmin)
    tmax = int(tmax)
    PrecSec = Prec[tmin:tmax]
    
    # Reduce measurements to relevant locations for the targets
    Z = []
    SP = []
    for kk in xrange(0,len(POI)):
        # Check if there are enough stations for interpolation, otherwise,
        # increase search radius
        POIDt = Dist.target(Loc, [POI[kk]])[0]
        TNS = 0
        MaxDist2 = MaxDist
        while TNS <= MinNumSt:
            #Minimum distance
            PlaceLoc = [i for i,v in enumerate(POIDt) if v > MaxDist2]
            TNS = len(Loc)-len(PlaceLoc)
            MaxDist2 += 1.0
        # Reduction of relevant stations (reduced data and cov matrices)            
        RedLoc = numpy.delete(Loc,PlaceLoc,0)
        RedPrec = numpy.delete(PrecSec,PlaceLoc,1)    
        RedCov = CovMea[:]    
        RedCov = numpy.delete(RedCov,PlaceLoc,0)
        RedCov = numpy.delete(RedCov,PlaceLoc,1)
        
        # Kriging interpolation
        TempRes = Kriging_core(ModOpt,POI[kk],RedLoc,VarFunArr, xopt,RedCov,RedPrec,KrigType)
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
    
    Loc, POI, Prec = DataLoad.lcsv(SiteInfo, XYTargets, DataRecord)
    SVExp, CovMea = exp_semivariogram(Prec, Loc)
    xopt, ModOpt, VarFunArr = theor_variogram(SVExp)
    Z, SP, ZAvg = Krig(xopt[0]/3, POI, Loc, Prec, CovMea, ModOpt, xopt,
         VarFunArr, tmin=0, tmax='def', MinNumSt=3, KrigType='Sim')
    

    DataSave.spkl(Z, 'VarField')
    DataSave.spkl(ZAvg, 'AvgVar')
    DataSave.spkl(SP, 'VarUnc')
    return Z, SP, ZAvg

def test():
    '''
    Module testing function
    '''
    Loc, POIC, Prec = DataLoad.lcsv('TestData\GaugeLoc.csv',
                                       'TestData\InterpPts.csv',
                                       'TestData\Dataset.csv')
    Loc = numpy.array(Loc)/1000.0
    POIC = numpy.array(POIC)/1000.0
    SVExp, CovMea = exp_semivariogram(Prec, Loc)
    xopt, ModOpt, VarFunArr = theor_variogram(SVExp)
    Z, SP, ZAvg = Krig(10.0, POIC, Loc, Prec, CovMea, ModOpt, xopt, 
                               VarFunArr,10 ,20, 'Ord')
    print 'Ordinary Kriging working fine'
    
    Z, SP, ZAvg = Krig(10.0, POIC, Loc, Prec, CovMea, ModOpt, xopt, 
                               VarFunArr,10 ,20, 'Sim')
    print 'Simple Kriging working fine'
    return 'Ended normally, module should be working properly'