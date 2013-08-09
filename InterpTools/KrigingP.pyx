# -*- coding: utf-8 -*-
"""
Kriging Interpolation for average areal precipitation
Implemented by Juan Chacon @ UNESCO-IHE
Integrated Water Systems and Governance Department
Hydroinformatics Core

-Use policy-
1- you should include the respective citation to the authors
2- if you find this tool usefull, you will give the main author a beer next
    time you see him :)

- Pre-requisites
you will need the following libraries, not coming alongside with the Anaconda
ditribution (it is recommended that you use this one, but any will do)
    
    -> pyx import - some functions are compiled in Cython to improve 
        performance, so you will need this to use them.
    -> csv - Read and write csv files    
    -> numpy - for obvious reasons
    -> cPickle - To store data in native python binary files.
    -> pyOpt - Optimisation engine, used to solve the semivariogram fitting  

- Suggestions
Make sure that you check which files are you going to keep, precipitation and
uncertainty maps can be quite large, so make sure you understand your needs and
the code as well before using these tools. File pickling is at the end of the 
script.
"""

#Libraries to be imported
#------------------------------------------------------------------------------
import pyximport
pyximport.install()
import VariogramFit
import Dist

import numpy
from numpy import linalg
import random
from pyOpt import ALHSO, Optimization
import DataSave
import DataLoad

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
except:
    raise ImportError('mpi4py is required for parallelization')

#------------------------------------------------------------------------------

def exp_semivariogram(Prec, Loc):
    '''
    Evaluation of experimental semivariogram    
    '''
    ## Removal of no precipitation events   
    WetMeas = []
    for i in xrange(0,len(Prec)):
        if numpy.max(Prec[i]) > 3:
            WetMeas.append(Prec[i])    
    
    ## Measurement covariance
    CovMea = numpy.cov(numpy.transpose(WetMeas))
    if linalg.det(CovMea) == 0:
        print 'Warning, Covariance Matrix is non-signular. Check your data'
    else:
        print 'Determinant Covariance Matrix', str(linalg.det(CovMea))
    print ''

    Dis = Dist.between(Loc)
    print 'Distance Matrix - Done'
    print ''
    
    ## Experimental Semivariogram
    SVExp = []
    for i in xrange(0,len(CovMea)-1):
        for j in xrange(i+1,len(CovMea)):
            Cov = CovMea[i][j]
            Lag = Dis[i][j]
            SVExp.append([Lag,Cov])
            
    print 'Experimental semivariogram - Done'
    print ''
    return SVExp, CovMea

def theor_variogram(SVExp, Sb=(0.01,400), Rb=(2,20), Nb=(0,400),
                  ab=(0,2), vb=(0,1000)):
    '''
    Fitting of theoretical variogram to experimental variogram. sample values
    for boundaries of variables are provided.
    '''                      
    
    # Array with functions to be called from the Variograms library
    VarFunArr = [VariogramFit.SVExponential, VariogramFit.SVGaussian, 
                 VariogramFit.SVSpherical, VariogramFit.SVCubic,
                 VariogramFit.SVPentaspherical, VariogramFit.SVSinehole, 
                 VariogramFit.SVPower]
    
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
        
    print 'Initialising Variogram fit'
    print ''
    
    # Optimisation starts to minimise differences between experimental and 
    # theoretical semivariograms
    for j in xrange(0,len(VarFunArr)):
        
        print 'Variogram Fitting ' + optFunNam[j]
        print ''
        
        VarProb = Optimization('Variogram Fitting: ' + optFunNam[j], OptFun)
        VarProb.addObj('RMSE')
        VarProb.addVar('Sill','c',lower=Sb[0],upper=Sb[1],value=sr)
        VarProb.addVar('Range','c',lower=Rb[0],upper=Rb[1],value=rr)
        VarProb.addVar('Nugget','c',lower=Nb[0],upper=Nb[1],value=nr)
        VarProb.addVar('Exponent (a)','c',lower=ab[0],upper=ab[1],value=ar)
        VarProb.addVar('Rank (v)','c',lower=vb[0],upper=vb[1],value=vr)
        
        args = (SVExp, j, VarFunArr, Var, Res, Mdl)
        optmz = ALHSO(pll_type='POA')
        optmz.setOption('fileout',0)
        optmz(VarProb)
    
        print VarProb.solution(0)
        print ''    
    
    # Get position of best semivariogram
    k = numpy.argmin(Res)
    xopt = Var[k]
    ModOpt = Mdl[k]
#    del Var
#    del Res
#    del Mdl
    
    print 'Theoretical variogram fit - Done!'
    print ''
    print 'range is equal to: ' + str(xopt[1])
    return xopt, ModOpt, VarFunArr

def Kriging_core(ModOpt,sPOI,Loc,VarFunArr,xopt,CovMea,Prec,KrigType):
    '''
    Kriging core where interpolating algorithms is taking place.    
    '''
    #Check for singularities in covariance matrix
    
    POID = Dist.target(Loc,[sPOI])[0]
#    SVm = numpy.zeros([len(POI),len(Loc)])    
    SVm = []
    for j in xrange(0,len(Loc)):
        # Calculate variance with between observation and targets
        SVm.append(VarFunArr[ModOpt](POID[j],xopt))

    if KrigType is 'Ord':
        #Ordinary Kriging
        CovMea = numpy.row_stack((CovMea,numpy.ones(len(CovMea))))
        CovMea = numpy.column_stack((CovMea,numpy.ones(len(CovMea))))
        CovMea[-1,-1] = 0.0
        SVm.append(1.0)
        
        SVr = numpy.array(CovMea)
        if linalg.det(CovMea) == 0:
            print('Non-singular covriance matrix - Sorry, cannot invert')
            return(9999,9999)
        InvSVr = linalg.inv(SVr)

        Z = []        
        WM= numpy.dot(InvSVr,SVm)
#        print WM
#        print Prec[0]
        for i in xrange(0,len(Prec)):
#            ResPrec = Prec[i]-numpy.average(Prec[i])
            Ztemp = numpy.dot(WM[:-1],Prec[i])
            Ztemp = numpy.clip(Ztemp,0,max(Prec[i])) ## cutoff at 0 and max prec
            Z.append(Ztemp)        
        S = xopt[0]*numpy.ones(len(SVm))-SVm
        SP = (numpy.dot(WM,numpy.transpose(S)))
        
    elif KrigType is 'Sim':
    #Simple Kriging
        SVr = numpy.array(CovMea)
        if linalg.det(CovMea) == 0:
            print('Non-singular covriance matrix - Sorry, cannot invert')
            return(9999,9999)
        InvSVr = linalg.inv(SVr)
        
        Z = []        
        WM= numpy.dot(InvSVr,numpy.transpose(numpy.array(SVm)))    
        for i in xrange(0,len(Prec)):
            Ztemp = numpy.dot(WM,Prec[i])
            Ztemp = numpy.clip(Ztemp,0,max(Prec[i])) ## cutoff at 0 and max prec
            Z.append(Ztemp)        
        S = xopt[0]*numpy.ones(len(SVm))-SVm
        SP = (numpy.dot(WM,numpy.transpose(S)))
    else:
        print 'I pity the fool for no chosing Kriging type'
        print 'only available Ord and Sim'
        Z = 9999*numpy.ones(len(Prec))
        SP = 9999*numpy.ones(len(Prec))
    
    return Z,SP

def Krig(MaxDist, POIC, Loc, Prec, CovMea, ModOpt, xopt,
         VarFunArr, tmin=0, tmax='def', MinNumSt=3, KrigType='Sim'):
    '''
    Kriging implementation in neighboorhood.
    MaxDist - Maximum radius for search of interpolation stations
    MinNumSt - Minimum number of stations to perform the interpolation
    '''
    if tmax == 'def':
        tmax = len(Prec)
    
    tmin = int(tmin)
    tmax = int(tmax)
    PrecSec = Prec[tmin:tmax]
    # Reduce measurements to relevant locations for the targets
    Z = []
    SP = []
    for kk in xrange(0,len(POIC)):

        # Check if there are enough stations for interpolation, otherwise, increase 
        # the search radius
        POIDt = Dist.target(Loc, [POIC[kk]])[0]
        
        TNS = 0
        MaxDist2 = MaxDist
        while TNS <= MinNumSt:
            #Minimum distance reduction
            PlaceLoc = [i for i,v in enumerate(POIDt) if v > MaxDist2]
            TNS = len(Loc)-len(PlaceLoc)
            MaxDist2 += 1.0
        print TNS
        # Trimming of Precipitation data and Covariance matrix (reduced matrices)            
        RedLoc = numpy.delete(Loc,PlaceLoc,0)
        RedPrec = numpy.delete(PrecSec,PlaceLoc,1)    
        RedCov = CovMea[:]    
        RedCov = numpy.delete(RedCov,PlaceLoc,0)
        RedCov = numpy.delete(RedCov,PlaceLoc,1)
        
        # Kriging interpolation
        TempRes = Kriging_core(ModOpt,POIC[kk],RedLoc,VarFunArr,
                                          xopt,RedCov,RedPrec,KrigType)
        if Z == []:
            Z = numpy.vstack(TempRes[0])
        else:
            temp = numpy.vstack(TempRes[0])
            Z = numpy.hstack((Z,temp))
        SP.append(TempRes[1])
#        print 'Interpolated precipitation in target '+ str(kk)

    ZAvg = numpy.average(Z,1)
    return Z, SP, ZAvg


def Simple_Krig(SiteInfo, XYTargets, DataRecord):
    '''
    Wrapper for Kriging interpolation of all data, and save in PKL format    
    '''
    Loc, POIC, Prec = DataLoad.lcsv(SiteInfo, XYTargets, DataRecord)
    SVExp, CovMea = exp_semivariogram(Prec, Loc)
    xopt, ModOpt, VarFunArr = theor_variogram(SVExp)
    Z, SP, ZAvg = Krig(10, POIC, Loc, Prec, CovMea, ModOpt, xopt, VarFunArr)

    DataSave.spkl(Z, 'PrecField')
    DataSave.spkl(ZAvg, 'AvgPrec')
    DataSave.spkl(SP, 'PrecUnc')
    
    