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

import csv
import numpy
from numpy import linalg
import cPickle
import random
from pyOpt import ALHSO, Optimization
#------------------------------------------------------------------------------

def load_data(SiteInfo, XYTargets, DataRecord):
    '''
    Load data from CSV files as presented in the description into working files
    inside the module. No data processing taking place here.    
    
    Input files description
    ---------------------------------------------------------------------------
     This script uses CSV (comma separated) files for import and ouput of data
     The necessary files are:
       -> SiteInfo.csv - File containing coordinates of sensors in the same
           order as the other files (This is mandatory).
          Format:
          [NAME, X, Y]
           
       -> XYTargets.csv - File containing the location of sample points inside
           the catchment that are going to be used to calculate the average of
           the given variable.
          Format:
           [X, Y, Catchment Number]
           
       -> DataRecord.csv - File containing the registers of the variable that
           is going to be interpolated. For historical reasons will be noted as
           Prec (precipitation).
          Format:
           [St1Data, St2Data, ..., StnData]
    '''
    Loc = []
    with open(SiteInfo,'rb') as SiteInfo:
        Lines = csv.reader(SiteInfo)
        Lines.next()
        for row in Lines:
            Loc.append([float(row[1]),float(row[2])])
    print 'Gauge Location, Imported'
    print '' 
    
    POI = []
    CN = []
    with open(XYTargets,'rb') as POIf:
        Lines = csv.reader(POIf)
        Lines.next()
        for row in Lines:
            POI.append([float(row[0]),float(row[1])])
            CN.append(int(row[2]))
    print 'Sampling Points, Imported'
    print ''    
            
    Prec = []
    with open(DataRecord,'rb') as Data:
        Lines = csv.reader(Data)
        Lines.next()
        for row in Lines:
            Prec.append([float(x) for x in row[:]])
    print 'Data, Imported'
    print ''
    return Loc, POI, CN, Prec

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
                 VariogramFit.SVPower, VariogramFit.SVMatern]
    
    # Names of functions for display only
    optFunNam = ['Exponential','Gaussian','Spherical','Cubic',
                 'Pentaspherical','Sinehole','Power','Matern']
    
    # Initial seed for variogram fit
    sr = random.uniform(Sb[0],Sb[1])
    rr = random.uniform(Rb[0],Rb[1])
    nr = random.uniform(Nb[0],Nb[1])
    ar = random.uniform(ab[0],ab[1])
    vr = random.uniform(vb[0],vb[1])
    return sr, rr, nr, ar, vr
    
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
        optmz = ALHSO()
        optmz(VarProb)
    
        print VarProb.solution(0)
        print ''    
    
    # Get position of best semivariogram
    k = numpy.argmin(Res)
    xopt = Var[k]
    ModOpt = Mdl[k]
    del Var
    del Res
    del Mdl
    
    print 'Theoretical variogram fit - Done!'
    print ''
    return xopt, ModOpt, VarFunArr

def Kriging_core(ModOpt,POI,Loc,VarFunArr,xopt,k,CovMea,MeaRow):
    '''
    Kriging core where interpolating algorithms is taking place.    
    '''
    #Check for singularities in covariance matrix
    SVr = numpy.array(CovMea)
    if linalg.det(CovMea) == 0:
        print('Non-singular covriance matrix - Sorry, cannot invert')
        return(9999,9999)
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
    return Z,SP

def Krig(MaxDist, MinNumSt=3, POI, Loc, Prec, CovMea, ModOpt, xopt,
         VarFunArr, ID, CN):
    '''
    Kriging implementation in neighboorhood.
    MaxDist - Maximum radius for search of interpolation stations
    MinNumSt - Minimum number of stations to perform the interpolation
    '''
    # Select only the targets within the catchment
    if ID == '':
        POIC = POI
    else:
        POIC = []
        for j in xrange(0,len(POI)):
            if ID == CN[j]:
                POIC.append(POI[j])
    
    # Reduce measurements to relevant locations for the targets
    POIDred = numpy.zeros([len(POIC),len(Loc)])
    for k in xrange(0,len(POIC)):
        for j in xrange(0,len(Loc)):
            # Calculate distance from target to stations
            POIDred[k][j] = numpy.sqrt((POIC[k][0]-Loc[j][0])**2 +
            (POIC[k][1]-Loc[j][1])**2 )
    
    # Check if there are enough stations for interpolation, otherwise, increase 
    # the search radius
    POIDt = numpy.transpose(POIDred)
    RedLoc = [] #initialise to jump for the first time in the cycle
    while len(RedLoc) <= MinNumSt:
       
        #Minimum distance reduction
        RedLoc = [] #Reduced Location of sensors
        RedPrec = [] #Reduced data from sensors
        PlaceLoc = [] #Station number list to be leftout
        for k in xrange(len(Loc)):
            if min(POIDt[k]) > MaxDist:
                PlaceLoc.append(k)
        RedLoc = numpy.delete(Loc,PlaceLoc,0)
        MaxDist = MaxDist + 1.0
    
    # Trimming of Precipitation data and Covariance matrix (reduced matrices)
    RedPrec = numpy.delete(Prec,PlaceLoc,1)    
    RedCov = CovMea[:]    
    RedCov = numpy.delete(RedCov,PlaceLoc,0)
    RedCov = numpy.delete(RedCov,PlaceLoc,1)
        
    # Kriging interpolation
    Z = []
    ZAvg = []
    SP = []
    f = numpy.zeros(len(RedPrec[0])) #Line of zeros if no data exists
    for ii in xrange(17136,17496):
        if max(RedPrec[ii]) == 0:
            Z.append(f)
            SP.append(f)
            ZAvg.append(0)
        else:        
            TempRes = Kriging_core(ModOpt,POIC,RedLoc,VarFunArr,
                                              xopt,k,RedCov,RedPrec[ii])
            Z.append(TempRes[0])
            SP.append(TempRes[1])
            ZAvg.append(numpy.average(TempRes[0]))
        if ii%100 == 0:
            print 'Interpolated precipitation register '+ str(ii)
    return Z, SP, ZAvg

def data_save_pkl(var, ID):
    '''
    Function used to save data from interpolation results into native python
    format (PKL)
    '''
    with open(ID +'.pkl','w') as DataFile:
        cPickle.dump(var,DataFile)
    print ID + ' Pickled'
    print ''
    return
    
def data_save_csv(var, ID):
    '''
    Function used to save data from interpolation results into comma separated
    value (CSV)
    ''' 
    with open(ID +'.csv', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                                quotechar=',', 
                                quoting=csv.QUOTE_MINIMAL)        
        for j in xrange(len(var)):
            spamwriter.writerow([var[j]])    
    print ID + ' CSV-ed'
    print ''    
    return