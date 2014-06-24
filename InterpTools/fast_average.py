# -*- coding: utf-8 -*-
"""
==========================
Fast Average Precipitation
==========================
Implemented by Juan Chacon @ UNESCO-IHE
Integrated Water Systems and Governance Department
Hydroinformatics Laboratory

Average precipitation from multiple stations

* Pre requisites
    you will need the following libraries, not coming alongside with the\ 
    Anaconda ditribution (recommended)

* Functions
    * Run: Calculate average precipitation from several stations, as shown in\
    Lindström et al [1997]

* Use policy
    * you should include the respective citation to the authors
    * if you find this tool usefull, you will give the main author a beer next\
    time you see him :)
    
* References
    * Lindström et al., “Development and Test of the Distributed HBV-96\
    Hydrological Model.”
"""
import numpy

def Run(Data, Covariance='def'):
    '''
    Calculates average precipitation from several stations. Data consists in\
    a matrix element containing all the precipitation from the stations. If\
    covariance matrix can be provided it can be passed as an optional argument.
    
    Paramters
    ----------
        **Data** -- Matrix of recordings from stations. Data has to be\
        oriented column wise
        
        **Covariance** -- Matrix of covariance between stations. If it is not\
        provided, Covariance data will be calculated out of the given dataset.
        
    Results
    -------
        **AData** -- Average precipitation based on punctual readings, for all\
        the precipitation events in the serie.
    '''
    if Covariance is 'def':    
        Covariance = numpy.cov(numpy.transpose(Data))
    
    Data = numpy.array(Data)
    if numpy.linalg.det(Covariance) == 0:
        print 'Singular covariance matrix... cannot make it'
        return 9999*numpy.ones(len(Data))
    PGuess = numpy.average(Data,1) #initial guess of average precipitation 
    WSt = []
    for i in xrange(0,len(Covariance)):
        WSt.append(numpy.cov(Data[:,i],PGuess)[0][1])
    WSI = numpy.dot(numpy.linalg.inv(Covariance),WSt)
    
    AData = []
    for i in xrange(0,len(Data)):
        AData.append(PGuess[i] + numpy.dot(WSI,Data[i,:]-PGuess[i]))
    AData = numpy.clip(AData,0,max(AData))
    
    return AData

