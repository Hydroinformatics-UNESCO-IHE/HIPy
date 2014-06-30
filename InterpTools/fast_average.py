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


def Run(data, covariance=None):
    '''
    Calculates spatial average from several stations. Data consists in\
    a matrix element containing all the recordas from the stations. If\
    covariance matrix can be provided it can be passed as an optional argument.

    Paramters
    ----------
        **data** -- Matrix of recordings from stations. data has to be\
        oriented column wise

        **covariance** -- Matrix of covariance between stations. If it is not\
        provided, covariance data will be calculated out of the given dataset.

    Results
    -------
        **spatial_average** -- Average precipitation based on punctual \
        readings, for all\
        the precipitation events in the serie.
    '''
    if covariance is None:
        covariance = numpy.cov(numpy.transpose(data))

    data = numpy.array(data)
    if numpy.linalg.det(covariance) == 0:
        print 'Singular covariance matrix... cannot make it'
        return 9999*numpy.ones(len(data))

    pre_average_guess = numpy.average(data, 1)
    station_weights = []
    for i in xrange(len(covariance)):
        station_weights.append(numpy.cov(data[:, i], pre_average_guess)[0][1])
    station_weights_updated = numpy.dot(numpy.linalg.inv(covariance),
                                        station_weights)

    spatial_average = []
    for i in xrange(len(data)):
        spatial_average.append(pre_average_guess[i] +
                               numpy.dot(station_weights_updated, data[i, :] -
                               pre_average_guess[i]))

    spatial_average = numpy.clip(spatial_average, 0, max(spatial_average))

    return spatial_average
