#!/usr/bin/env python
# -*- coding: utf-8 -*-
r'''
Tools
=====

smooth :
    TODO
windows_doy :
    TODO
aggregate :
    TODO
group_by :
    TODO
movmean : Moving average
    TODO
movmeannan : Moving average ignoring NaN values
    TODO
rolling_window : Generate a ndarray with the rolling window
    TODO

'''
import numpy
import scipy
import scipy.stats

def main():
    pass

def smooth(x, window_len=11, window='hanning', circular_reflect=True):
    r'''smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    Parameters
    ----------
        x : 
            the input signal 
        window_len : 
            the dimension of the smoothing window; should be an odd integer
        window : 
            the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 
            'blackman' flat window will produce a moving average smoothing.

    Returns
    -------
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    '''

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    ## Fix this to include not a reflection but a circular reflect
    if circular_reflect:
        s = numpy.r_[x[(x.shape[0]-window_len+1):], x, x[:(window_len-1)]]
    else:
        s = numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w = numpy.ones(window_len,'d')
    else:
        w = eval('numpy.'+window+'(window_len)')

    y = numpy.convolve(w/w.sum(),s,mode='valid')
    return y

def group_by(data, group_array, axis=0, accumulate=False):
    r'''
    Function to group a time series by either month or days or years.. etc or
    a combination any of these. If a combination is used it must pe provided
    already as the group array. For example to get all the months in all the 
    years you could use this grouping array
        group_array = (year * 100.0) + (month * 1.0)
        
    or for daily values
        group_array = (year * 1000.0) + (doy * 1.0), but in this case
        aggregate provides much faster results!
    
    Parameters
    ----------
    data : numpy array
        Climatic data 
    group_array : numpy array
        Date value of climatic data, must be same size, and can be month, year
        etc...
    accumulate : Boolean
        After grooping this will performa a NaNsum on the grouped data. This
        can be used to get for example yearly totals of rainfall 
    
    Returns
    -------
    data_grouped :
        asd
        
    Example
    -------
    
        >>> tool.group_by(rainfall, month)
    '''
    
    if data.ndim == 1:
        group_no_duplicates = numpy.unique(group_array)
        #group_no_duplicates = list(set(group_array))
    else:
        ## Decide what to do if data is not a 1D array
        pass
        group_no_duplicates = numpy.unique(group_array)

    if accumulate:
        data_grouped = numpy.zeros(group_no_duplicates.shape[0])
        index = 0
    else:
        data_grouped = []
    
    for group in group_no_duplicates:
        if accumulate:
            data_grouped[index] = numpy.nansum(data[group_array == group])
            #data_grouped[index] = np.nansum(data[group_array == group], axis)
            index = index + 1            
        else:
            temp = data[(group_array == (group))]
            temp_clean = temp[~numpy.isnan(temp)]
            data_grouped.append(temp_clean)
        
    return data_grouped
    
def aggregate(variable, ts, ts_new, accumulate=False, threshold=None):
    r'''Aggregate a variable based on a new time step

    Parameter
    ---------
    variable :
        The climatic variable to process. An array with CONTINUOUS time steps.
    ts :
        Time step of the array  hours, days, months, etc...
    ts_new :
        New time step to aggregate the data, is relative to the old time step
        and must be a multiple of it: 
            if ts = 1 hour, ts_new = 2 hours or 3 or 4 or n
            if ts = 6 hour, ts_new = 12 hours or 24 etc...
    accumulate :
        True if when aggregating the sum is wanted. Otherwise the mean is used
        on the aggregated values
    threshold :
        Defines the minimum ratio of non nanvalues in an aggregated group
        in order to accept or reject the group. If 0.5 at least half of the
        values in the aggregation must be non NaN. By default it uses the 
        condition where at least ONE value in the aggregation is non NaN
    
    Returns
    -------
    A new array aggregated to the specified new time step.
    '''
    if accumulate:
        # check in other numpy scipy versions
        # fun_acc = np.nansum
        fun_acc = scipy.nansum
    else:
        #fun_acc = np.nanmean    
        fun_acc = scipy.stats.stats.nanmean       
    
    n_ = variable.shape[0]
    fr = ts_new / (ts * 1.0) 
    m = numpy.floor(n_ / fr)
    var_reshape = variable[:m*fr].reshape(m, fr)
    
    if threshold == None:
        threshold_new = (fr - 1) / fr
    else:
        threshold_new = threshold
        
    index = numpy.sum(numpy.isnan(var_reshape)*1, 1) <= numpy.round(fr*threshold_new)

    if numpy.floor(fr) - fr <> 0:
        raise RuntimeError("""The desired time step is not a multiple of 
        the recorded time step""")

    X = fun_acc(var_reshape[index], 1)
    
    return X

def window_doy(doy, day_window, group_tail=True):
    r'''Generates a grouping 1d array that is taking into account a window
    in days over all the years in a time series. This grouping array can then
    be used with the group_by function. The idea was to generate different
    grouping that not necesarilly follow the traditional monthly convention.
    
    Parameters
    ----------
    doy : 1d numpy array
        as
    window : integer
        days
    group_tail : Bool
        group last group?
        
    Returns
    -------
    1d array
    '''
    n_group = numpy.ceil(366.0/day_window)
    last_group = n_group - 1
    doy_new = numpy.floor((doy-1)/day_window)
    
    if group_tail:
        doy_new[doy_new==last_group] = last_group-1

    return doy_new

def sample_properties(data, ts, ts_new, lag, accumulate=True, DEBUG=False,
                      MY_METHOD=True, THRESHOLD=None, dry_threshold=0.0):
    r'''Function to calculate statistical properties of rainfall series at 
    different aggregation times.
  
    Parameters
    ----------
    data : 
        Rainfall series as a numpy array
    ts : 
        time step of registration of the time series. Expected 1 hour
    ts_new : 
        desired time of step output of the time series [hours]
    lag : 
        lag step for autocorrelation calculation, defined based on 
        aggregation step.
    acc : 
        If working for example with rainfall, setting this to true will
        accumulate fo rthe ts_new instead of calculating the mean.
    DEBUG :
        
    MY_METHOD :
        
    THRESHOLD :
        
  
    Returns
    -------
    list :
        
    
    Eh : 
        Expected value (mean)
    VARh : 
        Sample variance (over n-1)
    CVh : 
        Coefficient of variation (sqrt(VAR)/mean)
    Rlh : 
        Autocorrelation coefficient for a given lag period
    SKh : 
        Skewness
    FFh : 
        Frequency of dry periods
    Fddh : 
        Frequency of dry dry periods
    Fwwh : 
        Frequency of wet wet periods
    '''  
    dof = 1
    ## Accumulate function: Sum or average
    if accumulate:
        fun_acc = numpy.nansum
    else:
        #fun_acc = numpy.nanmean
        fun_acc = scipy.stats.stats.nanmean
    
    if MY_METHOD:
        ## MY METHOD
        ## First group the time series to thenew time step and check ifnan
        ## There seems to be a bug with the np.all function. so a wat to do it is:    
        ## multiply the isnan vector by 1 and sum rowwise
        ## if a row has only nan the value of the sum in the row should be = to fr
        ## so the indexes to take should be those where the value is < fr
        ## after cleaning the series from grouped rows that only have nans
        ## a nansum or nanmean can be applied (depending if the varuiable)
        ## will be accumulated or not. THRESHOLD is a new variable that if defined
        ## would indicate the fraction of NaN accepted in aggregation... so 0.5
        ## implies that at least half of the aggregated elements need to be valid
        n_ = data.shape[0]
        fr = ts_new / (ts * 1.0) 
        m = numpy.floor(n_ / fr)
        var_reshape = data[:m*fr].reshape(m, fr)
        if THRESHOLD == None:
            #to generalize the THRESHOLD with none should evaluate to the 
            THRESHOLD_new = (fr-1) / fr            
            #var_index = np.sum(np.isnan(var_reshape) * 1, 1) < fr    
        else:
            THRESHOLD_new = THRESHOLD
        
        var_index = (numpy.sum(numpy.isnan(var_reshape) * 1, 1) <= numpy.round(fr * THRESHOLD_new))

        if numpy.floor(fr) - fr <> 0:
            raise RuntimeError("""The desired time step is not a multiple of 
            the recorded time step""")
    
        X = fun_acc(var_reshape[var_index], 1)
    else:
        ## Original Method... I think is wrong
        fr = ts_new / (ts * 1.0) 
        var_clean = data[~numpy.isnan(data).any(1)]
        n = var_clean.shape[0]
        m = numpy.floor(n / fr)
        var_reshape = var_clean[:m*fr].reshape(m, fr)
        X = fun_acc(var_reshape, 1)
        
        
    ## Main statistics  
    Eh = numpy.mean(X)           
    VARh = numpy.var(X, ddof=dof)    
    CVh = numpy.sqrt(VARh)/Eh

    temp_1 = X - Eh
    temp_2 = numpy.correlate(temp_1,temp_1,'valid')
    Rlh = (numpy.correlate(temp_1[lag:],temp_1[:-lag],'valid')/temp_2)[0]
    #Rlh = Rl[0]
    
    SKh = scipy.stats.skew(X)                           
    
    ## If accumulate is sum then calculate additional statistics
    #if acc:
        ## Additiona statistics only applicable to rainfall
    Xp = (X > dry_threshold) * 1.0;
    dXp = scipy.diff(Xp, n=1, axis=0)
    
    wd = numpy.sum(dXp==-1.0)
    dw = numpy.sum(dXp==1.0)
    Xpw = (dXp + 1.0) * (Xp[:-1])  
  
    ww = numpy.sum(Xpw)
    dd = m - wd - dw - ww

    FFh = 1.0 - (numpy.sum(Xp) / m)  
    Fwwh = ww / (wd + ww)
    Fddh = dd / (dd + dw)

    statistics = [Eh, VARh, CVh, Rlh, SKh, FFh, Fddh, Fwwh]
    
    if DEBUG:
        print "New time step: %5i" % ts_new
        print "--------------------"
        print "Eh   : %5.10f" % Eh
        print "VARh : %5.10f" % VARh
        print "CVh  : %5.10f" % CVh
        print "Rlh  : %5.10f" % Rlh
        print "SKh  : %5.10f" % SKh
        print "FFh  : %5.10f" % FFh
        print "Fddh : %5.10f" % Fddh
        print "Fwwh : %5.10f" % Fwwh
        print "\n"            
#    else:       
#        statistics = [Eh, VARh, CVh, Rlh, SKh]
#        if DEBUG:
#            print "Eh   : %f5.10" % Eh
#            print "VARh : %f5.10" % VARh
#            print "CVh  : %f5.10" % CVh
#            print "Rlh  : %f5.10" % Rlh
#            print "SKh  : %f5.10" % SKh   
#            print "\n"
    
    return statistics

def sample_correlation(rain, ts, aggregation_step, total_time):
    r'''Function to calculate statistical properties of rainfall series at 
    different aggregation times.
  
    Parameters
    ----------
    rain : 
        Rainfall series as a numpy array
    ts : 
        time of registration of the time series. Normally 1 hour
    aggregation_step : 
        normally a multiple of 1 hour [hours]
    lag : 
        lag step for autocorrelation calculation, defined based on 
        aggregation step. 
  
    Returns
    -------
    Rlh : Autocorrelation coefficient for a given lag period
    '''      
    ## Define common use functions
    mean = numpy.mean  

    fr = aggregation_step/ts  
    rain_clean = rain[~numpy.isnan(rain).any(1)]

    if int(fr) - fr <> 0:
        print('The aggregation period is not a multiple of the delta period')

    n = rain_clean.shape[0]
    m = numpy.int(n/fr)
    
    ## Remove NAN values
    rainfall_reshape = rain_clean[:m*fr].reshape(m,fr)
    ## Aggregate values
    X = numpy.sum(rainfall_reshape,1)
    
    n_ = X.shape[0]
  
    Eh = mean(X)
    
    steps = numpy.int(total_time/aggregation_step)
           
    temp_1 = X - Eh
    ## To make this go faster change the full, by a for that creates n loops
    ## for n steps, and using valid, this way will be much faster!!!!!!    
    
    auto_correlation = numpy.correlate(temp_1,temp_1,'valid') 
    Rl = (numpy.correlate(temp_1,temp_1,'full')/auto_correlation)[n_:n_+steps]
    Rlh = Rl
    
    X = [i*aggregation_step for i in range(1,steps+1)]
  
    return [X, Rlh]

def movmean(x, window_size):
    r'''
    Calculate the moving average for a given window size

    Obs:
    ----
    Faster and lower memory consumption but cannot ignore NaN values in the
    window  
    '''
    w = numpy.ones(window_size, 'd')
    y = numpy.convolve(w/w.sum(), x.flat, mode='valid')
    return y

def movnanmean(x, window_size, axis=0, accumulate=False):
    r'''
    Calculate the moving average for a given window size, without NaN

    Parameters
    ----------
    x : 1d numpy array
        .    
    window_size : int
        .
    axis : int
        .
        
    Returns
    -------
    
    Obs:
    ----
    Slower and higher memory consumption but can ignore NaN values in the
    window    
    '''
    if x.ndim == 1:
        y = rolling_window(x, window_size)    

    if accumulate:
        z = scipy.nansum(y, axis)
    else:
        z = scipy.stats.stats.nanmean(y, axis)
    return z
    

def rolling_window(a, window_size):
    '''
    Make an ndarray with a rolling window of the last dimension
  
    Parameters
    ----------
    a : array_like
       Array to add rolling window to
    window_size : int
        Size of rolling window

    Returns
    -------
    Array that is a view of the original array with a added dimension
    of size w.

    Examples
    --------
    >>> x=np.arange(10).reshape((2,5))
    >>> rolling_window(x, 3)
    array([[[0, 1, 2], [1, 2, 3], [2, 3, 4]],
           [[5, 6, 7], [6, 7, 8], [7, 8, 9]]])

    Calculate rolling mean of last dimension:
    >>> np.mean(rolling_window(x, 3), -1)
    array([[ 1.,  2.,  3.],
           [ 6.,  7.,  8.]])
         
    Developed 
    ---------
    Name : 
        Erik Rigtorp 
    E-mail : 
        erik@rigtorp.com
    URL : 
        http://mail.scipy.org/pipermail/numpy-discussion/2011-January/054401.html
    '''
    if window_size < 1:
        raise ValueError, "`window` must be at least 1."
        
    if window_size > a.shape[-1]:
        raise ValueError, "`window` is too long."

    shape = a.shape[:-1] + (a.shape[-1] - window_size + 1, window_size)
    strides = a.strides + (a.strides[-1],)

    return numpy.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def bootstrap(sample, resamples, verbose=False):
    r'''
    Function to perform a simple nonparametric bootstrap standard error
    
    Parameters
*    ----------
    sample : 1D Array or List
        TODO
    resamples : 
        number of resamples to draw from the sample
    verbose : False
        print results
        
    Returns
    -------
    
    '''
    Et = []
    VARt = []
    CVt = []
    STDt = []
    SKt = []
    
    N = resamples
    ddof = 1

    sample_clean = sample[(~numpy.isnan(sample))]
    n = sample_clean.shape[0]
    
    for i in range(N):
        index = numpy.random.randint(0, n, n)
        #print len(set(index))
        resample = sample_clean[index]
        E = numpy.mean(resample)
        VAR = numpy.var(resample, ddof=ddof)
        STD = numpy.sqrt(VAR)
        CV = STD/E
        SK = scipy.stats.skew(resample)
        
        Et.append(E)
        VARt.append(VAR)
        STDt.append(STD)
        CVt.append(CV)
        SKt.append(SK)
    
    SE_E = numpy.sqrt(numpy.var(Et, ddof=ddof))
    SE_VAR = numpy.sqrt(numpy.var(VARt, ddof=ddof))
    SE_STD = numpy.sqrt(numpy.var(STDt, ddof=ddof))
    SE_CV = numpy.sqrt(numpy.var(CVt, ddof=ddof))
    SE_SK = numpy.sqrt(numpy.var(SKt, ddof=ddof))
    
    if verbose:
        print('\tStandard errors\n\t---------------')
        print('\t\tSE-mean: ' + str(SE_E))
        print('\t\tSE-variance: ' + str(SE_VAR))
        print('\t\tSE-Coefficient of variation: ' + str(SE_CV))
        print('\t\tSE-Standard deviation: ' + str(SE_STD))
        print('\t\tSE-Skewness: ' + str(SE_SK))
        print('')

    # h = plt.hist(Et, bins=25, normed=True)    

    dic = {'E': SE_E, 
           'VAR': SE_VAR,
           'ST': SE_STD,
           'CV': SE_CV,
           'SK': SE_SK}
           
    return dic

    #return [SE_E, SE_VAR, SE_STD, SE_CV, SE_SK]

def mean_confidence_interval(data, confidence=0.95):
    r'''TODO
    
    Parameters
    ----------
    data : 
        TODO
    confidence :
        TODO
        
    Returns
    -------
    TODO
    '''
    a = 1.0 * numpy.array(data)
    n = a.shape[0]
    m = numpy.mean(a)
    se = numpy.std(data, ddof=1)/numpy.sqrt(n)
    #se = sp.stats.stderr(a)
    h = se * scipy.stats.t._ppf((1+confidence)/2., n-1)
    #return {'mean': m, 'lci': m - h, 'lci': m + h}
    return [m, m - h, m + h]

def series_statistics(dic, ts, aggregation, lag=1, group_type='m', accumulate=True):
    r"""Calculates statistics per period based on aggregation period. Usefull 
    for analyzing the trends of the statistics at different aggregation periods.
    
    Parameters
    ----------
    dic : Pythhon dictionary
        TODO
    odir : String
        TODO
    ts : Integer
        time series step [hours]
    aggregation : Python list or 1D numpy array of integers 
        TODO
    lag : Integer
        TODO
    group_type :
        TODO
    
    Returns
    -------
    None

    """ 
    data = dic['rainfall']

    # Depending on the grouping 
    if group_type == 'm':
        group_array = dic['month']
        grouped = group_by(data, group_array)
        group_title = ['January', 'February', 'March', 'April', 'May', 'June', 
                       'July', 'August', 'September', 'October', 'November', 
                       'December']
    elif group_type > 0:
        doy = group_array = dic['doy']
        group_tail = True
        group_array = window_doy(doy, group_type, group_tail=group_tail)
        grouped = group_by(data, group_array)
        interval = range(0, 365, group_type)        
        interval[-1] = 365
        group_title = [str(interval[i]+1) + '-' + str(interval[i+1]) for i in range(len(interval)-1)]

    n_group = len(grouped)

    # Empty periodly containers
    E_m, VAR_m, CV_m, Rl_m, SK_m, DSF_m = [], [], [], [], [], []

    # Empty periodly containers
    Ey, VARy, CVy, Rly, SKy, DSFy = [], [], [], [], [], []
         
    ## Find the yearly statistics     
    for ap in aggregation:           
        (Eh, VARh, CVh, Rlh, SKh, FFh, Fddh, Fwwh) = sample_properties(
        data, ts, ap, lag, accumulate=accumulate)

        Ey.append(Eh)
        VARy.append(VARh)
        CVy.append(CVh)
        Rly.append(Rlh)
        SKy.append(SKh)
        DSFy.append(FFh)

    del data

    ## numpy arrays holding the yearly statistics for different aggregations
    Ey = numpy.array(Ey)
    VARy = numpy.array(VARy)
    CVy = numpy.array(CVy)
    Rly = numpy.array(Rly)
    SKy = numpy.array(SKy)
    DSFy = numpy.array(DSFy)

    ## Find the statistics for the periods
    for group in range(n_group):
        E, VAR, CV, Rl, SK, DSF = [], [], [], [], [], []

        for ap in aggregation:           
            (Eh, VARh, CVh, Rlh, SKh, FFh, Fddh, Fwwh) = sample_properties(
            grouped[group], ts, ap, lag, accumulate=accumulate)

            E.append(Eh)
            VAR.append(VARh)
            CV.append(CVh)
            Rl.append(Rlh)
            SK.append(SKh)
            DSF.append(FFh)

        E = numpy.array(E)
        VAR = numpy.array(VAR)
        CV = numpy.array(CV)
        Rl = numpy.array(Rl)
        SK = numpy.array(SK)
        DSF = numpy.array(DSF)

        E_m.append(E)
        VAR_m.append(VAR)
        CV_m.append(CV)
        Rl_m.append(Rl)
        SK_m.append(SK)
        DSF_m.append(DSF)
                   
    ## Find max to adjust the y axis useful when ploting later on
    E_max, VAR_max, CV_max, Rl_max, SK_max= [], [], [], [], []

    for group in range(n_group):
        E_max.append(numpy.max(E_m[group]))
        VAR_max.append(numpy.max(VAR_m[group]))
        CV_max.append(numpy.max(CV_m[group]))
        Rl_max.append(numpy.max(Rl_m[group]))
        SK_max.append(numpy.max(SK_m[group]))

    E_max = numpy.max(E_max)
    VAR_max = numpy.max(VAR_max) 
    CV_max = numpy.max(CV_max)
    Rl_max = numpy.max(Rl_max)
    SK_max = numpy.max(SK_max)
    DSF_max = 1.0

    E_min = 0.0
    VAR_min = 0.0
    CV_min = 0.0
    Rl_min = -Rl_max
    SK_min = 0.0
    DSF_min = 0.0 

    new_dic = {}

    new_dic['labels'] =  ['Mean [mm]\n', 
                          r'Variance [${mm}^{2}$]', 
                          'Coefficient of variation [-]\n', 
                          'Skewness [-]\n', 
                          r'Lag ' + str(lag) + ' Autocovariance [-]', 
                          r'Dry Spell Fraction [-]']
    new_dic['group_title'] = group_title
    new_dic['group_type'] = 'm'
    new_dic['n'] = n_group
    new_dic['ts'] = ts
    new_dic['aggregation'] = aggregation
    new_dic['lag'] = lag
    new_dic['year'] = [Ey, VARy, CVy, SKy, Rly, DSFy]
    new_dic['period'] = [E_m, VAR_m, CV_m, SK_m, Rl_m, DSF_m]
    new_dic['period_max'] = [E_max, VAR_max, CV_max, SK_max, Rl_max, DSF_max]
    new_dic['period_min'] = [E_min, VAR_min, CV_min, SK_min, Rl_min, DSF_min]
    
    return new_dic

def get_max(rain, year, duration, accumulate=False):
    r'''Function to calculate the extreme value distribution of rainfall 
    for different rainfall duration and return periods. This can be then used
    to produce intensity-duration-frequency curves.
    
    Parameters
    ----------
    rain : 
        numpy array
    year : 
        numpy array
    duration : 
        list
    return_period :
                
    Returns
    -------
    A dictionary
    
    '''
    movmean = movnanmean
    year_max = int(numpy.amax(year))
    year_min = int(numpy.amin(year)) 
      
    rain_max_list = []
    #year = numpy.transpose(year)   
    
    for step in duration:
        rain_year_max = []   
        for y in range(year_min, (year_max+1)):           
            rain_per_year = rain[(year==y)]

            ##This is to comply with the rolling window function in case there
            ## are years with no data.
            if rain_per_year.shape[0] > step: 
                rain_year_m_ave = movmean(rain_per_year, step, 1, accumulate)    
                temp_calc = numpy.nanmax(rain_year_m_ave)
            
            ## If the moving average produces NaN when calculating the
            ## nanmax ignores the NaN. 
            ## Is this appropriate in every case?
                if not numpy.isnan(temp_calc): 
                    rain_year_max.append(temp_calc)
        
        size = len(rain_year_max)
        sorting =  numpy.sort(rain_year_max)
        rain_max_list.append(sorting)

    return {'rain': rain_max_list,
            'size': size}

def fix_name(n, digit='0', digits=3):
    r'''
    '''
    if n == 0:
        x = 1
    else:
        x = numpy.int(numpy.floor(numpy.log10(n))) + 1

    list_name = [digit]*(digits-x)
    list_name.append(str(n))
    name = ''.join(list_name)

    return name

if __name__ == '__main__':
    main()