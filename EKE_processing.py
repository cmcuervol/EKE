# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import scipy.stats as sts
from netCDF4 import Dataset

# function to read files in month order!!, taken from internet
import re
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts       = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
    
#function to find the position of a value in an array    
def find_nearest(array,value1):
    idx1 = (np.abs(array-value1)).argmin()
    return idx1  

dir  = '/scratch/uni/u234/u301103/' 
dir2 = '/home/zmaw/u301103/codes/'

#dir = '/media/esmontoyadu/ESTEFANIA/dropbox/estefania/ICCS/Schaffen/eke'
#dir2 = '/media/esmontoyadu/ESTEFANIA/dropbox/estefania/ICCS/Schaffen/'
#================================================================================================
# Read files
#================================================================================================   


def read_files(lat_area,lon_area,lev_area,mean_axis1,concat= False,mean_axis2=False):
    ''' As EKE is done with u and v fields
    as they are originally needed (i.e. 
    full field or mean over latitudes) in order
    to have correct magnitudes, EKE is computed
    every time depending on the analysis neded
    
    lat_area: list, min and max latitude to work
    lon_area: list, min and max longitude to work
    lev_area: list, min and max altitude to work
    mean_axis1: integer
    mean_axis2: Only used if a time series is needed
    concat: if the area needed is from east to west a concatanation is needed
    '''
    
    # general data
    # 1.1 geo location
    C   = Dataset(dir+'/wind/ERA_vwind_25_daily_19790101_19790131.nc')
    lat = C.variables['latitude'][:]
    lon = C.variables['longitude'][:]
    lev = C.variables['level'][:]

    C.close()

    date_era = pd.date_range('1979-01-01','2018-05-31 18:00',freq='6H')
    

    # 1.2 define area to read

    lat1,lat2 = find_nearest(lat,lat_area[0]),find_nearest(lat,lat_area[1])
    lon1,lon2 = find_nearest(lon,lon_area[0]),find_nearest(lon,lon_area[1])
    lev1,lev2 = find_nearest(lev,lev_area[0]),find_nearest(lev,lev_area[1])
    
    # cut lat, lon
    if concat == False:
        lat,lon,lev = lat[lat1:lat2+1],lon[lon1:lon2+1],lev[lev1:lev2+1]
    else:
        lon = np.concatenate([lon[:lon1+1],lon[lon2:]])
        lat,lev = lat[lat1:lat2+1],lev[lev1:lev2+1]

    if mean_axis1 == 1 and mean_axis2==False: # mean over latitudes
        axis = lon
    elif mean_axis1==2 and mean_axis2==False:
        axis = lat # mean over longitudes
        
    if mean_axis2<>False:
        axis = [1]

    # 1.3 variables to fill in 
    vwind_cycle  = np.empty([len(date_era),len(axis),len(lev[lev1:lev2+1])])
    uwind_cycle  = np.empty([len(date_era),len(axis),len(lev[lev1:lev2+1])])

    #~ # 1.4 filtering dates of strong and weak years
    #~ if strength <> None:
        #~ strength_dates = [i for i in date_era if i.year in strength]
        
    # 1.5 read files
    for ind,infile in enumerate(lev[lev1:lev2+1]):
        
        print infile,ind
        
        # read band-pass filtered data level by level
        if mean_axis2==False:
            if concat==False:
                A     = Dataset(dir+'/fulltime/u_f_%s.nc'%infile)
                u_aux = A.variables['var1'][:,lat1:lat2+1,lon1:lon2+1].mean(axis=mean_axis1)
                
                A.close()
                
                A     = Dataset(dir+'/fulltime/v_f_%s.nc'%infile)
                v_aux = A.variables['var2'][:,lat1:lat2+1,lon1:lon2+1].mean(axis=mean_axis1)
                
                A.close()
            else:
                A     = Dataset(dir+'/fulltime/u_f_%s.nc'%infile)
                u_aux = (np.concatenate([A.variables['var1'][:,lat1:lat2+1,:lon1+1],A.variables['var1'][:,lat1:lat2+1,lon2:]],axis=2)).mean(axis=mean_axis1)
                
                A.close()
                
                A     = Dataset(dir+'/fulltime/v_f_%s.nc'%infile)
                v_aux = (np.concatenate([A.variables['var2'][:,lat1:lat2+1,:lon1+1],A.variables['var2'][:,lat1:lat2+1,lon2:]],axis=2)).mean(axis=mean_axis1)
                
                A.close()
        else:
            if concat==False:
                # It is a time series
                A     = Dataset(dir+'/fulltime/u_f_%s.nc'%infile)
                u_aux = A.variables['var1'][:,lat1:lat2+1,lon1:lon2+1].mean(axis=mean_axis1).mean(axis=mean_axis2)
                
                A.close()
                
                A     = Dataset(dir+'/fulltime/v_f_%s.nc'%infile)
                v_aux = A.variables['var2'][:,lat1:lat2+1,lon1:lon2+1].mean(axis=mean_axis1).mean(axis=mean_axis2)
                
                A.close()
            else:
                A     = Dataset(dir+'/fulltime/u_f_%s.nc'%infile)
                u_aux = (np.concatenate([A.variables['var1'][:,lat1:lat2+1,:lon1+1],A.variables['var1'][:,lat1:lat2+1,lon2:]],axis=2)).mean(axis=mean_axis1).mean(axis=mean_axis2)
                
                A.close()
                
                A     = Dataset(dir+'/fulltime/v_f_%s.nc'%infile)
                v_aux = (np.concatenate([A.variables['var2'][:,lat1:lat2+1,:lon1+1],A.variables['var2'][:,lat1:lat2+1,lon2:]],axis=2)).mean(axis=mean_axis1).mean(axis=mean_axis2)
                
                A.close()
            
        vwind_cycle[:,:,ind] = v_aux
        uwind_cycle[:,:,ind] = u_aux
        
    return [uwind_cycle,vwind_cycle,lat,lon,lev,date_era]



def read_files2(lat_area,lon_area,lev_area,mean_axis1,concat= False,mean_axis2=False):
    ''' As EKE is done with u and v fields
    as they are originally needed (i.e. 
    full field or mean over latitudes) in order
    to have correct magnitudes, EKE is computed
    every time depending on the analysis neded
    
    lat_area: list, min and max latitude to work
    lon_area: list, min and max longitude to work
    lev_area: list, min and max altitude to work
    mean_axis1: integer
    mean_axis2: Only used if a time series is needed
    '''
    
    # general data
    # 1.1 geo location
    C   = Dataset(dir+'/wind/ERA_vwind_25_daily_19790101_19790131.nc')
    lat = C.variables['latitude'][:]
    lon = C.variables['longitude'][:]
    lev = C.variables['level'][:]

    C.close()

    date_era = pd.date_range('1979-01-01','2018-05-31 18:00',freq='6H')
    

    # 1.2 define area to read

    lat1,lat2 = find_nearest(lat,lat_area[0]),find_nearest(lat,lat_area[1])
    lon1,lon2 = find_nearest(lon,lon_area[0]),find_nearest(lon,lon_area[1])
    lev1,lev2 = find_nearest(lev,lev_area[0]),find_nearest(lev,lev_area[1])
    
    # cut lat, lon
    if concat == False:
        lat,lon,lev = lat[lat1:lat2+1],lon[lon1:lon2+1],lev[lev1:lev2+1]
    else:
        lon = np.concatenate([lon[:lon1+1],lon[lon2:]])
        lat,lev = lat[lat1:lat2+1],lev[lev1:lev2+1]

    if mean_axis1 == 1 and mean_axis2==False: # mean over latitudes
        axis = lon
    elif mean_axis1==2 and mean_axis2==False:
        axis = lat # mean over longitudes
        
    if mean_axis2<>False:
        axis = [1]

    # 1.3 variables to fill in 
    vwind_cycle  = np.empty([len(date_era),len(axis),len(lev[lev1:lev2+1])])
    uwind_cycle  = np.empty([len(date_era),len(axis),len(lev[lev1:lev2+1])])

    #~ # 1.4 filtering dates of strong and weak years
    #~ if strength <> None:
        #~ strength_dates = [i for i in date_era if i.year in strength]
        
    # 1.5 read files
    for ind,infile in enumerate(lev[lev1:lev2+1]):
        
        print infile
        
        # read band-pass filtered data level by level
        if mean_axis2==False:
            if concat==False:
                A     = Dataset(dir+'/fulltime/EKE_%s_fulltime.nc'%infile)
                u_aux = A.variables['u'][:,lat1:lat2+1,lon1:lon2+1].mean(axis=mean_axis1)
                v_aux = A.variables['v'][:,lat1:lat2+1,lon1:lon2+1].mean(axis=mean_axis1)
                
                A.close()
            else:
                A     = Dataset(dir+'/fulltime/EKE_%s_fulltime.nc'%infile)
                u_aux = (np.concatenate([A.variables['u'][:,lat1:lat2+1,:lon1+1],A.variables['u'][:,lat1:lat2+1,lon2:]],axis=2)).mean(axis=mean_axis1)
                v_aux = (np.concatenate([A.variables['v'][:,lat1:lat2+1,:lon1+1],A.variables['v'][:,lat1:lat2+1,lon2:]],axis=2)).mean(axis=mean_axis1)
                
                A.close()
        else:
            if concat==False:
                # It is a time series
                A     = Dataset(dir+'/fulltime/EKE_%s_fulltime.nc'%infile)
                u_aux = A.variables['u'][:,lat1:lat2+1,lon1:lon2+1].mean(axis=mean_axis1).mean(axis=mean_axis2)
                v_aux = A.variables['v'][:,lat1:lat2+1,lon1:lon2+1].mean(axis=mean_axis1).mean(axis=mean_axis2)
                
                A.close()
            else:
                A     = Dataset(dir+'/fulltime/EKE_%s_fulltime.nc'%infile)
                u_aux = (np.concatenate([A.variables['u'][:,lat1:lat2+1,:lon1+1],A.variables['u'][:,lat1:lat2+1,lon2:]],axis=2)).mean(axis=mean_axis1).mean(axis=mean_axis2)
                v_aux = (np.concatenate([A.variables['v'][:,lat1:lat2+1,:lon1+1],A.variables['v'][:,lat1:lat2+1,lon2:]],axis=2)).mean(axis=mean_axis1).mean(axis=mean_axis2)
                
                A.close()

        vwind_cycle[:,:,ind] = v_aux
        uwind_cycle[:,:,ind] = u_aux
        
    return [uwind_cycle,vwind_cycle,lat,lon,lev,date_era]


#================================================================================================
## all pre-processing functions are here to be called in other codes
#================================================================================================   

   
def fourier2D(lat,lon,y):
    ''' band pass filter
    execute fourier and leave only 3-10 days frequencies, 
    then compute the inverse to have the time series'''
    
    d0,d1,d2 = np.shape(y)
    
    inverse = np.empty([d0,d1,d2])  # empty matix to fill filtered data

    # as is an array i have to filter point by point
    for ind,la in enumerate(lat):
        for jnd,lo in enumerate(lon):
            N    = len(y[:,ind,jnd])
            Yf   = np.fft.fft(y[:,ind,jnd])
            freq = np.fft.fftfreq(N)
            
            # data is every 6hours and I need to only 3-10 days frequencies
            upper_limit = 10*24/6. # units every 6 hours
            lower_limit = 3*24/6.  # units evey 6 hours
            
            Yf[np.abs((1/freq))>upper_limit] = 0.0
            Yf[np.abs((1/freq))<lower_limit] = 0.0
            
            # Inverse..now all values outside 3-10 freq are deleted
            inv = np.real(np.fft.ifft(Yf))
            inverse[:,ind,jnd] = inv
            
    return inverse
    
from scipy.signal import butter, lfilter
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs,order=5):
    d0,d1,d2 = np.shape(data)
    
    y_filter = np.empty([d0,d1,d2])  # empty matix to fill filtered data

    # as is an array i have to filter point by point
    for ind in range(d1):
        for jnd in range(d2):
            b, a                = butter_bandpass(lowcut, highcut, fs, order=order)
            y_filter[:,ind,jnd] = lfilter(b, a, data[:,ind,jnd])
    return y_filter
    
    
def eke_space_full_time_2D(lat,lon,u,v):
    '''once band pass filter is done for each wind 
    component eke is computed in the array
    EKE dimensions t,lat,lon'''
    
    EKE = np.empty([len(u),len(lat)-1,len(lon)-1])
    

    for ind in range(1,len(lat)):
        for jnd in range(1,len(lon)):
            #note that: ind,jnd are array positions starting in one
            
            # monthly mean conditions in each point to make anomalies
            # zonal component
            umean2_i  = ((np.mean(u[:,ind,jnd],axis=0))**2.)
            umean2_i1 = ((np.mean(u[:,ind-1,jnd],axis=0))**2.) 
            # meridional component
            vmean2_j  = ((np.mean(v[:,ind,jnd],axis=0))**2.)
            vmean2_j1 = ((np.mean(v[:,ind,jnd-1],axis=0))**2.)
            
            # variables needed-> for the anomaly (all time record)
            # zonal component
            u2_i      = u[:,ind,jnd]**2. 
            u2_i1     = u[:,ind-1,jnd]**2
            # meridional component
            v2_j      = v[:,ind,jnd]**2.
            v2_j1     = v[:,ind,jnd-1]**2
            
            EKE[:,ind-1,jnd-1] = ((((u2_i - umean2_i) + (u2_i1 - umean2_i1))/2.) + (((v2_j - vmean2_j) + (v2_j1 - vmean2_j1))/2.))/2.         

    return EKE


def t_test(seriefull,seriehalf):
    ''' execute a t-test
    of two sets'''
    
    t,p    = sts.ttest_ind(seriefull, seriehalf)
    
    if p>0.05:
        result = -999
    if p<0.05:
        result = p
    return result
    

def linear_regression(xx,yy):
    '''perform a linear regression
    # y = b0 + b1x, x 
    and gives the R^2 and standard 
    deviation of residuals ro
    using pandas series'''
        
    b1_y = yy.apply(lambda x: x - yy.mean())
    b1_x = xx.apply(lambda x: x - xx.mean())
    b1   = ((b1_x*b1_y).sum())/ ((b1_x**2.).sum())
    
    b0   = yy.mean() - (b1*xx.mean())
    
    # squared error
    
    ss_x = ((b1_x)**2.).sum() 
    ss_y = ((b1_y)**2.).sum()  
    ss_l = ((yy - (b0 + b1*xx))**2.).sum()
    
    r2   = 1 - (ss_l/ss_y)
    
    # standard deviation of residuals or root mean squared error
    ro   = np.sqrt(ss_l/(len(yy)-2))
    
    # linear fit
    line = b0 + b1*xx
    
    return b0,b1,line,r2,ro
