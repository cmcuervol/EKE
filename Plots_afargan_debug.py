# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import calendar

from netcdftime import utime
import datetime as dt

from EKE_processing import *
from plot_functions import *

import glob


dir  = '/scratch/uni/u234/u301103/' 
dir2 = '/home/zmaw/u301103/codes/'

#dir = '/media/esmontoyadu/ESTEFANIA/dropbox/estefania/ICCS/Schaffen/eke'
#dir2 = '/media/esmontoyadu/ESTEFANIA/dropbox/estefania/ICCS/Schaffen/'

# 1.1 definition of weak and strong years: reference Afargan 2017
weak_yearsP   = [1982, 1989, 1990, 2006, 2009]
strong_yearsP = [1981, 1983, 1984, 1995, 2003]

weak_yearsA   = [1983, 1995, 1998, 2000, 2012]
strong_yearsA = [1989, 2010, 2014]


#2.1 Select two areas Atlantic and pacific

Area = {'Area1': {'title':'20 - 70 {}^oN',
                  'eke_weak':weak_yearsA,
                  'eke_strong':strong_yearsA,
                  'Data':read_files([70,20],[-180,180],[100,1000],1,concat= False,mean_axis2=False),
                  'u': read_files2([70,20],[-180,180],[100,1000],1,concat= False,mean_axis2=False)},
        'Area2': {'title':'40 - 60 {}^oN',
                  'eke_weak':weak_yearsA,
                  'eke_strong':strong_yearsA,                  
                  'Data':read_files([60,40],[-180,180],[100,1000],1,concat= False,mean_axis2=False),
                  'u': read_files2([70,20],[-180,180],[100,1000],1,concat= False,mean_axis2=False)}}
                 
lon0            = find_nearest(Area['Area1']['Data'][3],0)
lev250          = find_nearest(Area['Area1']['Data'][4],250)
      

# 2.2 vertically integrated EKE zonal average

for i in Area.keys():

    # first compute the integration over the vertical component (dp) of EKE
    eke              = np.empty([12,len(Area[i]['Data'][3]),len(Area[i]['Data'][-2])])
    # compute eke
    u_aux    = Area[i]['Data'][0]
    v_aux    = Area[i]['Data'][1]
    date_era = Area[i]['Data'][-1]
    u_jet    = Area[i]['u'][0]

    for ind,la in enumerate(Area[i]['Data'][-2]):
        print la,ind
        # integration in levels for every point
        # pandas allows to reduce time operations but with 2D arrays
        u_aux2   = pd.DataFrame(index=date_era,data=u_aux[:,:,ind])
        v_aux2   = pd.DataFrame(index=date_era,data=v_aux[:,:,ind]) 
        
        # annual cycle
        u_aux3 = (u_aux2**2.).groupby([u_aux2.index.month]).mean()
        v_aux3 = (v_aux2**2.).groupby([v_aux2.index.month]).mean()   
        
        eke[:,:,ind]  = ((u_aux3) + (v_aux3))/2. # Joules/m2

    EKE_integrated = (np.trapz(eke,x=Area[i]['Data'][-2]*100.,dx=2))/(9.8) # Joules
       
    # zonal everage wind (first was done over lon, here over lev and then annual cycle)
    to_plot_zonal = pd.DataFrame(index=date_era,data=u_jet.mean(axis=2))
    to_plot_zonal = (to_plot_zonal.groupby([to_plot_zonal.index.month]).mean()).values
    
    # the plot is from july to june
    to_plot2       = np.concatenate([EKE_integrated[6:,:],EKE_integrated[:6,:]],axis=0)
    # the plot is centered in 180 and not in zero!
    to_plot2 = np.concatenate([to_plot2[:,lon0:],to_plot2[:,:lon0]],axis=1)
    
    to_plot2_zonal = np.concatenate([to_plot_zonal[6:,:],to_plot_zonal[:6,:]],axis=0)
    to_plot2_zonal = np.concatenate([to_plot2_zonal[:,lon0:],to_plot2_zonal[:,:lon0]],axis=1)
        
    #fig,ax,cs = grafico_matrix(Area[i]['Data'][3],range(12),to_plot2.T/1000000.,\
    #            np.min(to_plot2)/1000000.,np.max(to_plot2)/1000000.,'jet',3,(12,7),norm=None) # MJ/m2

    fig,ax,cs = grafico_matrix(range(len(Area[i]['Data'][3])),range(12),to_plot2.T/1000000.,\
                0,0.3,'jet',3,(12,7),norm=None) 

    # add zonal wind lines
    #bounds = np.round( np.linspace(10,50,6, endpoint=True),1 )
    bounds = np.round( np.linspace(5,30,6, endpoint=True),1 )
    #z_cs   = plt.contour(range(12),Area[i]['Data'][3], to_plot2_zonal.T,levels=bounds,colors='k')
    z_cs   = plt.contour(range(12),range(len(Area[i]['Data'][3])), to_plot2_zonal.T,levels=bounds,colors='k')
    plt.clabel(z_cs, fontsize=9, fmt='%1.0f',inline=1)
    
    ax.set_title(r'EKE %s $(%s)$'%(i,Area[i]['title']),fontsize=24)

    cbar = fig.colorbar(cs,pad=0.1,fraction=0.046,orientation='vertical')
    cbar.ax.tick_params(labelsize=20)
    
    plt.ylabel('Longitude',fontsize=20)
    plt.xlabel('Month',fontsize=20)
    
    ax.set_xticks(range(12))      
    ax.set_xticklabels(['J','A','S','O','N','D','J','F','M','A','M','J'])
    
    ax.set_yticks(range(len(Area[i]['Data'][3]))[::35])
    ax.set_yticklabels(['0','87E','175E','95W','8W'])
    
    ax.tick_params(labelsize=20) 

    plt.savefig(dir2+'/Figures/'+i+'_vertically_integrated2.png',\
                dpi=300,transparent=True,bbox_inches='tight')
    plt.close('all')
    

# 2.3 latitude vs pressure

for i in Area.keys():

    # zonal mean (get rid of longitudes)    
    # compute eke
    u_aux    = Area[i]['Data'][0]
    v_aux    = Area[i]['Data'][1]
    date_era = Area[i]['Data'][-1]
    u_jet    = Area[i]['u'][0]
    
    to_plot = np.empty([12,len(Area[i]['Data'][3]),len(Area[i]['Data'][-2])])
    to_plot_zonal = np.empty([12,len(Area[i]['Data'][3]),len(Area[i]['Data'][-2])])
    for ind,la in enumerate(Area[i]['Data'][-2]):
        print ind,la
        # pandas allows to reduce time operations but with 2D arrays
        u_aux2   = pd.DataFrame(index=date_era,data=u_aux[:,:,ind])
        v_aux2   = pd.DataFrame(index=date_era,data=v_aux[:,:,ind])
        
        # annual cycle
        u_aux3 = (u_aux2**2.).groupby([u_aux2.index.month]).mean()
        v_aux3 = (v_aux2**2.).groupby([v_aux2.index.month]).mean()   
        
        eke  = ((u_aux3) + (v_aux3))/2. # Joules/m2

        to_plot[:,:,ind]       = eke.values
        
        u_jet_aux = pd.DataFrame(index=date_era,data=u_jet[:,:,ind])
        to_plot_zonal[:,:,ind] = (u_jet_aux.groupby([u_jet_aux.index.month]).mean()).values

    # month per month 
    for m in range(1,13):
        to_plot2       = to_plot[m-1,:,:]
        to_plot2 = np.concatenate([to_plot2[:,lon0:],to_plot2[:,:lon0]],axis=1)
        to_plot2_zonal = to_plot_zonal[m-1,:,:] 
        to_plot2_zonal = np.concatenate([to_plot2_zonal[:,lon0:],to_plot2_zonal[:,:lon0]],axis=1)   
    
        #~ fig,ax,cs = grafico_matrix(Area[i]['Data'][-2],Area[i]['Data'][2],to_plot2,\
                    #~ 0,140,'jet',1,(12,7),norm=None) # m2/s2

        fig,ax,cs = grafico_matrix(Area[i]['Data'][-2],range(len(Area[i]['Data'][3])),to_plot2.T,\
                    0,80,'jet',1,(12,7),norm=None) # m2/s2

        # add zonal wind
        bounds = np.round(np.linspace(5,30,20, endpoint=True),1 )
        #z_cs   = plt.contour(Area[i]['Data'][2],Area[i]['Data'][-2],to_plot2_zonal,levels=bounds,colors='k')
        z_cs   = plt.contour(range(len(Area[i]['Data'][3])),Area[i]['Data'][-2],to_plot2_zonal.T,levels=bounds,colors='k')
        # int labels
        # Recast levels to new class
        plt.clabel(z_cs,z_cs.levels[:10:2],fontsize=9, fmt='%1.0f',inline=1)

        ax.set_title(r'%s %s'%(i,calendar.month_abbr[m]),fontsize=24)
        
        ax.invert_yaxis()

        cbar = fig.colorbar(cs,pad=0.1,fraction=0.046,orientation='vertical')
        cbar.ax.tick_params(labelsize=20)
        
        plt.xlabel('Longitude',fontsize=20)
        plt.ylabel('Pressure [hPa]',fontsize=20)
        
        #ax.set_xticks(Area[i]['Data'][2][::5][::-1])
        
        ax.set_xticks(range(len(Area[i]['Data'][3]))[::35])
        ax.set_xticklabels(['0','87E','175E','95W','8W'])
                
        ax.tick_params(labelsize=20) 

        plt.savefig(dir2+'/Figures/'+i+'_latitude_vs_pressure_%s2.png'%m,\
                    dpi=300,transparent=True,bbox_inches='tight')
        plt.close('all')


# 2.5 strong - weak events

for i in Area.keys():
    
    for types in ['weak','strong']:

        # compute eke
        u_aux = Area[i]['Data'][0][:,:,lev250]
        v_aux = Area[i]['Data'][1][:,:,lev250]
        date_era = Area[i]['Data'][-1]
        strength_dates = [ii for ii in date_era if ii.year in Area[i]['eke_'+types]]

        # pandas allows to reduce time operations but with 2D arrays
        u_aux2   = pd.DataFrame(index=date_era,data=u_aux)
        v_aux2   = pd.DataFrame(index=date_era,data=v_aux)
        
        # take only weak or strong dates
        u_aux2 = u_aux2[u_aux2.index.isin(strength_dates)] 
        v_aux2 = v_aux2[v_aux2.index.isin(strength_dates)]
                
        # annual cycle
        u_aux3 = (u_aux2**2.).groupby([u_aux2.index.month]).mean()
        v_aux3 = (v_aux2**2.).groupby([v_aux2.index.month]).mean()   
        
        eke  = ((u_aux3) + (v_aux3))/2. # Joules/m2
        
        to_plot       = eke.values
        to_plot_zonal = (u_aux2.groupby([u_aux2.index.month]).mean()).values

        # the plot is from july to june
        to_plot2       = np.concatenate([to_plot[6:,:],to_plot[:6,:]],axis=0)
        to_plot2_zonal = np.concatenate([to_plot_zonal[6:,:],to_plot_zonal[:6,:]],axis=0)

        to_plot2 = np.concatenate([to_plot2[:,lon0:],to_plot2[:,:lon0]],axis=1)
        to_plot2_zonal = np.concatenate([to_plot2_zonal[:,lon0:],to_plot2_zonal[:,:lon0]],axis=1)

        #fig,ax,cs = grafico_matrix(Area[i]['Data'][3],range(12),to_plot2,\
        #            0,95,'jet',3,(12,7),norm=None) # MJ/m2

        fig,ax,cs = grafico_matrix(range(len(Area[i]['Data'][3])),range(12),to_plot2.T,\
                    0,95,'jet',3,(12,7),norm=None) # MJ/m2

        # add zonal wind
        bounds = np.round(np.linspace(10,80,12, endpoint=True),1 )
        #z_cs = plt.contour(range(12),Area[i]['Data'][2], to_plot2_zonal.T,levels=bounds,colors='k')
        z_cs = plt.contour(range(12),range(len(Area[i]['Data'][3])), to_plot2_zonal.T,levels=bounds,colors='k')
        plt.clabel(z_cs, z_cs.levels[:5:],fontsize=9, fmt='%1.0f',inline=1)
        

        ax.set_title(r'%s %s jet'%(i,types),fontsize=24)

        cbar = fig.colorbar(cs,pad=0.1,fraction=0.046,orientation='vertical')
        cbar.ax.tick_params(labelsize=20)

        plt.ylabel('Longitude',fontsize=20)
        plt.xlabel('Month',fontsize=20)

        ax.set_xticks(range(12))      
        ax.set_xticklabels(['J','A','S','O','N','D','J','F','M','A','M','J'])

        ax.set_yticks(range(len(Area[i]['Data'][3]))[::35])
        ax.set_yticklabels(['0','87E','175E','95W','8W'])

        ax.tick_params(labelsize=20) 

        plt.savefig(dir2+'/Figures/'+i+'_'+types+'2.png',\
                    dpi=300,transparent=True,bbox_inches='tight')
        plt.close('all')


# ask supervisor t-test!
#~ for i in Area.keys():
    
    #~ for types in ['weak','strong']:

        #~ to_plot        = Area[i]['Data']['eke_'+types][:,lev1,:,:].mean(axis=2) # Joules/m2
        #~ #to_plot_zonal = Area[i]['Data']['zonal_'+types][:,lev,:,:].mean(axis=2)
        #~ to_plot_pvalue = Area[i]['Data']['p_'+types]
        #~ to_plot_pvalue[to_plot_pvalue>0] = 1

        #~ # the plot is from july to june
        #~ to_plot2        = np.concatenate([to_plot[6:,:],to_plot[:6,:]],axis=0)
        #~ to_plot2_pvalue = np.concatenate([to_plot_pvalue[6:,:],to_plot_pvalue[:6,:]],axis=0)

        #~ fig,ax,cs = grafico_matrix(Area[i]['lat'],range(12),to_plot2.T,\
                    #~ 1.1,55,'jet',3,(12,7),norm=None) # MJ/m2

        #~ cs_p = plt.contour(range(12),Area[i]['lat'],to_plot_pvalue.T,1,colors='k')
        #~ cs_p = plt.contourf(range(12),Area[i]['lat'], to_plot_pvalue.T,1,colors='none',hatches=[None,'\\'])

        #~ ax.set_title(r'%s %s jet'%(i,types),fontsize=24)

        #~ cbar = fig.colorbar(cs,pad=0.1,fraction=0.046,orientation='vertical')
        #~ cbar.ax.tick_params(labelsize=20)

        #~ plt.ylabel('Latitude',fontsize=20)
        #~ plt.xlabel('Month',fontsize=20)

        #~ ax.set_xticks(range(12))      
        #~ ax.set_xticklabels(['J','A','S','O','N','D','J','F','M','A','M','J'])

        #~ ax.tick_params(labelsize=20) 

        #~ plt.savefig(dir2+'/Figures/'+i+'_'+types+'p_value.png',\
                    #~ dpi=300,transparent=True,bbox_inches='tight')
        #~ plt.close('all')

# 2.4 longitude temporal
# this one goes over northen hemisphere

Area_temporal = {'Area1': {'title':'20 - 70 {}^oN',\
                 'Data':read_files([70,20],[-180,180],[250,250],1,mean_axis2=False)},\
                 'Area2': {'title':'40 - 60 {}^oN',\
                 'Data':read_files([60,40],[-180,180],[250,250],1,mean_axis2=False)}}
                 
lon0            = find_nearest(Area_temporal['Area1']['Data'][3],0)

for i in Area_temporal.keys():
    
    # compute eke
    u_aux    = Area_temporal[i]['Data'][0][:,:,lev250]
    v_aux    = Area_temporal[i]['Data'][1][:,:,lev250]
    date_era = Area_temporal[i]['Data'][-1]

    # pandas allows to reduce time operations but with 2D arrays
    u_aux2   = pd.DataFrame(index=date_era,data=u_aux)
    v_aux2   = pd.DataFrame(index=date_era,data=v_aux) 
    
    # annual cycle
    u_aux3 = (u_aux2**2.).groupby([u_aux2.index.month]).mean()
    v_aux3 = (v_aux2**2.).groupby([v_aux2.index.month]).mean()   
    
    eke  = ((u_aux3) + (v_aux3))/2. # Joules/m2

    # the plot is from july to june
    to_plot = np.concatenate([eke.values[6:,:],eke.values[:6,:]],axis=0)
    
    # the plot is centered in 180 and not in zero!
    to_plot = np.concatenate([to_plot[:,lon0:],to_plot[:,:lon0]],axis=1)
    
    fig,ax,cs = grafico_matrix(range(12),range(len(Area_temporal[i]['Data'][3])),to_plot,\
                0,80,'jet',1,(12,7),norm=None) 

    ax.set_title(r'$(%s)$'%(Area_temporal[i]['title']),fontsize=24)

    cbar = fig.colorbar(cs,pad=0.1,fraction=0.046,orientation='vertical')
    cbar.ax.tick_params(labelsize=20)
    
    plt.xlabel('Longitude',fontsize=20)
    plt.ylabel('Month',fontsize=20)
    
    #ax.set_yticks([1,32,60,91,121,152,182,213,244,274,305,335])      
    ax.set_yticks(range(12))
    ax.set_yticklabels(['J','A','S','O','N','D','J','F','M','A','M','J'])
    
    ax.set_xticks(range(len(Area_temporal[i]['Data'][3]))[::35])
    ax.set_xticklabels(['0','87E','175E','95W','8W'])
    
    ax.tick_params(labelsize=20) 

    plt.savefig(dir2+'/Figures/'+i+'_longitude_time_evolution2.png',\
                dpi=300,transparent=True,bbox_inches='tight')
    plt.close('all')


#2.5 seasonal behaviour in the atlantic region EKE vs zonal wind
lon1 , lon2 = find_nearest(Area['Area2']['Data'][3],-70),find_nearest(Area['Area2']['Data'][3],-30)
u_jet = pd.DataFrame(index = Area['Area2']['u'][-1], data = Area['Area2']['u'][0][:,:,lev250]) 

# montly data to find max value
u_jet_month = u_jet.resample('M').mean()
# wind is at the latitude where maximum montlhy zonal wind takes place
timemax,lonmax = np.nonzero(u_jet_month.values>=u_jet_month.values.max())
u_jet = u_jet[lonmax]

# filtered data for seasonal behaviour
u_aux = Area['Area2']['Data'][0][:,lon1:lon2+1,lev250].mean(axis=1) # read_files does a mean over lat and here over lon
u_aux = pd.DataFrame(index = Area['Area2']['Data'][-1],data = u_aux)

# seasonal resample
eke_seasonal  = {}
wind_seasonal = {}
for y in range(1979,2019):
    eke_seasonal[str(y)] = {'JFM':float(u_aux[str(y)+'-01-01 00:00':str(y)+'-03-31 18:00'].mean()),\
                    'AMJ':float(u_aux[str(y)+'-04-01 00:00':str(y)+'-06-30 18:00'].mean()),\
                    'JAS':float(u_aux[str(y)+'-07-01 00:00':str(y)+'-09-30 18:00'].mean()),\
                     'OND':float(u_aux[str(y)+'-09-01 00:00':str(y)+'-12-31 18:00'].mean())}

    wind_seasonal[str(y)] = {'JFM':float(u_jet[str(y)+'-01-01 00:00':str(y)+'-03-31 18:00'].mean()),\
                    'AMJ':float(u_jet[str(y)+'-04-01 00:00':str(y)+'-06-30 18:00'].mean()),\
                    'JAS':float(u_jet[str(y)+'-07-01 00:00':str(y)+'-09-30 18:00'].mean()),\
                     'OND':float(u_jet[str(y)+'-10-01 00:00':str(y)+'-12-31 18:00'].mean())}
                     
wind_seasonal = pd.DataFrame(wind_seasonal).stack()
eke_seasonal  = pd.DataFrame(eke_seasonal).stack()

# linear regression
b0,b1,line,r2,ro = linear_regression(wind_seasonal['JFM'],eke_seasonal['JFM'])

plt.figure(figsize=(5,5))
plt.scatter(wind_seasonal['JFM'].values,eke_seasonal['JFM'].values,marker='o',c='royalblue',edgecolors='none',label='JFM')
plt.scatter(wind_seasonal['AMJ'].values,eke_seasonal['AMJ'].values,marker='o',c='seagreen',edgecolors='none',label='AMJ')
plt.scatter(wind_seasonal['JAS'].values,eke_seasonal['JAS'].values,marker='o',c='gold',edgecolors='none',label='JAS')
plt.scatter(wind_seasonal['OND'].values,eke_seasonal['OND'].values,marker='o',c='firebrick',edgecolors='none',label='OND')
plt.plot(wind_seasonal['JFM'].values,line.values,'--',color='gray')
leg=plt.legend(loc=1, scatterpoints = 1,bbox_to_anchor=(1.43,1),prop={'size':16})
plt.xlabel('U [m/s]',fontsize=20)
plt.ylabel(r'EKE [$\mathrm{{m}^2/{s}^2}$]',fontsize=20)
plt.tick_params(labelsize=20) 
plt.savefig(dir2+'/Figures/Atlantic_Area2_seasonal_relation2.png',\
            dpi=300,transparent=True,bbox_inches='tight',bbox_extra_artists=(leg,))
plt.close('all')

