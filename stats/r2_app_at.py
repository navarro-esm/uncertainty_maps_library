#-----------------------------------------------
# Compute R^2 statistic using statmodels
#-----------------------------------------------

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.graphics.gofplots import ProbPlot
import os
import sys


#-----------------------------------------------
# Read files
#-----------------------------------------------
file1= sys.argv[1] 	# model filename
case=sys.argv[2] # case (historical, ssp126, ssp245)

# reference data (CRU)
#file2="/media/andres/WD5TB/BACKUP_WD_CLIMATE_CLASS/CMIP6_TAPI/NAT_SCIDAT/data/r_2-new/with-fillvalues/app/app_class_Z-CRUTS.txt" # for APP
file2="/media/andres/WD5TB/BACKUP_WD_CLIMATE_CLASS/CMIP6_TAPI/NAT_SCIDAT/data/r_2-new/with-fillvalues/at/at_class_Z-CRUTS.txt" # for AT


data1=pd.read_csv(file1,header=None,dtype='float')
data2=pd.read_csv(file2,header=None,dtype='float')


file3="cell_area_1d_180x360-km.txt"         #weight
m3=np.loadtxt(file3)

m1=np.asarray(data1, dtype='float')
m2=np.asarray(data2, dtype='float')
m1=m1[:,0]
m2=m2[:,0]


#-----------------------------------------------
# masking  fill values
#-----------------------------------------------
m1masked =np.ma.masked_where(m1==-9999,m1)
m2masked =np.ma.masked_where(m1==-9999,m2)

m3masked =np.ma.masked_where(m1==-9999,m3)
m3masked=np.ma.masked_where(m2==-9999,m3masked)


#-----------------------------------------------
# New arrays (excl. fillvalues)
#-----------------------------------------------
m1a=m1masked.compressed()			#new array without fillvalues
m2a=m2masked.compressed()

m3a=m3masked.compressed()
m3asum=np.sum(m3a)
m3a=m3a/m3asum


xfill=np.where(m1a ==-9999)
yfill=np.where(m2a ==-9999)

m1a[yfill]=-9999
m2a[xfill]=-9999

#-----------------------------------------------
# statsmodels (weighted)
#-----------------------------------------------
m1a1 = sm.add_constant(m1a)
mod_wls = sm.WLS(m2a,m1a1, m3a)
mod_res = mod_wls.fit()
mod_res_r2=mod_res.rsquared


#-----------------------------------------------
#writing results to txt file
#-----------------------------------------------
l2=[None]
l2[0]=mod_res_r2
du2 = pd.DataFrame(data=[l2])

# for APP
#du2.to_csv("rsquared_app_results_"+case+".txt", sep=' ', index=False, header=False, mode='a', encoding='utf-8',na_rep='NA', float_format='%1.3f')

# for AT
du2.to_csv("rsquared_abt_results_"+case+".txt", sep=' ', index=False, header=False, mode='a', encoding='utf-8',na_rep='NA', float_format='%1.3f')
