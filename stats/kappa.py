#--------------------------------------------------------------------------------
# Kappa statistic
# Based on Cohen's Kappa
# Coefficient of Agreement for Nominal Scales. Educ. Psychol. Measurement 20, 37–46 (1960).
#
# 0 no agreement
# 1 perfect agreement
#--------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import sklearn.metrics
import os
import sys

#--------------------------------------------------------------------------------
# Filenames and open files
#--------------------------------------------------------------------------------
file1= sys.argv[1] 	# OUTER
file2= sys.argv[2] 	# INNER
case=sys.argv[3] # case
cclas=sys.argv[4] #climate classification


dtype = np.dtype([
    ("PosX", np.float32),
    ])

f = open(file1, "rb")
data1 = np.fromfile(f, dtype=dtype)
f2= open(file2, "rb")
data2=np.fromfile(f2, dtype=dtype)

#--------------------------------------------------------------------------------
#Reading latitude weights from external file
#--------------------------------------------------------------------------------
file3="/cell-area/cell_area_1d_180x360-km.txt"
m3=np.loadtxt(file3)

#--------------------------------------------------------------------------------
#Storing data into 1D array m1=reference m2=model
#--------------------------------------------------------------------------------
m1=data1.astype(int)	#convert to integer
m2=data2.astype(int)



#--------------------------------------------------------------------------------
#masking data: thornthwaite, whittaker, holdridge and koppen
#--------------------------------------------------------------------------------
m1masked =np.ma.masked_where(m1==-9999,m1) #mask no fill values
m2masked =np.ma.masked_where(m1==-9999,m2)
m1masked=np.ma.masked_where(m2==-9999,m1masked)
m2masked=np.ma.masked_where(m2==-9999,m2masked)
m3masked =np.ma.masked_where(m1==-9999,m3)          #Apply the same mask to latitude weights
m3masked=np.ma.masked_where(m2==-9999,m3masked)



m1a=m1masked.compressed()			#new array excluding nofill values
m2a=m2masked.compressed()

m3a=m3masked.compressed()       #selected pixels
m3asum=np.sum(m3a)              # sum of the area of selected pixels
m3a=m3a/m3asum                  # weights


#--------------------------------------------------------------------------------
# Apply's Kappa
#--------------------------------------------------------------------------------
k=sklearn.metrics.cohen_kappa_score(m1a,m2a, sample_weight=m3a)         #Cohen's Kappa

print(file1+" "+ file2+ " "+str(k))


#--------------------------------------------------------------------------------
# Save to file txt
#--------------------------------------------------------------------------------
l=[None]* 3
l[0]=file1
l[1]=file2
l[2]=k
# Writing data to a common txt file
du = pd.DataFrame(data=[l])
du.to_csv("kappa_results_SCIDAT_"+cclas+"_"+case+".txt", sep=' ', index=False, header=False, mode='a', encoding='utf-8',na_rep='NA', float_format='%1.5f')
