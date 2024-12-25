#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# Scatter plot
# input: (1) precipitation.txt; (2) temperature.txt; (3) modelnames.txt
#        (4) Kappa-scores.txt
# output: (1) Scatter plot. y=precip(r2) x=class-agreement(kappa) z=temp(r2).
#         (2) Table of the best models (top-10) with the modelname, kappa, r2 pr, r2 ts.
# Includes zoom-in
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

from scipy import ndimage
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from netCDF4 import Dataset

#---------------------------------------------------------------------------------------------------------------
#reading data (size: 1x54 [kappa] 1x57 [precipitation and temperature])
#---------------------------------------------------------------------------------------------------------------
CNAME='whittaker'            #change the name (koppen, holdridge, whittaker, thornthwaite)

inputpr='/data/precipitation/pattern_corr_pramean_cc.txt'      # precipitation scores
inputts='/data/temperature/pattern_corr_tsamean_cc.txt'       #temperature scores
inputlist='/data/Z_filename_pr.txt'     #modelname list

#---------------------------------------------------------------------------------------------------------------
# Read model names (read from precipitation filename. Includes optimal ensembles [Holdridge, Thornthwaite, Koppen, Whittaker])
#---------------------------------------------------------------------------------------------------------------
fln_mod=np.genfromtxt(inputlist, delimiter='',dtype='str')

if CNAME=='koppen':
    fln_mod=np.delete(fln_mod, (53), axis=0)        #remove T10 holdgridge's ensemble
    fln_mod=np.delete(fln_mod, (54), axis=0)        #remove T10 thornthwaite's ensemble
    fln_mod=np.delete(fln_mod, (54), axis=0)        #remove T10 whittaker's ensemble
elif CNAME=='holdridge':
    fln_mod=np.delete(fln_mod, (54), axis=0)        #remove T10 koppen's ensemble
    fln_mod=np.delete(fln_mod, (54), axis=0)        #remove T10 thornthwaite's ensemble
    fln_mod=np.delete(fln_mod, (54), axis=0)        #remove T10 whittaker's ensemble
elif CNAME=='whittaker':
    fln_mod=np.delete(fln_mod, (53), axis=0)        #remove T10 holdridge's ensemble
    fln_mod=np.delete(fln_mod, (53), axis=0)        #remove T10 koppen's ensemble
    fln_mod=np.delete(fln_mod, (53), axis=0)        #remove T10 thornthwaite's ensemble
else:
    fln_mod=np.delete(fln_mod, (53), axis=0)        #remove T10 holdridge
    fln_mod=np.delete(fln_mod, (53), axis=0)        #remove T10 koppen
    fln_mod=np.delete(fln_mod, (54), axis=0)        #remove T10 whittaker

#---------------------------------------------------------------------------------------------------------------
# read climate classification (include only k for individual models, MME (52) and its own T10 ensemble
#---------------------------------------------------------------------------------------------------------------
if CNAME=='koppen':
    inputcclas='/data/koppen/kappa-koppen.txt'
elif CNAME=='holdridge':
    inputcclas='/data/holdridge/kappa-holdridge.txt'
elif CNAME=='whittaker':
    inputcclas='/data/whittaker/kappa-whittaker.txt'
else:
    inputcclas='/data/thornthwaite/kappa-thornthwaite.txt'

datin_cclas=np.genfromtxt(inputcclas, delimiter='',dtype='float')


#---------------------------------------------------------------------------------------------------------------
# read R2 for precipitation and temperature
#---------------------------------------------------------------------------------------------------------------
datin_pr=np.genfromtxt(inputpr, delimiter='',dtype='float')
datin_pr=datin_pr**2        #R2
datin_ts=np.genfromtxt(inputts, delimiter='',dtype='float')
datin_ts=datin_ts**2        #R2

if CNAME=='koppen':
    datin_pr=np.delete(datin_pr, (53), axis=0)        #remove holdgridge pr
    datin_pr=np.delete(datin_pr, (54), axis=0)        #remove thornthwaite pr
    datin_pr=np.delete(datin_pr, (54), axis=0)        #remove whittaker pr
    datin_ts=np.delete(datin_ts, (53), axis=0)        #remove holdgridge ts
    datin_ts=np.delete(datin_ts, (54), axis=0)        #remove thornthwaite ts
    datin_ts=np.delete(datin_ts, (54), axis=0)        #remove whittaker ts
elif CNAME=='holdridge':
    datin_pr=np.delete(datin_pr, (54), axis=0)        #remove koppen pr
    datin_pr=np.delete(datin_pr, (54), axis=0)        #remove thornthwaite pr
    datin_pr=np.delete(datin_pr, (54), axis=0)        #remove whittaker pr
    datin_ts=np.delete(datin_ts, (54), axis=0)        #remove koppen ts
    datin_ts=np.delete(datin_ts, (54), axis=0)        #remove thornthwaite ts
    datin_ts=np.delete(datin_ts, (54), axis=0)        #remove whittaker ts
elif CNAME=='whittaker':
    datin_pr=np.delete(datin_pr, (53), axis=0)        #remove holdridge pr
    datin_pr=np.delete(datin_pr, (53), axis=0)        #remove koppen pr
    datin_pr=np.delete(datin_pr, (53), axis=0)        #remove thornthwaite pr
    datin_ts=np.delete(datin_ts, (53), axis=0)        #remove holdridge ts
    datin_ts=np.delete(datin_ts, (53), axis=0)        #remove koppen ts
    datin_ts=np.delete(datin_ts, (53), axis=0)        #remove thornthwaite ts
else:
    datin_pr=np.delete(datin_pr, (53), axis=0)        #remove holdridge pr
    datin_pr=np.delete(datin_pr, (53), axis=0)        #remove koppen pr
    datin_pr=np.delete(datin_pr, (54), axis=0)        #remove whittaker pr
    datin_ts=np.delete(datin_ts, (53), axis=0)        #remove holdridge ts
    datin_ts=np.delete(datin_ts, (53), axis=0)        #remove koppen ts
    datin_ts=np.delete(datin_ts, (54), axis=0)        #remove whittaker ts


#---------------------------------------------------------------------------------------------------------------
# Labeling dots
#---------------------------------------------------------------------------------------------------------------
dsize=datin_pr.shape
dsize=dsize[0]

#label for scatters
datin_cnames=[]
for i in range(1,dsize+1):
    datin_cnames.append(str(i))

#---------------------------------------------------------------------------------------------------------------
# put all data into single array
#---------------------------------------------------------------------------------------------------------------
data= np.zeros((dsize,3),dtype=float)
data[:,0]=datin_cclas
data[:,1]=datin_pr
data[:,2]=datin_ts

#---------------------------------------------------------------------------------------------------------------
# make a copy of all data for the last table (exclude info about ensembles)
#---------------------------------------------------------------------------------------------------------------
data2=data
data2=np.delete(data2, (52), axis=0)        #remove MME (52)
data2=np.delete(data2, (52), axis=0)        #remove T10 ensemble

#---------------------------------------------------------------------------------------------------------------
# descriptive stats (take into account only individual models, not MME (52) and optimal ensemble (t10))
#---------------------------------------------------------------------------------------------------------------
cclas_median=np.nanmedian(data[0:51,0]) #dotted lines
pr_median=np.nanmedian(data[0:51,1])

cclas_q1=np.nanquantile(data[0:51,0],0.25) #Q1
cclas_q3=np.nanquantile(data[0:51,0],0.75) #Q3

pr_q1=np.nanquantile(data[0:51,1],0.25)
pr_q3=np.nanquantile(data[0:51,1],0.75)

#---------------------------------------------------------------------------------------------------------------
# plotting data (include all individual models and MME (52) and optimal ensemble (t10))
#---------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(7., 7.))
ax.set_xlabel('Class agreement [k]', fontsize=25,labelpad=7.5)
ax.set_ylabel('Precipitation [R²]', fontsize=25,labelpad=7.5)
ax.set_title('', fontsize=30)
plt.xticks([0.0,0.25,0.5,0.75,1.0], fontsize = 15)
plt.yticks([0.0,0.25,0.5,0.75,1.0], fontsize = 15)
plt.xlim(0.0, 1.0)
plt.ylim(0.0,1.0)

#line plots
plt.plot([cclas_median,cclas_median],[1,0], linewidth=1, color='black' ,linestyle=':')
plt.plot([1,0],[pr_median,pr_median], linewidth=1, color='black' ,linestyle=':')


#plotting scatters
PCM=ax.scatter(data[:,0], data[:,1],c=data[:,2],cmap='viridis_r',edgecolors='none',vmin=0.9,vmax=1.0, s=60)

#annotating
for i, txt in enumerate(datin_cnames):
    ax.annotate(txt, (data[i,0], data[i,1]))


#adding colorbar
fig_coord = [0.17,0.78,0.25,0.025]

cbar_ax = fig.add_axes(fig_coord)

cbar=plt.colorbar(PCM,orientation='horizontal',cax=cbar_ax)
cbar.ax.set_title("Temperature [R²]", pad=7.5 , fontsize=15)

cbar.ax.tick_params(labelsize=10)

#-----------------------------
# add inner plot
#-----------------------------
if CNAME=='koppen' or CNAME=='whittaker' or CNAME=='holdridge':
    axin = ax.inset_axes([0.10, 0.05, 0.43, 0.43])# Plot the data on the inset axis and zoom in on the important part
else:
    axin = ax.inset_axes([0.55, 0.05, 0.43, 0.43])# Plot the data on the inset axis and zoom in on the important part

axin.scatter(data[:,0], data[:,1],c=data[:,2], cmap='viridis_r',edgecolors='none',vmin=0.9,vmax=1.0, s=60)
axin.set_xlim(cclas_median, np.amax(data[:,0]+0.01))
axin.set_ylim(pr_median, np.amax(data[:,1]+0.01))# Add the lines to indicate where the inset axis is coming from
ax.indicate_inset_zoom(axin)

plt.draw()

plt.savefig('scatter-'+CNAME+'.pdf', format='pdf')

#-------------------------------------------------------------------
#get data from top 10 models
#-------------------------------------------------------------------
kcoeff=data2[:,0]
prcoeff=data2[:,1]
tscoeff=data2[:,2]

#-------------------------------------------------------------------
#check if top 10 models are located in upper-right quadrant (if not, substitute
# with other models using precipitation as reference)
#-------------------------------------------------------------------


kcoeffidx = np.argwhere( kcoeff >= cclas_median)            #idx kappa great than median
prcoeffidx =np.argwhere( prcoeff >= pr_median)              #idx r2 precip great than median

urq = np.unique(kcoeffidx[np.isin(kcoeffidx, prcoeffidx)])      #compare both and get shared idx (pr r2>median and k>median)

kselcoeff=kcoeff[urq]           #kappa of models with score great than meadian
prselcoeff=prcoeff[urq]         #r2 of models with score great than meadian
tsselcoeff=tscoeff[urq]         #r2 of ts (same models as r2 pr. Only for including rscore into a new table)
topselmod=fln_mod[urq]          #model names with score great than meadian

ind = np.argpartition(kselcoeff, -10)[-10:]        #get idx of selected top 10 models
ind=np.sort(ind)            #sort by index
topcoeff=kselcoeff[ind]    #selected kappa values
topprcoeff=prselcoeff[ind]      #name of selected top10 models (pr)
toptscoeff=tsselcoeff[ind]      #name of selected top10 models (ts)
topmod=topselmod[ind]      #name of selected top10 models

isin_result = np.isin(fln_mod, topmod)
shared_idx = []
for i, val in enumerate(isin_result):
    if val:
        shared_idx.append(i)

shared_idx=np.array(shared_idx)
topnumber=shared_idx+1         #get the number from the plot


#-------------------------------------------------------------------
#write top10 list to txt
#-------------------------------------------------------------------
names  = topmod
floats1 = topcoeff
floats2 = topprcoeff
floats3 = toptscoeff
ab = np.zeros(names.size, dtype=[('var1', 'U16'), ('var2', float),('var3', float),('var4', float)  ])
ab['var1'] = names
ab['var2'] = floats1
ab['var3'] = floats2
ab['var4'] = floats3

np.savetxt('top10-models_'+CNAME+'.txt', ab, fmt="%10s %10.3f" "%10.3f" "%10.3f")
