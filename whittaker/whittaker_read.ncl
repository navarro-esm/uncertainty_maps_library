; Whittaker (part01)
;----------------------------------------------------------------------------------
; Author: Andres Navarro
; e-mail:andres.navarro@uclm.es
;
; Create intermediate files (pr.csv, ts.csv) to be read by the R script (whittaker.R)
;----------------------------------------------------------------------------------

debug=False    ;True or False
;modelname="Z-CRUTS"
;case="- Historical"
modelname=getenv("MODELNAME")
case=getenv("CASE")
print(modelname)
case2 = str_split(case, "-")  ;strip after -
case3 =str_left_strip(case2)  ;remove white spaces
case4=str_lower(case3)      ;lowercase

if (case2 .eq. " Historical") then
  if(modelname .eq. "Z-CRUTS") then
    fts=addfile("ts_amon_Z-CRUTS_historical_r1i1p1_190101-201412.nc","r")
    fpr=addfile("pr_amon_Z-CRUTS_historical_r1i1p1_190101-201412.nc","r")
    pr=fpr->pr(948:,:,:)
    ts=fts->ts(948:,:,:)
  else
    ftsin="ts_amon_"+modelname+"_historical_r1i1p1_185001-201412.nc"
    fprin="pr_amon_"+modelname+"_historical_r1i1p1_185001-201412.nc"
    fts=addfile(ftsin,"r")
    fpr=addfile(fprin,"r")
    pr=fpr->pr(1560:,:,:)
    ts=fts->ts(1560:,:,:)
  end if
else
  ftsin="ts_amon_"+modelname+"_"+case4+"_r1i1p1_201501-210012.nc"
  fprin="pr_amon_"+modelname+"_"+case4+"_r1i1p1_201501-210012.nc"
  fts=addfile(ftsin,"r")
  fpr=addfile(fprin,"r")
  pr=fpr->pr(:,:,:)
  ts=fts->ts(:,:,:)
end if


;----------------------------------------------------------------------------------
;Read CRU data to mask OCN & Antarctica
;----------------------------------------------------------------------------------
inmaskpath=addfile("ts_amon_Z-CRUTS_historical_r1i1p1_190101-201412.nc","r")
inmask=inmaskpath->ts(948,:,:)
inmask@_FillValue=-9999
;------------------------------------

pr=pr*2592000 ;to mm/month

tsClm = clmMonTLL( ts )  ;monthly climatology
tsClm@_FillValue=-9999
tsClm= tsClm-273.15
prClm = clmMonTLL(pr)    ;monthly climatology
prClm@_FillValue=-9999


;weighting
nlon=  360
nlat = 180
gwt  = latGauWgt(nlat, "lat", "gaussian weights", "")


;----------------------------------------------------------------------------------
;Annual precipitation (APP)
;----------------------------------------------------------------------------------
P=new((/nlat,nlon/),"float") ;annual precipitation
P=dim_sum_n(prClm,0) ;annual Precipitation
P@_FillValue=-9999
APP=dim_sum_n(prClm,0) ;annual precipitation (APP)
;APP=APP/10            ;from mm/year to cm/year
APP@_FillValue=-9999
APP!0="lat"
APP!1="lon"
APP&lat=ts&lat
APP&lon=ts&lon
APP=where(APP .eq. 0,0.001,APP) ;avoid div. by zero

;----------------------------------------------------------------------------------
;Annual mean Temperature (AT)
;----------------------------------------------------------------------------------
AT=new((/nlat,nlon/),"float")
AT@_FillValue=-9999

TTdummy=new((/12,nlat,nlon/),"float")     ;temporal temperature (month,lat,lon)
TTdummy=tsClm
TTdummy@_FillValue=-9999

AT=dim_sum_n(TTdummy/12,0) ;annual mean temperature  (AT)
AT!0="lat"
AT!1="lon"
AT&lat=ts&lat
AT&lon=ts&lon


;----------------------------------------------------------------------------------
;write to csv (APP & AT) to be read by the r script (whittaker.R)-> part02
;units mm/yr and ºC
;----------------------------------------------------------------------------------
asciiwrite("ts.csv",ndtooned(AT))
asciiwrite("pr.csv",ndtooned(APP))
