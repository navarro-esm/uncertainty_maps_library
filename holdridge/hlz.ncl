;---------------------------------------------------------------------
; Holdridge Life Zone Scheme
; Based on Sisneros et al. (2011) Visualizing Life Zone Boundary Sensitivities
; Across Climate Models and Temporal Spans.
; Procedia Computer Science 4, pp: 1582-1591.
; Aggregated into 13 categories using Monserud & Leemans (1992)
;
; Author: Andres Navarro
; e-mail: andres.navarro@uclm.es
;---------------------------------------------------------------------
;IN: netcdf files of precipitation and temperature for each CMIP6 model
; e.g. ts_amon_CESM2_historical_r1i1p1_185001-201412.nc
;      pr_amon_CESM2_historical_r1i1p1_185001-201412.nc
;
;OUT: Holdridge's Life Zones (Dims. 180x360)
;     netcdf file (x1)
;     bil file (x1)
;     hdr file (x1)
;     Plot (x4) HLZ,APP,ABT,PER
;---------------------------------------------------------------------



debug=False    ;True or False
modelname=getenv("MODELNAME")
case=getenv("CASE")
print(modelname)

case2 = str_split(case, "-")  ;strip after -
case3 =str_left_strip(case2)  ;remove white spaces
case4=str_lower(case3)      ;lowercase
cr = tochar(10)

;---------------------------------------------------------------------
; Read netcdf files
;---------------------------------------------------------------------

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

;---------------------------------------------------------------------
;Read CRU data to mask OCN
;---------------------------------------------------------------------
inmaskpath=addfile("ts_amon_Z-CRUTS_historical_r1i1p1_190101-201412.nc","r")
inmask=inmaskpath->ts(948,:,:)
inmask@_FillValue=-9999
;---------------------------------------------------------------------

pr=pr*2592000 ;to mm/month

tsClm = clmMonTLL( ts )  ;monthly climatology (extraer antes 1980-2014)
tsClm@_FillValue=-9999
tsClm= tsClm-273.15
prClm = clmMonTLL(pr)    ;monthly climatology (extraer antes 1980-2014)
prClm@_FillValue=-9999

;---------------------------------------------------------------------
;Computing:
;       Annual precipitation (APP)
;       Annual mean Biotemperature (ABT)
;       Potential evapotranspiration (PER)
;---------------------------------------------------------------------


;weighting
nlon=  360
nlat = 180
gwt  = latGauWgt(nlat, "lat", "gaussian weights", "")


;---------------------------------------------------------------------
;Annual precipitation (APP)
;---------------------------------------------------------------------
P=new((/nlat,nlon/),"float") ;annual precipitation
P=dim_sum_n(prClm,0) ;annual Precipitation
P@_FillValue=-9999
APP=dim_sum_n(prClm,0) ;annual precipitation (APP)

;Coordinates
APP!0="lat"
APP!1="lon"
APP&lat=ts&lat
APP&lon=ts&lon
APP@_FillValue=-9999
APP@name="Precipitation"
APP@long_name="Annual precipitation"
APP@units="mm/yr"
APP=where(APP .eq. 0,0.001,APP) ;avoid div. by zero

;---------------------------------------------------------------------
;Annual mean Biotemperature (ABT)
;---------------------------------------------------------------------
ABT=new((/nlat,nlon/),"float") ;annual mean biotemperature
ABT@_FillValue=-9999

TTdummy=new((/12,nlat,nlon/),"float")     ;temporal biotemperature (month,lat,lon)
TTdummy=tsClm
TTdummy@_FillValue=-9999
TTdummy=where(TTdummy .le. 0.0, 0.0, TTdummy)   ; set <0.0 degrees to 0 degrees
TTdummy=where(TTdummy .ge. 30.0, 30.0, TTdummy)  ; set >30 degree to 0 degrees

ABT=dim_sum_n(TTdummy/12,0) ;annual mean biotemperature (ABT)

;Coordinates
ABT!0="lat"
ABT!1="lon"
ABT&lat=ts&lat
ABT&lon=ts&lon
ABT@name="Biotemperature"
ABT@long_name="Mean annual biotemperature"
ABT@units="Degrees Celsius"

;---------------------------------------------------------------------
;Potential evapotranspiration ratio (PER)
;---------------------------------------------------------------------
APE=new((/nlat,nlon/),"float") ;Annual potential evapotranspiration
APE@_FillValue=-9999
APE=ABT*58.93

PER=new((/nlat,nlon/),"float") ;Potential evapotranspiration ratio
PER@_FillValue=-9999
do ilat=0,nlat-1
  do jlon=0,nlon-1
    if (ismissing(APE(ilat,jlon))) then
      PER(ilat,jlon)=PER@_FillValue
    else
      PER(ilat,jlon) = APE(ilat,jlon)
      PER(ilat,jlon) = PER(ilat,jlon)/APP(ilat,jlon)
    end if
  end do
end do

;Coordinates
PER!0="lat"
PER!1="lon"
PER&lat=ts&lat
PER&lon=ts&lon
PER@name="PER"
PER@long_name="Potential evapotranspiration ratio (PET/APP)"
PER@units=""

;---------------------------------------------------------------------
;Initialize HLZ array
;---------------------------------------------------------------------
HLZCL3=new((/nlat,nlon/),"integer")
HLZCL3!0="lat"
HLZCL3!1="lon"
HLZCL3&lat=ts&lat
HLZCL3&lon=ts&lon
HLZCL3@_FillValue=-9999
HLZCL3=-9999


; Usamos categorizacion Monserud & Leemans 1992 modificada (13 cats.)
; Version usando los centroides sacados de Sisneros et al. 2011.
; computing the euclidean distance sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

APPr=(/88.39,176.78,353.55,88.39,176.78,353.55,707.11,88.39,176.78,353.55,707.11,1414.21, \
      88.39,176.78,353.55,707.11,1414.21,2828.43,88.39,176.78,353.55,707.11,1414.21,2828.43,\
      5656.85,88.39,176.78,353.55,707.11,1414.21,2828.43,5656.85,11313.71/)                      ;centroids Precipitation
ABTr=(/0.0,0.0,0.0,2.12,2.12,2.12,2.12,4.24,4.24,4.24,4.24,4.24,8.49,8.49,8.49,8.49,8.49,8.49, \
      16.97,16.97,16.97,16.97,16.97,16.97,16.97,26.83,26.83,26.83,26.83,26.83,26.83,26.83,26.83/)  ;centroids biotemperature
PERr=(/0.71,0.35,0.18,1.41,0.71,0.35,0.18,2.83,1.41,0.71,0.35,0.18,5.66,2.83,1.41,0.71,0.35, \
      0.18, 11.31,5.66,2.83,1.41,0.71,0.35,0.18,22.63,11.31,5.66,2.83,1.41,0.71,0.35,0.18/)        ;Centroids Pot. evapotrans
HLZCodes=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33/) ;codes
HLZNames=(/"PD","PD","PD","DT","MT","WT","RT","BD","BDS","BMF","BWF","BRF","CD","CDS","STP","CMF", \
         "CWF","CRF","SbD","SbDS","ThS","SbDF","SbMF","SbWF","SbRF","TD","TDS","ThW","TVDF","TDF",\
         "TMF","TWF","TRF"/)                                                                           ;Names HLZ
  ; Computing Euclidean distance and assign climate category
 if (modelname .eq. "Z-CRUTS") then                         ;Check special case (CRU)
   do ilat=0,nlat-1
     do jlon=0,nlon-1
       if (ismissing(APP(ilat,jlon))) then
         HLZCL3(ilat,jlon)=-9999
       else
         ncentroids=33
         edist=new((/ncentroids/),"float")
         do a=0,ncentroids-1
           x=(APP(ilat,jlon)-APPr(a))^2
           y=(ABT(ilat,jlon)-ABTr(a))^2
           z=(PER(ilat,jlon)-PERr(a))^2
           edist(a)=sqrt((x+y+z))
         end do
         edismin=min(edist)
         idx=ind(edist .eq. edismin)
         HLZCL3(ilat,jlon)=HLZCodes(idx(0))     ;use idx(0)
         delete(edist)
         delete(edismin)
         delete(idx)
       end if
     end do
   end do
 else                                                         ; Do this for other cases
  do ilat=0,nlat-1
    do jlon=0,nlon-1
      ncentroids=33
      edist=new((/ncentroids/),"float")
      do a=0,ncentroids-1
        x=(APP(ilat,jlon)-APPr(a))^2
        y=(ABT(ilat,jlon)-ABTr(a))^2
        z=(PER(ilat,jlon)-PERr(a))^2
        edist(a)=sqrt((x+y+z))
      end do
      edismin=min(edist)
      idx=ind(edist .eq. edismin)
      HLZCL3(ilat,jlon)=HLZCodes(idx(0))     ;use idx(0)
      delete(edist)
      delete(edismin)
      delete(idx)
    end do
  end do
 end if

;re-assign values
HLZCL3!0="lat"
HLZCL3!1="lon"
HLZCL3&lat=ts&lat
HLZCL3&lon=ts&lon
HLZCL3@_FillValue=-9999

;---------------------------------------------------------------------
;Create new HLZ to save 13 categories instead of 33 of standard Holdridge
;---------------------------------------------------------------------
HLZCL4=HLZCL3
HLZCL4!0="lat"
HLZCL4!1="lon"
HLZCL4&lat=ts&lat
HLZCL4&lon=ts&lon
HLZCL4@_FillValue=-9999
HLZCL4=-9999
HLZCL4@name="HLZ"
HLZCL4@long_name="Holdridge Life Zones"
HLZCL4@comment= "Code 1: Tundra. Code 2: Cold Parklands. Code 3: Forest Tundra." +cr+ \
                "Code 4: Boreal Forest. Code 5: Cool Desert. Code 6: Steppe." +cr+  \
                "Code 7: Cool Forest. Code 8: Hot Desert. Code 9: Chaparral." +cr+  \
                "Code 10: Temperate Forest. Code 11: Tropical Semi-arid." +cr+  \
                "Code 12: Tropical dry forest. Code 13: Tropical Rain Forest."

HLZCL4=where(HLZCL3 .eq. 1 .or. HLZCL3 .eq. 2 .or.HLZCL3 .eq. 3, 1, HLZCL4)     ; tundra
HLZCL4=where(HLZCL3 .eq. 4 .or. HLZCL3 .eq. 8 .or.HLZCL3 .eq. 9, 2, HLZCL4)     ; Cold Parklands
HLZCL4=where(HLZCL3 .eq. 5 .or. HLZCL3 .eq. 6 .or.HLZCL3 .eq. 7, 3, HLZCL4)     ; Forest tundra
HLZCL4=where(HLZCL3 .eq. 10 .or. HLZCL3 .eq. 11 .or.HLZCL3 .eq. 12, 4, HLZCL4)  ; Boreal Forest
HLZCL4=where(HLZCL3 .eq. 13 .or. HLZCL3 .eq. 14, 5, HLZCL4)                     ; Cool Desert
HLZCL4=where(HLZCL3 .eq. 15 , 6, HLZCL4)                                        ; Steppe
HLZCL4=where(HLZCL3 .eq. 16 .or. HLZCL3 .eq. 17 .or.HLZCL3 .eq. 18, 7, HLZCL4)  ; Cool Forest
HLZCL4=where(HLZCL3 .eq. 19 .or. HLZCL3 .eq. 20  .or.HLZCL3 .eq. 26 .or. HLZCL3 .eq. 27, 8, HLZCL4)  ; Hot desert
HLZCL4=where(HLZCL3 .eq. 21 .or. HLZCL3 .eq. 22, 9, HLZCL4)                     ; Chapparal
HLZCL4=where(HLZCL3 .eq. 23 .or. HLZCL3 .eq. 24 .or.HLZCL3 .eq. 25, 10, HLZCL4)  ; Temperate forest
HLZCL4=where(HLZCL3 .eq. 28 .or. HLZCL3 .eq. 29, 11, HLZCL4)                     ; Tropical Semi-arid
HLZCL4=where(HLZCL3 .eq. 30 ,12 , HLZCL4)                                        ; Tropical dry forest
HLZCL4=where(HLZCL3 .eq. 31 .or. HLZCL3 .eq. 32 .or.HLZCL3 .eq. 33, 13, HLZCL4)  ; Tropical Rain forest

HLZCodes2=(/1,2,3,4,5,6,7,8,9,10,11,12,13/) ;codes
HLZNames2=(/"Tundra","Cold Parklands","Forest Tundra","Boreal Forest","Cool Desert","Steppe","Cool Forest","Hot Desert",\
            "Chaparral","Temperate Forest","Trop. Semi-arid","Trop. dry","Trop. Rain"/)


;Masking OCN for plotting
;---------------------------------------------------------------------
HLZCL4=where(ismissing(inmask),-9999,HLZCL4)



;---------------------------------------------------------------------
;plotting HLZ
;---------------------------------------------------------------------
HLZCL3pl=HLZCL4
long_name = (/"(13) Holdridge ", modelname, case/)
HLZCL3pl@long_name = str_concat(long_name)

wks = gsn_open_wks("pdf", modelname)

res = True
res@cnFillOn = True
res@cnLinesOn = False
;res@cnFillMode = "RasterFill"
res@cnFillMode = "CellFill"
;res@cnRasterSmoothingOn = True
res@gsnAddCyclic = True
res@lbLabelBarOn   = False
res@cnLineLabelsOn=False
;for robinson  activate the following two lines and disable pmTickMarkDisplayMode
res@mpProjection="Robinson"
res@mpGridAndLimbOn       =  False
res@mpPerimOn =   False
res@mpLimitMode ="LatLon"
res@mpMinLatF   =-60.0
res@mpMaxLatF   = 90.0

res@mpDataBaseVersion = "MediumRes"
res@mpOutlineOn = True
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnFillPalette = "material_navarro_01" ;custom palette
labels = (/"T","CP","FT","BF","CD","ST","CF","HD","CH","TF","TS","TD","TR"/)
res@lbLabelStrings = labels
res@cnLevels = (/1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5/)
res@cnFillColors = (/243,246,233,234,111,127,235,178,193,236,208,218,237/) ;new 15/07/22 (material navarro)
res@lbLabelAlignment = "BoxCenters"
res@lbOrientation = "vertical"
res@lbAutoManage = True

map = gsn_csm_contour_map(wks,HLZCL3pl,res); plot HLZ

;---------------------------------------------------------
;Plotting APP,ABT & PER
;---------------------------------------------------------
;APP
;---------------------------------------
;Masking OCN for plotting
APPpl=APP
APPpl=where(ismissing(inmask),-9999,APPpl)

wks2 = gsn_open_wks("pdf", "Z-APP_"+modelname)
res2 = True
res2@cnFillOn = True
res2@cnLinesOn = False
res2@cnFillMode = "RasterFill"
;res2@cnRasterSmoothingOn = True
res2@gsnAddCyclic = True
res2@lbLabelBarOn   = False        ;quitar labelbar
res2@cnLineLabelsOn=False
res2@mpProjection="Robinson"
res2@mpGridAndLimbOn       =  False
res2@mpPerimOn =   False
res2@mpLimitMode ="LatLon"
res2@mpMinLatF   =-60.0
res2@mpMaxLatF   = 90.0

long_name = (/"APP ", modelname, case/)
APPpl@long_name = str_concat(long_name)
res2@mpDataBaseVersion = "MediumRes"
res2@mpOutlineOn = True
;res2@pmTickMarkDisplayMode = "Always"
res2@cnLevelSelectionMode = "ExplicitLevels"
res2@cnFillPalette = "GMT_drywet"

labels2 = (/"0","125","250","500","1000","2000","4000","8000"/)
res2@lbLabelStrings = labels2
res2@cnLevels = (/0,125,250,500,1000,2000,4000,8000/)
res2@cnFillColors = (/2,6,9,17,26,33,40,49,61/)
res2@lbOrientation = "vertical"
res2@lbAutoManage = True
map2 = gsn_csm_contour_map(wks2,APPpl,res2); plot APP

;---------------------------------------
; ABT
;---------------------------------------
;Masking OCN for plotting
ABTpl=ABT
ABTpl=where(ismissing(inmask),-9999,ABTpl)

wks3 = gsn_open_wks("pdf", "Z-ABT_"+modelname)
res3 = True
res3@cnFillOn = True
res3@cnLinesOn = False
res3@cnFillMode = "RasterFill"
;res3@cnRasterSmoothingOn=True
res3@gsnAddCyclic = True
res3@lbLabelBarOn   = False
res3@cnLineLabelsOn=False
res3@mpProjection="Robinson"
res3@mpGridAndLimbOn       =  False
res3@mpPerimOn =   False
res3@mpLimitMode ="LatLon"
res3@mpMinLatF   =-60.0
res3@mpMaxLatF   = 90.0

long_name = (/"ABT ", modelname, case/)
ABTpl@long_name = str_concat(long_name)
res3@mpDataBaseVersion = "MediumRes"
res3@mpOutlineOn = True
;res3@pmTickMarkDisplayMode = "Always"
res3@cnLevelSelectionMode = "ExplicitLevels"
res3@cnFillPalette = "nice_gfdl"

labels3 = (/"0.375","0.75","1.5","3","6","12","24"/)
res3@lbLabelStrings = labels3
res3@cnLevels = (/0.375,0.75,1.5,3.0,6.0,12.0,24.0/)
res3@cnFillColors = (/49,65,81,97,113,146,154,171/)
res3@lbOrientation = "vertical"
res3@lbAutoManage = True
map3 = gsn_csm_contour_map(wks3,ABTpl,res3); plot APP

;---------------------------------------
; PER
;---------------------------------------
;Masking OCN for plotting
PERpl=PER
PERpl=where(ismissing(inmask),-9999,PERpl)
wks4 = gsn_open_wks("pdf", "Z-PER_"+modelname)
res4 = True
res4@cnFillOn = True
res4@cnLinesOn = False
res4@cnFillMode = "RasterFill"
;res4@cnRasterSmoothingOn=True
res4@gsnAddCyclic = True
res4@lbLabelBarOn   = False
res4@cnLineLabelsOn=False
res4@mpProjection="Robinson"
res4@mpGridAndLimbOn       =  False
res4@mpPerimOn =   False
res4@mpLimitMode ="LatLon"
res4@mpMinLatF   =-60.0
res4@mpMaxLatF   = 90.0

long_name = (/"PER ", modelname, case/)
PERpl@long_name = str_concat(long_name)
res4@mpDataBaseVersion = "MediumRes"
res4@mpOutlineOn = True
;res4@pmTickMarkDisplayMode = "Always"
res4@cnLevelSelectionMode = "ExplicitLevels"
res4@cnFillPalette = "OceanLakeLandSnow"

labels4 = (/"0.125","0.25","0.5","1.0","2.0","4.0","8.0","16.0","32.0"/)
res4@lbLabelStrings = labels4
res4@cnLevels = (/0.125,0.25,0.5,1.0,2.0,4.0,8.0,16.0,32.0/)
res4@cnFillColors = (/5,17,33,65,81,97,129,161,193,225/)
res4@lbOrientation = "vertical"
res4@lbAutoManage = True
map4 = gsn_csm_contour_map(wks4,PERpl,res4); plot APP

;---------------------------------------------------------------------
;Flip longitude to -180 +180
;---------------------------------------------------------------------
HLZCL4=lonFlip(HLZCL4)
delete(HLZCL4@lonFlip)
PER=lonFlip(PER)
delete(PER@lonFlip)
APP=lonFlip(APP)
delete(APP@lonFlip)
ABT=lonFlip(ABT)
delete(ABT@lonFlip)
;---------------------------------------------------------------------
;write to netcdf file
;---------------------------------------------------------------------
qq=addfile("holdridge_class_"+case4+"_"+modelname+".nc","c")
HLZCL4@long_name="Holdridge Life Zones"
qq->APP=APP
qq->ABT=ABT
qq->PER=PER
qq->HLZ=HLZCL4

;---------------------------------------------------------------------
;write to bil file (and .hdr)
;---------------------------------------------------------------------
cbinwrite("holdridge_class_"+case4+"_"+modelname+".bil", HLZCL4(::-1,:))  ;reverse lat

;adding data to the hdr file
strings = (/"BYTEORDER I","LAYOUT BIL","NBANDS 1","NROWS 180","NCOLS 360", \
"ULXMAP -180.0","ULYMAP 89.5","NBITS 32","NODATA -9999","PROJECTION GEOGRAPHIC", \
"DATUM WGS84","SPHEROID EPSG4326"/)
alist   = [/strings/]
write_table("holdridge_class_"+case4+"_"+modelname+".hdr", "w", alist, "%s") ;write .hdr
