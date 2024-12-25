;Whittaker (part03)
;----------------------------------------------------------------------------------
; Author: Andres Navarro
; e-mail:andres.navarro@uclm.es
;----------------------------------------------------------------------------------
; IN: CSV file from Whittaker part02 (whittaker.csv)
;     Dimensions: 3x180x360 (APP,AT,Class).
;
;
; OUT: Whittaker's biomes (Dims. 180x360)
;     netcdf file (x1)
;     bil files (x1)
;     hdr files (x1)
;     Plot (x1)
;----------------------------------------------------------------------------------

debug=False    ;True or False
;modelname="Z-CRUTS"
;case="- Historical"
modelname=getenv("MODELNAME")       ;modelname
case=getenv("CASE")                 ;Historical, SSP126, SSP245,SSP585
print(modelname)
cr = tochar(10)


case2 = str_split(case, "-")  ;strip after -
case3 =str_left_strip(case2)  ;remove white spaces
case4=str_lower(case3)      ;lowercase

;----------------------------------------------------------------------------------
;Reading csv from Whittaker part02
;----------------------------------------------------------------------------------
;color code
; 'Tundra'(1),'Boreal forest'(2),'Temperate grassland/desert'(3),'Woodland/shrubland' (4),
; 'Temperate seasonal forest'(5),'Temperate rain forest' (6),'Tropical rain forest'(7)
; 'Tropical seasonal forest/savanna'(8), 'Subtropical desert' (9)
; code biomes 1= T, 2= BF, 3=TG-D, 4=W-Sh, 5=TsF, 6=TF, 7=TrF, 8=TsF, 9=SbTD
;----------------------------------------------------------------------------------
whit_in=asciiread("whittaker.csv",(/64800,3/),"float")
whit1d=whit_in(:,2)   ;only whittaker data
ts1d=whit_in(:,0)     ; only ts
pr1d=whit_in(:,1)     ; only pr
whit1d@_FillValue=-9999
ts1d@_FillValue=-9999
pr1d@_FillValue=-9999

;reshape to 180x360
whit1d:=toint(whit1d) ;convert to integer
whitdata=onedtond(whit1d,(/180,360/))
tsdata=onedtond(ts1d,(/180,360/))
prdata=onedtond(pr1d,(/180,360/))

;----------------------------------------------------------------------------------
;Read CRU data to mask OCN in plots
;----------------------------------------------------------------------------------
inmaskpath=addfile("ts_amon_Z-CRUTS_historical_r1i1p1_190101-201412.nc","r")
inmask=inmaskpath->ts(948,:,:)
inmask@_FillValue=-9999

;----------------------------------------------------------------------------------
;Defining coordinates
;----------------------------------------------------------------------------------
whitdata@_FillValue=-9999
whitdata!0="lat"
whitdata!1="lon"
whitdata&lat=inmask&lat
whitdata&lon=inmask&lon
whitdata@name="Whittaker"
whitdata@long_name="Whittaker biomes"
whitdata@comment="Code 1: Tundra. Code 2: Boreal forest. Code 3: Temperate grassland/desert." +cr+ \
                "Code 4: Woodland/Shrubland. Code 5: Temperate seasonal forest. Code 6: Temperate rain forest." +cr+  \
                "Code 7: Tropical rain forest. Code 8: Tropical seasonal forest/Savanna. Code 9: Subtropical desert."


;precip data
prdata!0="lat"
prdata!1="lon"
prdata&lat=inmask&lat
prdata&lon=inmask&lon
prdata=lonFlip(prdata) ;flip from 0-360 to -180 180
delete(prdata@lonFlip)  ; delete lonflip attribute
prdata@name="precipitation"
prdata@long_name="Annual precipitation"
prdata@units="cm/yr"


;temperature data
tsdata!0="lat"
tsdata!1="lon"
tsdata&lat=inmask&lat
tsdata&lon=inmask&lon
tsdata=lonFlip(tsdata) ;flip from 0-360 to -180 180
delete(tsdata@lonFlip)  ; delete lonflip attribute
tsdata@name="Temperature"
tsdata@long_name="Mean annual temperature"
tsdata@units="Degrees Celsius"


;----------------------------------------------------------------------------------
;plotting parameters
;----------------------------------------------------------------------------------
;masking data with cru mask
whitdata=where(ismissing(inmask),-9999,whitdata)
whitdata=lonFlip(whitdata) ;flip from 0-360 to -180 180
delete(whitdata@lonFlip)  ; delete lonflip attribute

wks = gsn_open_wks("pdf", modelname)

  ;; Plot resources
  res = True
  res@cnFillOn = True
  res@cnLinesOn = False
  ;res@cnFillMode = "RasterFill"
  res@cnFillMode = "CellFill"
  ;res@cnFillMode = "AreaFill"
  ;res@cnFixFillBleed =True
  ;res@cnRasterSmoothingOn = True
  res@gsnAddCyclic = True
  res@lbLabelBarOn   = False        ;quitar labelbar
  res@cnLineLabelsOn=False
  ;for robinson  activate the following two lines and disable pmTickMarkDisplayMode
  res@mpProjection="Robinson"
  res@mpGridAndLimbOn       =  False
  res@mpPerimOn =   False
  res@mpLimitMode ="LatLon"           ;Limitar Latitude (Remove Antarctica)
  res@mpMinLatF   =-60.0
  res@mpMaxLatF   = 90.0

  ;; Map
  long_name = (/"Whittaker ", modelname, case/)
  whitdata@long_name = str_concat(long_name)
  res@mpDataBaseVersion = "MediumRes"
  res@mpOutlineOn = True
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnFillPalette = "whittaker"
  labels = (/"T","BF","TG-D","W-Sh","TsF","TF","TrF","TsF","SbTD"/)
  res@lbLabelStrings = labels
  res@cnLevels = (/1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5/)
  res@cnFillColors = (/10,6,3,4,5,7,8,9,2,1/)
  res@lbLabelAlignment = "BoxCenters"
  res@lbOrientation = "vertical"
  res@lbAutoManage = True

 map = gsn_csm_contour_map(wks,whitdata,res); plot HLZ


;----------------------------------------------------------------------------------
;write to netcdf file
;----------------------------------------------------------------------------------
qq=addfile("whittaker_class_"+case4+"_"+modelname+".nc","c")
whitdata@long_name="Whittaker biomes"
qq->APP=prdata
qq->AT=tsdata
qq->whittaker=whitdata

;----------------------------------------------------------------------------------
;write to bil file (and .hdr)
;----------------------------------------------------------------------------------
cbinwrite("whittaker_class_"+case4+"_"+modelname+".bil", whitdata(::-1,:))  ;reverse lat

;adding data to the hdr file
strings = (/"BYTEORDER I","LAYOUT BIL","NBANDS 1","NROWS 180","NCOLS 360", \
"ULXMAP -180.0","ULYMAP 89.5","NBITS 32","NODATA -9999","PROJECTION GEOGRAPHIC", \
"DATUM WGS84","SPHEROID EPSG4326"/)
alist   = [/strings/]
write_table("whittaker_class_"+case4+"_"+modelname+".hdr", "w", alist, "%s") ;write .hdr
