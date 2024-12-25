;---------------------------------------------------------------------
; Modified version of Thornthwaite climate classification
; Based on Feddema (2005) A revised Thornthwaite-type global climate classification
; Physical Geography, 26, pp 442-466.
; Author: Andres Navarro
; e-mail:andres.navarro@uclm.es
;---------------------------------------------------------------------
;IN: NC files of precipitation and temperature of a CMIP6 model
; e.g. ts_amon_CESM2_historical_r1i1p1_185001-201412.nc
;      pr_amon_CESM2_historical_r1i1p1_185001-201412.nc
;
;OUT: Thornthwaite's climate types (Dims. 180x360)
;     netcdf file (x1)
;     bil file (x1)
;     hdr file (x1)
;     Plot (x1)
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
;reading data
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

pr=pr*2592000 ;to mm/month

tsClm = clmMonTLL( ts )  ;monthly climatology
tsClm@_FillValue=-9999
tsClm= tsClm-273.15
prClm = clmMonTLL(pr)    ;monthly climatology
prClm@_FillValue=-9999


;coordinate variables and weighting
;---------------------------------------------------------------------
nlon=  360
nlat = 180
gwt  = latGauWgt(nlat, "lat", "gaussian weights", "")

;---------------------------------------------------------------------
;Potential evapotranspiration
;---------------------------------------------------------------------
MEP=new((/12,nlat,nlon/),"float") ;Monthly Potential evapotranspiration
MEP@_FillValue=-9999
P=new((/nlat,nlon/),"float") ;annual precipitation
P=dim_sum_n(prClm,0) ;annual Precipitation
P@_FillValue=-9999
EP=dim_sum_n(tsClm,0) ;annual PE

do ilat=0,nlat-1
  do jlon=0,nlon-1
    MEP(:,ilat,jlon) = thornthwaite(tsClm(:,ilat,jlon), gwt(ilat), False, 0)
  end do
end do

EP=dim_sum_n(MEP,0) ;annual PE
EP!0="lat"
EP!1="lon"
EP&lat=ts&lat
EP&lon=ts&lon
EP@_FillValue=-9999
EP@name="EP"
EP@long_name="Potential evapotranspiration"
EP@comment="Based on Thornthwaite (1948)"
EP@units="mm/yr"

;---------------------------------------------------------------------
;Moisture index (from -1 to 1)
;---------------------------------------------------------------------
MI=new((/nlat,nlon/),"float")
MI@_FillValue=-9999
do ilat=0,nlat-1
  do jlon=0,nlon-1
    if (ismissing(P(ilat,jlon))) then
      MI(ilat,jlon)=-9999
    else if( P(ilat,jlon) .lt. EP(ilat,jlon)) then
      MI(ilat,jlon)=(P(ilat,jlon)/EP(ilat,jlon)-1)
    else if (P(ilat,jlon) .ge. EP(ilat,jlon)) then
      MI(ilat,jlon)=(1-EP(ilat,jlon)/P(ilat,jlon))
    else
      MI(ilat,jlon)=0.0
    end if
    end if
    end if
  end do
end do

MI!0="lat"
MI!1="lon"
MI&lat=ts&lat
MI&lon=ts&lon
MI@name="MoistIndex"
MI@long_name="Moisture Index"
MI@comment="Based on Willmott & Feddema (2005)"
MI@units="dimensionless"

;---------------------------------------------------------------------
;classification
;---------------------------------------------------------------------

;Moisture Index
;---------------------------------------------------------------------
; Saturated 0.66 - 1.00 (1)
; Wet 0.33-0.66         (2)
; Moist 0.00-0.33       (3)
; Dry -0.33 0.00        (4)
; Semiarid -0.33 -0.66  (5)
; Arid -0.66 -1.00      (6)


;thermal (annual data)
;---------------------------------------------------------------------
; Torrid >1500          (1)
; Hot 1200-1500         (2)
; Warm 900-1200         (3)
; Cool 600-900          (4)
; Cold 300-600          (5)
; Frost 0-300           (6)


;---------------------------------------------------------------------
; Create array for Thornthwaite's
;---------------------------------------------------------------------
ThorCl=new((/nlat,nlon/),"integer")
ThorCl!0="lat"
ThorCl!1="lon"
ThorCl&lat=ts&lat
ThorCl&lon=ts&lon
ThorCl@_FillValue=-9999
ThorCl@name="Thornthwaite"
ThorCl@long_name="Revised Thornthwaite"
ThorCl@comment= "Code 1: Torrid-Saturated." +cr+  \
                "Code 2: Hot-Saturated." +cr+  \
                "Code 3: Warm-Saturated." +cr+  \
                "Code 4: Cool-Saturated." +cr+  \
                "Code 5: Cold-Saturated." +cr+  \
                "Code 6: Frigid-Saturated." +cr+  \
                "Code 7: Torrid-Wet." +cr+  \
                "Code 8: Hot-Wet."+cr+  \
                "Code 9: Warm-Wet." +cr+  \
                "Code 10: Cool-Wet." +cr+  \
                "Code 11: Cold-Wet." +cr+  \
                "Code 12: Frigid-Wet." +cr+  \
                "Code 13: Torrid-Moist." +cr+ \
                "Code 14: Hot-Moist." +cr+  \
                "Code 15: Warm-Moist." +cr+  \
                "Code 16: Cool-Moist." +cr+  \
                "Code 17: Cold-Moist." +cr+  \
                "Code 18: Frigid-Moist." +cr+  \
                "Code 19: Torrid-Dry." +cr+  \
                "Code 20: Hot-Dry." +cr+  \
                "Code 21: Warm-Dry." +cr+  \
                "Code 22: Cool-Dry." +cr+  \
                "Code 23: Cold-Dry." +cr+  \
                "Code 24: Frigid-Dry." +cr+  \
                "Code 25: Torrid-Semiarid." +cr+  \
                "Code 26: Hot-Semiarid." +cr+  \
                "Code 27: Warm-Semiarid." +cr+  \
                "Code 28: Cool-Semiarid." +cr+  \
                "Code 29: Cold-Semiarid." +cr+  \
                "Code 30: Frigid-Semiarid." +cr+  \
                "Code 31: Torrid-Arid." +cr+  \
                "Code 32: Hot-Arid." +cr+  \
                "Code 33: Warm-Arid." +cr+  \
                "Code 34: Cool-Arid." +cr+  \
                "Code 35: Cold-Arid." +cr+  \
                "Code 36: Frigid-Arid."

;---------------------------------------------------------------------
; Classifiying and coding
;---------------------------------------------------------------------
ThorCl=where(MI .ge. 0.66 .and. EP .ge. 1500, 1,ThorCl) ;Torrid and saturated
ThorCl=where(MI .lt. 0.66 .and. MI .ge.0.33 .and. EP .ge. 1500, 7,ThorCl) ;Torrid and wet
ThorCl=where(MI .lt. 0.33 .and. MI .ge.0.00 .and. EP .ge. 1500, 13,ThorCl) ;Torrid and moist
ThorCl=where(MI .lt. 0.00 .and. MI .ge.-0.33 .and. EP .ge. 1500, 19,ThorCl) ;Torrid and dry
ThorCl=where(MI .lt. -0.33 .and. MI .ge.-0.66 .and. EP .ge. 1500, 25,ThorCl) ;Torrid and semiarid
ThorCl=where(MI .lt. -0.66 .and. EP .ge. 1500, 31,ThorCl) ;Torrid and arid

ThorCl=where(MI .ge. 0.66 .and. EP .ge. 1200 .and. EP .lt. 1500, 2,ThorCl) ;Hot and saturated
ThorCl=where(MI .lt. 0.66 .and. MI .ge.0.33 .and. EP .ge. 1200 .and. EP .lt. 1500, 8,ThorCl) ;Hot and wet
ThorCl=where(MI .lt. 0.33 .and. MI .ge.0.00 .and. EP .ge. 1200 .and. EP .lt. 1500, 14,ThorCl) ;Hot and moist
ThorCl=where(MI .lt. 0.00 .and. MI .ge.-0.33 .and. EP .ge. 1200 .and. EP .lt. 1500, 20,ThorCl) ;Hot and dry
ThorCl=where(MI .lt. -0.33 .and. MI .ge.-0.66 .and. EP .ge. 1200 .and. EP .lt. 1500, 26,ThorCl) ;Hot and semiarid
ThorCl=where(MI .lt. -0.66 .and. EP .ge. 1200 .and. EP .lt. 1500, 32,ThorCl) ;Hot and arid

ThorCl=where(MI .ge. 0.66 .and. EP .ge. 900 .and. EP .lt. 1200, 3,ThorCl) ;Warm and saturated
ThorCl=where(MI .lt. 0.66 .and. MI .ge.0.33 .and. EP .ge. 900 .and. EP .lt. 1200, 9,ThorCl) ;Warm and wet
ThorCl=where(MI .lt. 0.33 .and. MI .ge.0.00 .and. EP .ge. 900 .and. EP .lt. 1200, 15,ThorCl) ;Warm and moist
ThorCl=where(MI .lt. 0.00 .and. MI .ge.-0.33 .and. EP .ge. 900 .and. EP .lt. 1200, 21,ThorCl) ;Warm and dry
ThorCl=where(MI .lt. -0.33 .and. MI .ge.-0.66 .and. EP .ge. 900 .and. EP .lt. 1200, 27,ThorCl) ;Warm and semiarid
ThorCl=where(MI .lt. -0.66 .and. EP .ge. 900 .and. EP .lt. 1200, 33,ThorCl) ;Warm and arid

ThorCl=where(MI .ge. 0.66 .and. EP .ge. 600 .and. EP .lt. 900, 4,ThorCl) ;Cool and saturated
ThorCl=where(MI .lt. 0.66 .and. MI .ge.0.33 .and. EP .ge. 600 .and. EP .lt. 900, 10,ThorCl) ;Cool and wet
ThorCl=where(MI .lt. 0.33 .and. MI .ge.0.00 .and. EP .ge. 600 .and. EP .lt. 900, 16,ThorCl) ;Cool and moist
ThorCl=where(MI .lt. 0.00 .and. MI .ge.-0.33 .and. EP .ge. 600 .and. EP .lt. 900, 22,ThorCl) ;Cool and dry
ThorCl=where(MI .lt. -0.33 .and. MI .ge.-0.66 .and. EP .ge. 600 .and. EP .lt. 900, 28,ThorCl) ;Cool and semiarid
ThorCl=where(MI .lt. -0.66 .and. EP .ge. 600 .and. EP .lt. 900, 34,ThorCl) ;Cool and arid

ThorCl=where(MI .ge. 0.66 .and. EP .ge. 300 .and. EP .lt. 600, 5,ThorCl) ;Cold and saturated
ThorCl=where(MI .lt. 0.66 .and. MI .ge.0.33 .and. EP .ge. 300 .and. EP .lt. 600, 11,ThorCl) ;Cold and wet
ThorCl=where(MI .lt. 0.33 .and. MI .ge.0.00 .and. EP .ge. 300 .and. EP .lt. 600, 17,ThorCl) ;Cold and moist
ThorCl=where(MI .lt. 0.00 .and. MI .ge.-0.33 .and. EP .ge. 300 .and. EP .lt. 600, 23,ThorCl) ;Cold and dry
ThorCl=where(MI .lt. -0.33 .and. MI .ge.-0.66 .and. EP .ge. 300 .and. EP .lt. 600, 29,ThorCl) ;Cold and semiarid
ThorCl=where(MI .lt. -0.66 .and. EP .ge. 300 .and. EP .lt. 600, 35,ThorCl) ;Cold and arid

ThorCl=where(MI .ge. 0.66 .and.  EP .lt. 300, 6,ThorCl) ;Frigid and saturated
ThorCl=where(MI .lt. 0.66 .and. MI .ge.0.33 .and. EP .lt. 300, 12,ThorCl) ;Frigid and wet
ThorCl=where(MI .lt. 0.33 .and. MI .ge.0.00 .and. EP .lt. 300, 18,ThorCl) ;Frigid and moist
ThorCl=where(MI .lt. 0.00 .and. MI .ge.-0.33 .and. EP .lt. 300, 24,ThorCl) ;Frigid and dry
ThorCl=where(MI .lt. -0.33 .and. MI .ge.-0.66 .and. EP .lt. 300, 30,ThorCl) ;Frigid and semiarid
ThorCl=where(MI .lt. -0.66 .and. EP .lt. 300, 36,ThorCl) ;Frigid and arid


;---------------------------------------------------------------------
;Read CRU data to mask OCN in plots
;---------------------------------------------------------------------
inmaskpath=addfile("ts_amon_Z-CRUTS_historical_r1i1p1_190101-201412.nc","r")
inmask=inmaskpath->ts(948,:,:)
inmask@_FillValue=-9999
ThorCl=where(ismissing(inmask),-9999,ThorCl)


;---------------------------------------------------------------------
;plotting parameters
;---------------------------------------------------------------------

wks = gsn_open_wks("pdf", modelname)

  ;; Plot resources
  res = True
  res@cnFillOn = True
  res@cnLinesOn = False
  ;res@cnFillMode = "RasterFill"
  res@cnFillMode = "CellFill"
  ;res@cnRasterSmoothingOn = True
  res@gsnAddCyclic = True
  res@lbLabelBarOn   = False
  res@cnLineLabelsOn=False
  res@mpProjection="Robinson"
  res@mpGridAndLimbOn       =  False
  res@mpPerimOn =   False
  res@mpLimitMode ="LatLon"
  res@mpMinLatF   =-60.0
  res@mpMaxLatF   = 90.0
  long_name = (/"Thornthwaite ", modelname, case/)
  ThorCl@long_name = str_concat(long_name)
  res@mpDataBaseVersion = "MediumRes"
  res@mpOutlineOn = True
  ;res@pmTickMarkDisplayMode = "Always"
  colors=(/"gray20","gray30","gray40","gray60","gray70","gray90", \
        "steelblue4","steelblue3","steelblue2","steelblue1","lightskyblue","lightskyblue1",\
        "darkolivegreen","olivedrab4","olivedrab3","olivedrab2","olivedrab1","darkolivegreen1", \
        "lightgoldenrod4","lightgoldenrod3","lightgoldenrod2","lightgoldenrod","lightgoldenrod1","lightgoldenrodyellow", \
        "sienna4","sienna3","sienna2","sienna1","tan1","peachpuff1", \
        "red4","red3","red2","red1","indianred1","lightpink"/)
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnFillPalette = colors
  labels = (/"Ts","Hs","Ws","Cs","CDs","Fs","Tw","Hw","Ww","Cw","CDw","Fw","Tm","Hm","Wm","Cm","CDm","Fm", \
            "Td","Hd","Wd","Cd","CDd","Fd","Tsa","Hsa","Wsa","Csa","CDsa","Fsa","Ta","Ha","Wa","Ca","CDa","Fa"/)
  res@lbLabelStrings = labels
  res@cnLevels = (/1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5, \
                  21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5/)
  res@lbLabelAlignment = "BoxCenters"
  res@lbOrientation = "vertical"
  res@lbAutoManage = True

 map = gsn_csm_contour_map(wks,ThorCl,res);


 ;---------------------------------------------------------------------
 ;Flip longitude to -180 +180
 ;---------------------------------------------------------------------
 ThorCl=lonFlip(ThorCl)
 delete(ThorCl@lonFlip)
 EP=lonFlip(EP)
 delete(EP@lonFlip)
 MI=lonFlip(MI)
 delete(MI@lonFlip)

 ;---------------------------------------------------------------------
 ;write to netcdf file
 ;---------------------------------------------------------------------
 qq=addfile("thornthwaite_class_"+case4+"_"+modelname+".nc","c")
 ThorCl@long_name="Revised Thornthwaite"
 qq->MI=MI
 qq->EP=EP
 qq->thorn=ThorCl

 ;---------------------------------------------------------------------
 ;write to bil file (and .hdr)
 ;---------------------------------------------------------------------
 cbinwrite("thornthwaite_class_"+case4+"_"+modelname+".bil", ThorCl(::-1,:))  ;reverse lat

 ;adding data to the hdr file
 strings = (/"BYTEORDER I","LAYOUT BIL","NBANDS 1","NROWS 180","NCOLS 360", \
 "ULXMAP -180.0","ULYMAP 89.5","NBITS 32","NODATA -9999","PROJECTION GEOGRAPHIC", \
 "DATUM WGS84","SPHEROID EPSG4326"/)
 alist   = [/strings/]
 write_table("thornthwaite_class_"+case4+"_"+modelname+".hdr", "w", alist, "%s") ;write .hdr
