;--------------------------------------------------------------------------------
; Uncertainty analysis for top-10 CMIP6 models.
; Climate classifications: Holdridge, Köppen, Thornthwaite & Whittaker.
; Author: A.Navarro.
; e-mail: andres.navarro@uclm.es
;--------------------------------------------------------------------------------
;IN:  Bil files of top-10 models (10x180x360).-> Read from filenames: Z-per80.txt
;     Bil file reference data (1x180x360)
;       Present: CRU
;       Future: Top-10 ensemble mean.
;
;OUT: Consensus maps (Dims. 180x360)
;     netcdf file (x1)
;     bil files (x2): confidence & modvar
;     Plot (x1)
;--------------------------------------------------------------------------------


;--------------------------------------------------------------------------------
;Parameters
;--------------------------------------------------------------------------------
case=getenv("CASE")         ; historical, ssp126, ssp245, ssp585
tclass=getenv("CLASS")      ; holdridge, köppen,thornthwaite, whittaker
;print(case)
;print(tclass)
oplot=1                    ; 0= do no overplot inter-model variability ; 1= overplot
cr = tochar(10)

;dimensions
nlat=180
nlon=360

;coordinates
lat = fspan(89.5,-89.5,nlat)
lon = fspan(-180.0,179.0,nlon)
;lat = fspan(-90,90,nlat)
;lon = fspan(0,359.5,nlon)
lat!0="lat"
lon!0="lon"
lat&lat=lat
lon&lon=lon
lat@units = "degrees_north"
lat@standard_name = "latitude" ;
lat@long_name = "latitude" ;
lat@axis = "Y" ;
lon@units = "degrees_east"
lon@standard_name = "longitude" ;
lon@long_name = "longitude" ;
lon@axis = "X" ;



;--------------------------------------------------------------------------------
; Reading reference data
; CRU-> for historical, TOP10-ENSEMBLE for SSP585, SSP245 and SSP126)
;--------------------------------------------------------------------------------
if (case .eq. "ssp245") then
  if (tclass .eq. "holdridge") then
    crupath=("/case_ssp245/holdridge/holdridge_class_ssp245_Z-ENS-HOLDRIDGE.bil")
  else if (tclass .eq. "thornthwaite")
    crupath=("/case_ssp245/thornthwaite/thornthwaite_class_ssp245_Z-ENS-THORNTHWAITE.bil")
  else if (tclass .eq. "whittaker")
    crupath=("/case_ssp245/whittaker/whittaker_class_ssp245_Z-ENS-WHITTAKER.bil")
  else
    crupath=("/case_ssp245/koppen/koppen_class_ssp245_Z-ENS-KOPPEN.bil")
  end if
  end if
  end if
else if (case .eq. "ssp585")
  if (tclass .eq. "holdridge") then
    crupath=("/case_ssp585/holdridge/holdridge_class_ssp585_Z-ENS-HOLDRIDGE.bil")
  else if (tclass .eq. "thornthwaite")
    crupath=("/case_ssp585/thornthwaite/thornthwaite_class_ssp585_Z-ENS-THORNTHWAITE.bil")
  else if (tclass .eq. "whittaker")
    crupath=("/case_ssp585/whittaker/whittaker_class_ssp585_Z-ENS-WHITTAKER.bil")
  else
    crupath=("/case_ssp585/koppen/koppen_class_ssp585_Z-ENS-KOPPEN.bil")
  end if
  end if
  end if
else if (case .eq. "ssp126")
  if (tclass .eq. "holdridge") then
    crupath=("/case_ssp126/holdridge/holdridge_class_ssp126_Z-ENS-HOLDRIDGE.bil")
  else if (tclass .eq. "thornthwaite")
    crupath=("/case_ssp126/thornthwaite/thornthwaite_class_ssp126_Z-ENS-THORNTHWAITE.bil")
  else if (tclass .eq. "whittaker")
    crupath=("/case_ssp126/whittaker/whittaker_class_ssp126_Z-ENS-WHITTAKER.bil")
  else
    crupath=("/case_ssp126/koppen/koppen_class_ssp126_Z-ENS-KOPPEN.bil")
  end if
  end if
  end if
else    ;for historical
  if (tclass .eq. "holdridge") then
    crupath=("/case_historical/holdridge/holdridge_class_historical_Z-CRUTS.bil")
  else if (tclass .eq. "thornthwaite")
    crupath=("/case_historical/thornthwaite/thornthwaite_class_historical_Z-CRUTS.bil")
  else if (tclass .eq. "whittaker")
    crupath=("/case_historical/whittaker/whittaker_class_historical_Z-CRUTS.bil")
  else
    crupath=("/case_historical/koppen/koppen_class_historical_Z-CRUTS.bil")
  end if
  end if
  end if
end if
end if
end if

;--------------------------------------------------------------------------------
;Reference data (Historical= CRU ; FUTURE= Top-10 Ensemble)
;--------------------------------------------------------------------------------
cruin=cbinread(crupath, (/nlat , nlon /), "integer")
cruin!0 = "lat"
cruin!1 = "lon"
cruin@_FillValue = -9999


;--------------------------------------------------------------------------------
;BEST MODELS (Top- 10 models)
;--------------------------------------------------------------------------------

if (case .eq. "ssp585") then
  if (tclass .eq. "holdridge") then
    modpath="/case_ssp585/holdridge/"
  else if (tclass .eq. "thornthwaite")
    modpath="/case_ssp585/thornthwaite/"
  else if (tclass .eq. "whittaker")
    modpath="/case_ssp585/whittaker/"
  else
    modpath="/case_historical/koppen/"
  end if
  end if
  end if
else if (case .eq. "ssp245")
  if (tclass .eq. "holdridge") then
    modpath="/case_ssp245/holdridge/"
  else if (tclass .eq. "thornthwaite")
    modpath="/case_ssp245/thornthwaite/"
  else if (tclass .eq. "whittaker")
    modpath="/case_ssp245/whittaker/"
  else
    modpath="/case_ssp245/koppen/"
  end if
  end if
  end if
else if (case .eq. "ssp126")
  if (tclass .eq. "holdridge") then
    modpath="/case_ssp126/holdridge/"
  else if (tclass .eq. "thornthwaite")
    modpath="/case_ssp126/thornthwaite/"
  else if (tclass .eq. "whittaker")
    modpath="/case_ssp126/whittaker/"
  else
    modpath="/case_ssp126/koppen/"
  end if
  end if
  end if
else
  if (tclass .eq. "holdridge") then
    modpath="/case_historical/holdridge/"
  else if (tclass .eq. "thornthwaite")
    modpath="/case_historical/thornthwaite/"
  else if (tclass .eq. "whittaker")
    modpath="/case_historical/whittaker/"
  else
    modpath="/case_historical/koppen/"
  end if
  end if
  end if
end if
end if
end if


;--------------------------------------------------------------------------------
;Reading list top-10 models for each CCS
;--------------------------------------------------------------------------------
per80in=asciiread(modpath+"Z-per80.txt",-1,"string")
nper80=dimsizes(per80in)            ;nmodels

;--------------------------------------------------------------------------------
;create new 3D array to store information
;--------------------------------------------------------------------------------

ens80=new((/nper80,180,360/),"integer")    ; Climate classification for all models in one array
do nf=0,nper80-1
  fi=cbinread(modpath+per80in(nf), (/nlat , nlon /), "integer")
  fi!0 = "lat"
  fi!1 = "lon"
  fi@_FillValue = -9999
  fi&lat = lat
  fi&lon = lon
  fi=where(ismissing(cruin),-9999,fi)
  ens80(nf,:,:)=fi
end do

;--------------------------------------------------------------------------------
;check (dis)agreement between CRU and Models
;--------------------------------------------------------------------------------
agreens80=new((/nper80,180,360/),"integer")
agreens80!0 = "nmodels"
agreens80!1 = "lat"
agreens80!2 = "lon"
agreens80&lat = lat
agreens80&lon = lon
agreens80@_FillValue = -9999

do nf=0,nper80-1
  agreens80(nf,:,:)=where(cruin .eq. ens80(nf,:,:),1,-1)      ;1 agreement -1 disagreement, -9999 _FillValue
end do

;--------------------------------------------------------------------------------
;ratio of agreement (0 to 1)
;--------------------------------------------------------------------------------
peragreens80=new((/180,360/),"float")
peragreens80!0 = "lat"
peragreens80!1 = "lon"
peragreens80&lat = lat
peragreens80&lon = lon
peragreens80@_FillValue = -9999
peragreens80@units="%"
peragreens80@name="Confidence"
peragreens80@long_name="Confidence levels"

if (case .eq. "historical") then
  peragreens80@comment= "Percent agreement CRU & Top-10 models"
else
  peragreens80@comment= "Percent agreement T10 Ensemble & Top-10 models"
end if


do ilat=0,nlat-1
  do jlon=0,nlon-1
    peragreens80(ilat,jlon)=num(agreens80(:,ilat,jlon).eq.1 )/tofloat(nper80)
    ;peragreens80(ilat,jlon)=num(agreens80(:,ilat,jlon).eq.1 )/nper80
  end do
end do
peragreens80=where(ismissing(cruin),-9999,peragreens80)   ;set missing vals

;--------------------------------------------------------------------------------
;COMPUTING CONFIDENCE
;--------------------------------------------------------------------------------
; >=80% [dark blue] color code=1
; 79.99-60% [light blue] color code=2
; 59.99-40% [yellow] color code =3
; 39.99-20% [light red] color code=4
; <20% [dark-red] color code=5
;--------------------------------------------------------------------------------

;Creating new array
reg80=new((/180,360/),"integer")
reg80!0 = "lat"
reg80!1 = "lon"
reg80&lat = lat
reg80&lon = lon
reg80@_FillValue = -9999

;Adding description for netcdf file
if (case .eq. "historical") then
  reg80@name="Confidence"
  reg80@long_name="Confidence levels"
  reg80@comment="Confidence level of the top-10 models for historical climate (1980-2014)"+cr+ \
  "Code 1: >=80% confidence."+cr+ \
  "Code 2: 79.99-60% confidence."+cr+ \
  "Code 3: 59.99-40% confidence."+cr+ \
  "Code 4: 39.99-20% confidence."+cr+ \
  "Code 5: <20% confidence."
  reg80name="top10_consensus_historical_"+tclass
else if (case .eq. "ssp126")
  reg80@name="Confidence"
  reg80@long_name="Confidence levels"
  reg80@comment="Confidence level of the top-10 models for SSP1-2.6 (2015-2100)." +cr+ \
  "Code 1: >=80% confidence."+cr+ \
  "Code 2: 79.99-60% confidence."+cr+ \
  "Code 3: 59.99-40% confidence."+cr+ \
  "Code 4: 39.99-20% confidence."+cr+ \
  "Code 5: <20% confidence."
  reg80name="top10_consensus_ssp126_"+tclass
else if (case .eq. "ssp245")
  reg80@name="Confidence"
  reg80@long_name="Confidence levels"
  reg80@comment="Confidence level of the top-10 models for SSP2-4.5 (2015-2100)."+cr+ \
  "Code 1: >=80% confidence."+cr+ \
  "Code 2: 79.99-60% confidence."+cr+ \
  "Code 3: 59.99-40% confidence."+cr+ \
  "Code 4: 39.99-20% confidence."+cr+ \
  "Code 5: <20% confidence."
  reg80name="top10_consensus_ssp245_"+tclass
else
  reg80@name="Confidence"
  reg80@long_name="Confidence levels"
  reg80@comment="Confidence level of the top-10 models for SSP5-8.5 (2015-2100)."+cr+ \
  "Code 1: >=80% confidence."+cr+ \
  "Code 2: 79.99-60% confidence."+cr+ \
  "Code 3: 59.99-40% confidence."+cr+ \
  "Code 4: 39.99-20% confidence."+cr+ \
  "Code 5: <20% confidence."
  reg80name="top10_consensus_ssp585_"+tclass
end if
end if
end if

reg80=where(peragreens80 .ge. 0.8, 1,reg80)                                    ; totally agree [>=80%]
reg80=where(peragreens80 .lt. 0.8 .and. peragreens80 .ge. 0.6 , 2,reg80)       ;  agree [79.99-60%]
reg80=where(peragreens80 .lt. 0.6 .and. peragreens80 .ge. 0.4 , 3,reg80)       ;  partially agree [59.99-40%]
reg80=where(peragreens80 .lt. 0.4 .and. peragreens80 .ge. 0.2 , 4,reg80)       ;  disagree [39.99-20%]
reg80=where(peragreens80 .le. 0.2, 5,reg80)                                     ;totally disagree [<20%]

;--------------------------------------------------------------------------------
;Inter-model variability (modvar).
;--------------------------------------------------------------------------------

climtypes=new((/180,360/),"integer")             ; class agreement
climtypes!0 = "lat"
climtypes!1 = "lon"
climtypes&lat = lat
climtypes&lon = lon
climtypes@_FillValue = -9999
climtypes@name="modvar"
climtypes@long_name= "Count of unique climate zones"
climtypes@comment="Multiple climate types suggest high variation in model predictions."


do i=0,nlat-1
  do j=0,nlon-1
    if (.not. ismissing(ens80(0,i,j))) then
      countclimtypes = get_unique_values(ens80(:,i,j))
      nsizeclimtypes=dimsizes(countclimtypes)
      climtypes(i,j)=nsizeclimtypes
      delete(countclimtypes)
    else
      climtypes(i,j)=-9999     ;FillValue
    end if
 end do
end do

;--------------------------------------------------------------------------------
;Plotting consensus maps
;--------------------------------------------------------------------------------

wks = gsn_open_wks("pdf", reg80name)

;COMMON
res                = True
res@gsnFrame       = False
res@gsnDraw        = False
res@gsnLeftString  = ""
res@gsnRightString = ""
res@gsnTickMarksOn = True
res@lbLabelBarOn   = False
res@cnLineLabelsOn = False
res@mpProjection="Robinson"
res@mpGridAndLimbOn       =  False
res@mpPerimOn =   False
res@mpLimitMode ="LatLon"
res@mpMinLatF   =-60.0
res@mpMaxLatF   = 90.0
res@gsnAddCyclic=True

;MAP 1
res1 = res
res1@cnFillOn = True
res1@cnLinesOn = False
res1@cnFillMode = "RasterFill"
res1@mpDataBaseVersion = "MediumRes"
res1@mpOutlineOn = True
res1@cnLevelSelectionMode = "ExplicitLevels"
res1@cnFillPalette = (/"steelblue3","steelblue1","yellow","salmon","tomato3"/)
labels2 = (/"","","","",""/)
res1@lbLabelStrings = labels2
res1@cnLevels = (/1.5,2.5,3.5,4.5/)

plot = gsn_csm_contour_map(wks,reg80,res1);

if (oplot .eq. 1) then        ;check if we use overplot

  ;MAP 2
  res2                      = res
  res2@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
  if (tclass .eq. "thornthwaite" ) then
      res2@cnLevels             =fspan(0,9,3)
  else
    res2@cnLevels             =fspan(0,5,3)
  end if
  res2@cnLineThicknessF    =0.1
  res2@cnLinesOn    =False
  res2@cnInfoLabelOn = False
  plot2   = gsn_csm_contour(wks,climtypes, res2)

  opt     = True
  opt@cnFillOn=True
  opt@gsnShadeFillType = "pattern"
  opt@gsnShadeMid = 4 ;10 dash 17 dots (between lower and maximum value)
  opt@gsnShadeFillScaleF=0.35 ;density
  opt@gsnShadeFillDotSizeF = 0.00      ; dot size
  opt@cnLineLabelsOn       = False
  if (tclass .eq. "thornthwaite" ) then
    plot2   = gsn_contour_shade(plot2, 5,9, opt)  ; plot values >5 and < 9
  else
    plot2   = gsn_contour_shade(plot2, 3,5, opt)  ; plot values >3 and < 5
  end if
  overlay(plot,plot2)
  draw(plot)
  frame(wks)
else
  draw(plot)
  frame(wks)
end if


;--------------------------------------------------------------------------------
;Percent confidence
;--------------------------------------------------------------------------------
peragreens80=peragreens80*100
peragreens80=decimalPlaces(peragreens80,2,True) ;round
;--------------------------------------------------------------------------------
;Create netcdf file
;--------------------------------------------------------------------------------
gg=addfile(tclass+"_consensus_"+case+".nc","c")
gg->lat=lat
gg->lon=lon
gg->confidence=peragreens80  ;accuracy (confidence). In percentage
gg->modvar=climtypes  ;precision (modvar)

;--------------------------------------------------------------------------------
;Create bil files
;--------------------------------------------------------------------------------
cbinwrite(tclass+"_confidence_"+case+".bil", peragreens80)   ;accuracy (confidence)
cbinwrite(tclass+"_modvar_"+case+".bil",climtypes)    ;precision (modvar)

;adding data to the hdr file
strings1 = (/"BYTEORDER I","LAYOUT BIL","NBANDS 1","NROWS 180","NCOLS 360", \
"ULXMAP -180.0","ULYMAP 89.5","NBITS 32","PIXELTYPE FLOAT","NODATA -9999",  \
"PROJECTION GEOGRAPHIC","DATUM WGS84","SPHEROID EPSG4326"/)
alist1   = [/strings1/]
write_table(tclass+"_confidence_"+case+".hdr", "w", alist1, "%s") ;write .hdr

strings2 = (/"BYTEORDER I","LAYOUT BIL","NBANDS 1","NROWS 180","NCOLS 360", \
"ULXMAP -180.0","ULYMAP 89.5","NBITS 32","NODATA -9999",  \
"PROJECTION GEOGRAPHIC","DATUM WGS84","SPHEROID EPSG4326"/)
alist2   = [/strings2/]
write_table(tclass+"_modvar_"+case+".hdr", "w", alist2, "%s") ;write .hdr
