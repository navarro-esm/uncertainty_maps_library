;-------------------------------------------------------------------------------
; Köppen climate types.
;-------------------------------------------------------------------------------
; Based on Lohmann, U., Sausen, R., Bengtsson, L., Cubasch, U., Perlwitz, J.
; and Roeckner, E. (1993) The Köppen climate classification as a diagnostic tool
; for general circulation models. Climate Research, 3, 177–193.
;
; To be used after classifier.out. This is not the classifier code.
; Author: Andres Navarro
; e-mail:andres.navarro@uclm.es
;-------------------------------------------------------------------------------
;IN: intermediate bil files (Fortran direct access write)
;    created by the classifier code (classifier.out)
;    e.g. koppen_class_ACCESS-CM2.bil
;
;
;OUT: Köppen climate classification (Dims. 180x360)
;     netcdf file (x1)
;     bil file (x1) ->standard format (180x360)
;     hdr file (x1)
;     Plot (x1)
;-------------------------------------------------------------------------------


debug=False    ;True or False
modelname=getenv("MODELNAME")
case=getenv("CASE")
print(modelname)
case2 = str_split(case, "-")  ;strip after -
case3 =str_left_strip(case2)  ;remove white spaces
case4=str_lower(case3)      ;lowercase
cr = tochar(10)


;coordinate variables and weighting
;---------------------------------------------------------------------
nlon=  360
nlat = 180
gwt  = latGauWgt(nlat, "lat", "gaussian weights", "")


;---------------------------------------------------------------------
; Reading data from intermediate bil files
;---------------------------------------------------------------------

if (case4 .eq. "historical") then
  if(modelname .eq. "Z-CRUTS") then
    filepath="/case_historical/koppen_source/"
    data = fbindirread(filepath+"koppen_class_"+modelname+".bil", 0, (/nlat + nlon + (nlat * nlon)/), "float")
    lat = data(0:(nlat-1))
    lon = data(nlat:nlat + (nlon-1))
    lat@units = "degrees_north"
    lon@units = "degrees_east"
    cls = reshape(data((nlat + nlon):), (/nlat,nlon/))
    cls!0 = "lat"
    cls!1 = "lon"
    cls@_FillValue = 99
    cls&lat = lat
    cls&lon = lon
    cls@_FillValue=-9999
  else
    ;model
    filepath="/case_historical/koppen_source/"
    data = fbindirread(filepath+"koppen_class_"+modelname+".bil", 0, (/nlat + nlon + (nlat * nlon)/), "float")
    lat = data(0:(nlat-1))
    lon = data(nlat:nlat + (nlon-1))
    lat@units = "degrees_north"
    lon@units = "degrees_east"
    cls = reshape(data((nlat + nlon):), (/nlat,nlon/))
    cls!0 = "lat"
    cls!1 = "lon"
    cls@_FillValue = 99
    cls&lat = lat
    cls&lon = lon

    ;cru mask
    data2 = fbindirread(filepath+"koppen_class_Z-CRUTS.bil", 0, (/nlat + nlon + (nlat * nlon)/), "float")
    refdata = reshape(data2((nlat + nlon):), (/nlat,nlon/))
    refdata!0 = "lat"
    refdata!1 = "lon"
    refdata@_FillValue = 99
    refdata&lat = lat
    refdata&lon = lon

    cls=where(refdata .eq. 99,99,cls)   ;masking data
    cls@_FillValue=-9999                ;new fillvalue
    refdata@_FillValue=-9999
  end if
end if

if (case4 .eq. "ssp585") then
  ;model
  filepath="/case_ssp585/koppen_source/"
  data = fbindirread(filepath+"koppen_class_"+modelname+".bil", 0, (/nlat + nlon + (nlat * nlon)/), "float")
  lat = data(0:(nlat-1))
  lon = data(nlat:nlat + (nlon-1))
  lat@units = "degrees_north"
  lon@units = "degrees_east"
  cls = reshape(data((nlat + nlon):), (/nlat,nlon/))
  cls!0 = "lat"
  cls!1 = "lon"
  cls@_FillValue = 99
  cls&lat = lat
  cls&lon = lon

  ;cru mask
  filepath2="/case_ssp585/koppen_source/"
  data2 = fbindirread(filepath2+"koppen_class_Z-CRUTS.bil", 0, (/nlat + nlon + (nlat * nlon)/), "float")
  refdata = reshape(data2((nlat + nlon):), (/nlat,nlon/))
  refdata!0 = "lat"
  refdata!1 = "lon"
  refdata@_FillValue = 99
  refdata&lat = lat
  refdata&lon = lon

  cls=where(refdata .eq. 99,99,cls)   ;masking data
  cls@_FillValue=-9999                ;new fillvalue
  refdata@_FillValue=-9999
end if

if (case4 .eq. "ssp126") then
  ;model
  filepath="/case_ssp126/koppen_source/"
  data = fbindirread(filepath+"koppen_class_"+modelname+".bil", 0, (/nlat + nlon + (nlat * nlon)/), "float")
  lat = data(0:(nlat-1))
  lon = data(nlat:nlat + (nlon-1))
  lat@units = "degrees_north"
  lon@units = "degrees_east"
  cls = reshape(data((nlat + nlon):), (/nlat,nlon/))
  cls!0 = "lat"
  cls!1 = "lon"
  cls@_FillValue = 99
  cls&lat = lat
  cls&lon = lon

  ;cru mask
  filepath2="/case_ssp585/koppen_source/"
  data2 = fbindirread(filepath2+"koppen_class_Z-CRUTS.bil", 0, (/nlat + nlon + (nlat * nlon)/), "float")
  refdata = reshape(data2((nlat + nlon):), (/nlat,nlon/))
  refdata!0 = "lat"
  refdata!1 = "lon"
  refdata@_FillValue = 99
  refdata&lat = lat
  refdata&lon = lon

  cls=where(refdata .eq. 99,99,cls)   ;masking data
  cls@_FillValue=-9999                ;new fillvalue
  refdata@_FillValue=-9999
end if

if (case4 .eq. "ssp245") then
  ;model
  filepath="/case_ssp245/koppen_source/"
  data = fbindirread(filepath+"koppen_class_"+modelname+".bil", 0, (/nlat + nlon + (nlat * nlon)/), "float")
  lat = data(0:(nlat-1))
  lon = data(nlat:nlat + (nlon-1))
  lat@units = "degrees_north"
  lon@units = "degrees_east"
  cls = reshape(data((nlat + nlon):), (/nlat,nlon/))
  cls!0 = "lat"
  cls!1 = "lon"
  cls@_FillValue = 99
  cls&lat = lat
  cls&lon = lon

  ;cru mask
  filepath2="/case_ssp585/koppen_source/"
  data2 = fbindirread(filepath2+"koppen_class_Z-CRUTS.bil", 0, (/nlat + nlon + (nlat * nlon)/), "float")
  refdata = reshape(data2((nlat + nlon):), (/nlat,nlon/))
  refdata!0 = "lat"
  refdata!1 = "lon"
  refdata@_FillValue = 99
  refdata&lat = lat
  refdata&lon = lon

  cls=where(refdata .eq. 99,99,cls)   ;masking data
  cls@_FillValue=-9999                ;new fillvalue
  refdata@_FillValue=-9999
end if

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
  res@gsnAddCyclic = True
  res@lbLabelBarOn   = False
  res@cnLineLabelsOn=False
  res@mpProjection="Robinson"
  res@mpGridAndLimbOn       =  False
  res@mpPerimOn =   False
  res@mpLimitMode ="LatLon"
  res@mpMinLatF   =-60.0
  res@mpMaxLatF   = 90.0

  ;; Map
  long_name = (/"koppen ", modelname," ", case/)
  cls@long_name = str_concat(long_name)
  res@mpDataBaseVersion = "MediumRes"
  res@mpOutlineOn = True
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnFillPalette = "navarro_200c"      ;this is a custom palette
  labels = (/"Af","Aw", "BS","BW", "Cs","Cw","Cf", "Dw","Df", "ET","EF"/)
  res@lbLabelStrings = labels
  res@cnLevels = (/0.5,1.5, 2.5,3.5, 4.5,5.5,6.5, 7.5,8.5, 9.5/)
  res@cnFillColors = (/73,71, 219,105, 237,247,227, 205,207, 193,192/)
  res@lbLabelAlignment = "BoxCenters"
  res@lbOrientation = "vertical"
  res@lbAutoManage = True

  map = gsn_csm_contour_map(wks,cls,res); plot koppen


;----------------------------------------------------------------------------------
; Common attributes for koppen climate classification.
; Convert cls to integer. Conversion does not preserve attributes.
;----------------------------------------------------------------------------------
koppen=toint(cls)
koppen!0="lat"
koppen!1="lon"
koppen&lat=cls&lat
koppen&lon=cls&lon
koppen@long_name="Köppen climate classification"
koppen@name="Köppen"
koppen@comment= "Code 0: Af." +cr+  \
                "Code 1: Aw." +cr+  \
                "Code 2: BS." +cr+  \
                "Code 3: BW." +cr+  \
                "Code 4: Cs." +cr+  \
                "Code 5: Cw." +cr+  \
                "Code 6: Cf." +cr+  \
                "Code 7: Dw."+cr+  \
                "Code 8: Df." +cr+  \
                "Code 9: ET." +cr+  \
                "Code 10: EF."


;----------------------------------------------------------------------------------
;Flip longitude to -180 +180 (reshape coordinate variables)
;----------------------------------------------------------------------------------
koppen=lonFlip(koppen)
delete(koppen@lonFlip)

;----------------------------------------------------------------------------------
;write to netcdf
;----------------------------------------------------------------------------------
zz=addfile("koppen_class_"+case4+"_"+modelname+".nc","c")
;zz->koppen=cls
zz->koppen=koppen
zz->lat=koppen&lat
zz->lon=koppen&lon

;---------------------------------------------------------------------
;write to bil file (and .hdr)
;---------------------------------------------------------------------
;cbinwrite("koppen_class_"+case4+"_"+modelname+".bil", cls(::-1,:))  ;reverse lat
cbinwrite("koppen_class_"+case4+"_"+modelname+".bil", koppen(::-1,:))  ;reverse lat

;adding data to the hdr file
strings = (/"BYTEORDER I","LAYOUT BIL","NBANDS 1","NROWS 180","NCOLS 360", \
"ULXMAP -180.0","ULYMAP 89.5","NBITS 32","NODATA -9999","PROJECTION GEOGRAPHIC", \
"DATUM WGS84","SPHEROID EPSG4326"/)
alist   = [/strings/]
write_table("koppen_class_"+case4+"_"+modelname+".hdr", "w", alist, "%s") ;write .hdr
