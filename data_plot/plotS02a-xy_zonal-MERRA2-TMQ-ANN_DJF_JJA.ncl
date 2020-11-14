load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "/home/yihsuan/Downloads/yihsuan//script/ncl//ignorantNCL.ncl"
load "./ignorantNCL.ncl"
;load "/home/yihsuan/Downloads/yihsuan//script/ncl//ignorantNCL-colormaps.ncl"
;load "zz-read_parameters.ncl"

begin

;------------------
; general settings
;------------------

  filename_ANN = "./data_plotS02a-MERRA2_zonal_TCWV_2006_2015_Annual.txt"
  filename_DJF = "./data_plotS02a-MERRA2_zonal_TCWV_2006_2015_DJF.txt"
  filename_JJA = "./data_plotS02a-MERRA2_zonal_TCWV_2006_2015_JJA.txt"

  wks_name = "figS02a-xy_zonal-MERRA2-TMQ-ANN_DJF_JJA"
  wks_type = "eps"  ; or x11,pdf,ncgm

  nlat = 361
  nvar = 2

  data_ANN = asciiread(filename_ANN, (/nlat, nvar/), "float")
  TMQ_ANN = data_ANN(:,1)
  TMQ_ANN = multiply(TMQ_ANN,10.) ; cm to mm

  data_DJF = asciiread(filename_DJF, (/nlat, nvar/), "float")
  TMQ_DJF = data_DJF(:,1)
  TMQ_DJF = multiply(TMQ_DJF,10.) ; cm to mm

  data_JJA = asciiread(filename_JJA, (/nlat, nvar/), "float")
  TMQ_JJA = data_JJA(:,1)
  TMQ_JJA = multiply(TMQ_JJA,10.) ; cm to mm

  lat = data_ANN(:,0)
  lat@units = "degrees_north"

  TMQ_all = (/TMQ_ANN,TMQ_DJF,TMQ_JJA/)

;--------
; plot
;--------
  opt = True
    opt@legend_mean = False

  labels = (/"ANN","DJF","JJA"/)

  res = True
  res@tiMainString = "MERRA-2 10-year mean (2006-2015)"
  res@tiYAxisString = "Total column water vapor (mm)"
  res@tiXAxisString = "Latitude"

  res@xyMonoLineThickness = True
  res@xyLineThicknessF = 3.

  wks = gsn_open_wks(wks_type, wks_name)

  plot = ignorant_xy(wks, lat, TMQ_all, labels, res, opt)

  opt1 = True
  trim_fig(wks_type, wks_name, opt1)

end
