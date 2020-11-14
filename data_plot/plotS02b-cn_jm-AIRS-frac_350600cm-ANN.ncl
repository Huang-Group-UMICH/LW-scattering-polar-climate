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

;-----------------------------
; open files & read varibles
;-----------------------------

  wks_name = "figS02b-cn_jm-AIRS-frac_350600cm"
  wks_type = "eps"  ; or x11,pdf,ncgm

  file_AIRS = addfile("./data_plotS02b-AIRS2CERES_zonal_flux_RRTMG_LW_2006_2015_monthly.nc", "r")

  airs_flux_all = file_AIRS->airs_flux_all
  airs_flux_all@_FillValue = 9.96921e+36 

  airs_flux_broad = dim_sum_n_Wrap(airs_flux_all,1)
  ;printVarSummary(airs_flux_broad)

  airs_flux_350630 = dim_sum_n_Wrap(airs_flux_all(:,1:2,:),1)
  ;print(airs_flux_350630)

  airs_flux_frac = divide(airs_flux_350630,airs_flux_broad)
  airs_flux_frac = multiply(airs_flux_frac,100.)

  airs_flux_frac@units = "%"
  ;print(airs_flux_350630+" , "+airs_flux_broad+", "+airs_flux_frac)
  ;printVarSummary(airs_flux_frac)

;*** contour data ***
  ;cndata = airs_flux_frac(:,:)
  cndata = airs_flux_frac(:,:)
  cndata&lat@units = "degrees_north"

  printVarSummary(cndata)
  ;printVarSummary(cndata&lat)

;*** xy plot 1 - x_avg.vs.y *** 
; y-axis data  
  dim0  = cndata!0
  y1data = cndata&$dim0$

; x-axis data
  x1 = dim_avg_n_Wrap(cndata,1)
  x1data = (/x1/)

;*** xy plot 2 - x .vs.y_avg *** 
; y-axis data  
  dim1  = cndata!1
  x2data = cndata&$dim1$

; x-axis data
  y2 = dim_avg_n_Wrap(cndata,0)
  y2data = (/y2/)

;printVarSummary(varvar_)
;print(x2data)
;print(y2data)
;print(x1data)
;print(y1data)

;-------------------
; set work station
;-------------------
  wks = gsn_open_wks(wks_type,wks_name)
;
;  gsn_define_colormap(wks,"gui_default")
;  gsn_draw_colormap(wks)
;  gsn_merge_colormaps(wks,"amwg","gui_default")

;drawNDCGrid(wks)


;***********
; share res
;***********
  res = True
  res@gsnFrame = False
;  res@trYReverse = True

;**************
; contour res
;**************

  cnres = res
  cnres@vpHeightF = 0.5
  cnres@vpWidthF = 0.6
  cnres@vpXF = 0.25
  cnres@vpYF = 0.8

  ;cnres@trXReverse = True
  ;cnres@trYReverse = True
  ;cnres@gsnYAxisIrregular2Linear = True

;  cnres@lbOrientation = "Vertical"
  cnres@tiYAxisString = ""
  cnres@tmYLLabelsOn = False
  cnres@tmYRLabelsOn = True 
  cnres@tmXTOn = False
  ;cnres@tmXBLabelsOn = False
  cnres@cnLinesOn = False

  cnres@tmXBMode = "Explicit"
  cnres@tmXBValues = cndata&time
  cnres@trXMinF = 0.5
  cnres@trXMaxF = 12.5
  cnres@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)

;*** adjust contour levels ***
;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnMinLevelValF = 
;  res@cnMaxLevelValF = 
;  res@cnLevelSpacingF = 
;
;  res@cnLevelSelectionMode = "ExplicitLevels"
;  res@cnLevels = (//)

;*** adjust contour fill ***
  cnres@cnFillOn = True
  cnres@cnFillMode = "RasterFill"
;  res@cnMissingValFillColor = 28  ; light gray when color map is "default"
;  res@cnFillColors = 
;  res@gsnSpreadColors = True
;  res@gsnSpreadColorStart = 1
;  res@gsnSpreadColorEnd = -1
;
  cnres@cnFillPalette = "rainbow"
;  res@cnSpanFillPalette = False

;*** set title string ***
  cnres@tiMainString = "Percent of 350-630 cm~S~-1~N~ band over the broadband OLR" 
;  res@tiXAxisString =
;  res@tiYAxisString =
;  res@gsnLeftString =
;  res@gsnLeftStringFontHeightF =
;  res@gsnCenterString =
;  res@gsnRightString =


;**********
; xy1 res
;**********
  xy1res = res
  xy1res@vpWidthF = 0.1
  xy1res@vpXF = 0.1
  xy1res@tmXTOn = False
  xy1res@tmYROn = False
  xy1res@tmLabelAutoStride = True

  xy1res@xyLineThicknessF = 3.

;*** adjust xy dash pattern ***
;  xy1res@xyMonoDashPattern = False
;  xy1res@xyDashPattern = 1
;  xy1res@xyDashPatterns = (/1,2,3/)
;
;*** adjust xy line color ***
;  xy1res@xyMonoLineColor = False
;  xy1res@xyLineColors = (/1,3,3,14,14/)
;  xy1res@xyLineColor = 10

;*** adjust string ***
;  xy1res@tiMainString =
;  xy1res@tiXAxisString = 
;  xy1res@tiYAxisString = 
  xy1res@tiYAxisString = "Latitude"
;  xy1res@gsnLeftString =
;  xy1res@gsnLeftStringFontHeightF =
;  xy1res@gsnCenterString =
;  xy1res@gsnRightString =

;**********
; xy1 res
;**********
  xy2res = res

  xy2res@vpHeightF = 0.15
  xy2res@vpYF = 0.28

  xy2res@tmXTOn = False
  xy2res@tmYLOn = False
  xy2res@tmYRLabelsOn = True

  xy2res@trXMinF = min(x2data)
  xy2res@trXMaxF = max(x2data)

;*** adjust xy dash pattern ***
;  xy1res@xyMonoDashPattern = False
;  xy1res@xyDashPattern = 1
;  xy1res@xyDashPatterns = (/1,2,3/)
;
;*** adjust xy line color ***
;  xy1res@xyMonoLineColor = False
;  xy1res@xyLineColors = (/1,3,3,14,14/)
;  xy1res@xyLineColor = 10

;*** adjust string ***
;  xy1res@tiMainString =
;  xy1res@tiXAxisString = 
;  xy1res@tiYAxisString = 
  xy2res@tiXAxisString = dim1
;  xy1res@gsnLeftString =
;  xy1res@gsnLeftStringFontHeightF =
;  xy1res@gsnCenterString =
;  xy1res@gsnRightString =


;------
; plot
;------

; plot contour 
  plot = gsn_csm_contour(wks,cndata,cnres) 

; get height,width NDC coordinate
  getvalues plot
    "vpHeightF" : cn_height
    "vpWidthF"  : cn_width
    "vpXF"      : cn_xf
    "vpYF"      : cn_yf
  end getvalues

; plot xy  
  xy1res@vpHeightF = cn_height
  xy1res@vpYF = cn_yf
  plot = gsn_csm_xy(wks,x1data,y1data,xy1res)

;  xy2res@vpWidthF = cn_width
;  xy2res@vpXF = cn_xf
;  plot = gsn_csm_xy(wks,x2data,y2data,xy2res)

; plot text
;  tires = True
;  tires@txFontHeightF = 0.015

;  tires@txFontColor = 15
;  gsn_text_ndc(wks,"x1 ",0.15,0.77,tires)

;  tires = True
;  tires@txFontColor = 1
;  gsn_text_ndc(wks,"x2",0.15,0.74,tires)

  frame(wks)

  opt1 = True
  trim_fig(wks_type, wks_name, opt1)

end
