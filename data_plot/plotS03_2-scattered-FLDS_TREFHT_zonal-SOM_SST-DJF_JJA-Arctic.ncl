load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "/glade/u/home/yihsuan//script/ncl//ignorantNCL.ncl"
load "./ignorantNCL.ncl"
;load "/glade/u/home/yihsuan//script/ncl//ignorantNCL-colormaps.ncl"
;load "zz-read_parameters.ncl"

begin

;---------
; setting
;---------

  wks_type = "eps"
  ;wks_type = "ps"
  ;wks_type@wkOrientation = "landscape"

  plotname_region = "Arctic"

  wks_name = "figS03-scattered-FLDS_TREFHT_TMQ_FLDS_zonal-SOM_SST-DJF_JJA-"+plotname_region
  wks = gsn_open_wks(wks_type, wks_name)

  ;outfile_scattered = "data_plotS03-scattered-FLDS_TREFHT_zonal-SOM_SST-DJF_JJA-"+plotname_region+".nc"
  outfile_scattered = "data_plotS03-scattered-FLDS_TREFHT_TMQ_FLDS_zonal-SOM_SST-DJF_JJA-"+plotname_region+".nc"

  VAR1x_name = "FLDS"
  VAR1y_name = "TREFHT"

  VAR2x_name = "TMQ"
  VAR2y_name = "FLDS"

  VAR1_xstring = "~F33~D~F21~FLDS (W m~S~-2~N~)"
  VAR1_ystring = "~F33~D~F21~SAT (K)"

  VAR2_xstring = "~F33~D~F21~TCWV (mm)"
  VAR2_ystring = VAR1_xstring

  suffix_DJF = "DJF"
  suffix_JJA = "JJA"

  suffix_SOM = "SOM"
  suffix_SST = "SST"

;------------------
; open files and read variables
;------------------

;**********************
  out_scattered = addfile(outfile_scattered,"r")

  ;*** VAR1, DJF ***
  VAR1x_ensmean_SOM_ScatMnoScat_DJF = out_scattered->VAR1x_ensmean_SOM_ScatMnoScat_DJF
  VAR1y_ensmean_SOM_ScatMnoScat_DJF = out_scattered->VAR1y_ensmean_SOM_ScatMnoScat_DJF
  VAR1x_ensmean_SST_ScatMnoScat_DJF = out_scattered->VAR1x_ensmean_SST_ScatMnoScat_DJF
  VAR1y_ensmean_SST_ScatMnoScat_DJF = out_scattered->VAR1y_ensmean_SST_ScatMnoScat_DJF

  ;*** VAR1, JJA ***
  VAR1x_ensmean_SOM_ScatMnoScat_JJA = out_scattered->VAR1x_ensmean_SOM_ScatMnoScat_JJA
  VAR1y_ensmean_SOM_ScatMnoScat_JJA = out_scattered->VAR1y_ensmean_SOM_ScatMnoScat_JJA
  VAR1x_ensmean_SST_ScatMnoScat_JJA = out_scattered->VAR1x_ensmean_SST_ScatMnoScat_JJA
  VAR1y_ensmean_SST_ScatMnoScat_JJA = out_scattered->VAR1y_ensmean_SST_ScatMnoScat_JJA

  ;*** VAR2, DJF ***
  VAR2x_ensmean_SOM_ScatMnoScat_DJF = out_scattered->VAR2x_ensmean_SOM_ScatMnoScat_DJF
  VAR2y_ensmean_SOM_ScatMnoScat_DJF = out_scattered->VAR2y_ensmean_SOM_ScatMnoScat_DJF
  VAR2x_ensmean_SST_ScatMnoScat_DJF = out_scattered->VAR2x_ensmean_SST_ScatMnoScat_DJF
  VAR2y_ensmean_SST_ScatMnoScat_DJF = out_scattered->VAR2y_ensmean_SST_ScatMnoScat_DJF

  ;*** VAR2, JJA ***
  VAR2x_ensmean_SOM_ScatMnoScat_JJA = out_scattered->VAR2x_ensmean_SOM_ScatMnoScat_JJA
  VAR2y_ensmean_SOM_ScatMnoScat_JJA = out_scattered->VAR2y_ensmean_SOM_ScatMnoScat_JJA
  VAR2x_ensmean_SST_ScatMnoScat_JJA = out_scattered->VAR2x_ensmean_SST_ScatMnoScat_JJA
  VAR2y_ensmean_SST_ScatMnoScat_JJA = out_scattered->VAR2y_ensmean_SST_ScatMnoScat_JJA

;----------
; regline
;----------

  opt_regline = True
    string_input = ""
    opt_regline@printout = False

  string_input = ""

;****************
; VAR1, ensmean
;****************

  ;*** VAR1, SOM, DJF ***
  regline_ensmean_SOM_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_ensmean_SOM_ScatMnoScat_DJF,VAR1y_ensmean_SOM_ScatMnoScat_DJF,string_input,opt_regline)
  string_ensmean_SOM_VAR1_ScatMnoScat_DJF = regline_ensmean_SOM_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_ensmean_SOM_VAR1_ScatMnoScat_DJF)
  ;print(regline_ensmean_SOM_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_ensmean_SOM_VAR1_DJF = regline_ensmean_SOM_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_ensmean_SOM_VAR1_DJF = regline_ensmean_SOM_VAR1_ScatMnoScat_DJF@line_yy

  ;*** VAR1, SST, DJF ***
  regline_ensmean_SST_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_ensmean_SST_ScatMnoScat_DJF,VAR1y_ensmean_SST_ScatMnoScat_DJF,string_input,opt_regline)
  string_ensmean_SST_VAR1_ScatMnoScat_DJF = regline_ensmean_SST_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_ensmean_SST_VAR1_ScatMnoScat_DJF)
  ;print(regline_ensmean_SST_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_ensmean_SST_VAR1_DJF = regline_ensmean_SST_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_ensmean_SST_VAR1_DJF = regline_ensmean_SST_VAR1_ScatMnoScat_DJF@line_yy

  VAR1_plot_xx_DJF = (/VAR1x_ensmean_SOM_ScatMnoScat_DJF, VAR1x_ensmean_SST_ScatMnoScat_DJF/) 
  VAR1_plot_yy_DJF = (/VAR1y_ensmean_SOM_ScatMnoScat_DJF, VAR1y_ensmean_SST_ScatMnoScat_DJF/)

  labels_VAR1_DJF = (/suffix_SOM+": "+regline_ensmean_SOM_VAR1_ScatMnoScat_DJF@string_fit, \
                      suffix_SST+": "+regline_ensmean_SST_VAR1_ScatMnoScat_DJF@string_fit/)

  ;*** VAR1, SOM, JJA ***
  regline_ensmean_SOM_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_ensmean_SOM_ScatMnoScat_JJA,VAR1y_ensmean_SOM_ScatMnoScat_JJA,string_input,opt_regline)
  string_ensmean_SOM_VAR1_ScatMnoScat_JJA = regline_ensmean_SOM_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_ensmean_SOM_VAR1_ScatMnoScat_JJA)
  ;print(regline_ensmean_SOM_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_ensmean_SOM_VAR1_JJA = regline_ensmean_SOM_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_ensmean_SOM_VAR1_JJA = regline_ensmean_SOM_VAR1_ScatMnoScat_JJA@line_yy

  ;*** VAR1, SST, JJA ***
  regline_ensmean_SST_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_ensmean_SST_ScatMnoScat_JJA,VAR1y_ensmean_SST_ScatMnoScat_JJA,string_input,opt_regline)
  string_ensmean_SST_VAR1_ScatMnoScat_JJA = regline_ensmean_SST_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_ensmean_SST_VAR1_ScatMnoScat_JJA)
  ;print(regline_ensmean_SST_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_ensmean_SST_VAR1_JJA = regline_ensmean_SST_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_ensmean_SST_VAR1_JJA = regline_ensmean_SST_VAR1_ScatMnoScat_JJA@line_yy

  VAR1_plot_xx_JJA = (/VAR1x_ensmean_SOM_ScatMnoScat_JJA, VAR1x_ensmean_SST_ScatMnoScat_JJA/) 
  VAR1_plot_yy_JJA = (/VAR1y_ensmean_SOM_ScatMnoScat_JJA, VAR1y_ensmean_SST_ScatMnoScat_JJA/)

  VAR1_plot_xx_JJA = (/VAR1x_ensmean_SOM_ScatMnoScat_JJA, VAR1x_ensmean_SST_ScatMnoScat_JJA/) 
  VAR1_plot_yy_JJA = (/VAR1y_ensmean_SOM_ScatMnoScat_JJA, VAR1y_ensmean_SST_ScatMnoScat_JJA/)

  labels_VAR1_JJA = (/suffix_SOM+": "+regline_ensmean_SOM_VAR1_ScatMnoScat_JJA@string_fit, \
                      suffix_SST+": "+regline_ensmean_SST_VAR1_ScatMnoScat_JJA@string_fit/)

  ;*** VAR2, SOM, DJF ***
  regline_ensmean_SOM_VAR2_ScatMnoScat_DJF = regline_print(VAR2x_ensmean_SOM_ScatMnoScat_DJF,VAR2y_ensmean_SOM_ScatMnoScat_DJF,string_input,opt_regline)
  string_ensmean_SOM_VAR2_ScatMnoScat_DJF = regline_ensmean_SOM_VAR2_ScatMnoScat_DJF@string_fit
  ;print(regline_ensmean_SOM_VAR2_ScatMnoScat_DJF)
  ;print(regline_ensmean_SOM_VAR2_ScatMnoScat_DJF@string_fit)
  line_xx_ensmean_SOM_VAR2_DJF = regline_ensmean_SOM_VAR2_ScatMnoScat_DJF@line_xx
  line_yy_ensmean_SOM_VAR2_DJF = regline_ensmean_SOM_VAR2_ScatMnoScat_DJF@line_yy

  ;*** VAR2, SST, DJF ***
  regline_ensmean_SST_VAR2_ScatMnoScat_DJF = regline_print(VAR2x_ensmean_SST_ScatMnoScat_DJF,VAR2y_ensmean_SST_ScatMnoScat_DJF,string_input,opt_regline)
  string_ensmean_SST_VAR2_ScatMnoScat_DJF = regline_ensmean_SST_VAR2_ScatMnoScat_DJF@string_fit
  ;print(regline_ensmean_SST_VAR2_ScatMnoScat_DJF)
  ;print(regline_ensmean_SST_VAR2_ScatMnoScat_DJF@string_fit)
  line_xx_ensmean_SST_VAR2_DJF = regline_ensmean_SST_VAR2_ScatMnoScat_DJF@line_xx
  line_yy_ensmean_SST_VAR2_DJF = regline_ensmean_SST_VAR2_ScatMnoScat_DJF@line_yy

  VAR2_plot_xx_DJF = (/VAR2x_ensmean_SOM_ScatMnoScat_DJF, VAR2x_ensmean_SST_ScatMnoScat_DJF/) 
  VAR2_plot_yy_DJF = (/VAR2y_ensmean_SOM_ScatMnoScat_DJF, VAR2y_ensmean_SST_ScatMnoScat_DJF/)

  labels_VAR2_DJF = (/suffix_SOM+": "+regline_ensmean_SOM_VAR2_ScatMnoScat_DJF@string_fit, \
                      suffix_SST+": "+regline_ensmean_SST_VAR2_ScatMnoScat_DJF@string_fit/)

  ;*** VAR2, SOM, JJA ***
  regline_ensmean_SOM_VAR2_ScatMnoScat_JJA = regline_print(VAR2x_ensmean_SOM_ScatMnoScat_JJA,VAR2y_ensmean_SOM_ScatMnoScat_JJA,string_input,opt_regline)
  string_ensmean_SOM_VAR2_ScatMnoScat_JJA = regline_ensmean_SOM_VAR2_ScatMnoScat_JJA@string_fit
  ;print(regline_ensmean_SOM_VAR2_ScatMnoScat_JJA)
  ;print(regline_ensmean_SOM_VAR2_ScatMnoScat_JJA@string_fit)
  line_xx_ensmean_SOM_VAR2_JJA = regline_ensmean_SOM_VAR2_ScatMnoScat_JJA@line_xx
  line_yy_ensmean_SOM_VAR2_JJA = regline_ensmean_SOM_VAR2_ScatMnoScat_JJA@line_yy

  ;*** VAR2, SST, JJA ***
  regline_ensmean_SST_VAR2_ScatMnoScat_JJA = regline_print(VAR2x_ensmean_SST_ScatMnoScat_JJA,VAR2y_ensmean_SST_ScatMnoScat_JJA,string_input,opt_regline)
  string_ensmean_SST_VAR2_ScatMnoScat_JJA = regline_ensmean_SST_VAR2_ScatMnoScat_JJA@string_fit
  ;print(regline_ensmean_SST_VAR2_ScatMnoScat_JJA)
  ;print(regline_ensmean_SST_VAR2_ScatMnoScat_JJA@string_fit)
  line_xx_ensmean_SST_VAR2_JJA = regline_ensmean_SST_VAR2_ScatMnoScat_JJA@line_xx
  line_yy_ensmean_SST_VAR2_JJA = regline_ensmean_SST_VAR2_ScatMnoScat_JJA@line_yy

  VAR2_plot_xx_JJA = (/VAR2x_ensmean_SOM_ScatMnoScat_JJA, VAR2x_ensmean_SST_ScatMnoScat_JJA/) 
  VAR2_plot_yy_JJA = (/VAR2y_ensmean_SOM_ScatMnoScat_JJA, VAR2y_ensmean_SST_ScatMnoScat_JJA/)

  VAR2_plot_xx_JJA = (/VAR2x_ensmean_SOM_ScatMnoScat_JJA, VAR2x_ensmean_SST_ScatMnoScat_JJA/) 
  VAR2_plot_yy_JJA = (/VAR2y_ensmean_SOM_ScatMnoScat_JJA, VAR2y_ensmean_SST_ScatMnoScat_JJA/)

  labels_VAR2_JJA = (/suffix_SOM+": "+regline_ensmean_SOM_VAR2_ScatMnoScat_JJA@string_fit, \
                      suffix_SST+": "+regline_ensmean_SST_VAR2_ScatMnoScat_JJA@string_fit/)

;-----------------
; plot scattered 
;-----------------
  opt_xy_scat = True
    opt_xy_scat@legend_mean = False
    labels = (/"SOM","SST"/)

  lnres = True
  lnres@gsLineThicknessF = 3.

  lnres_SOM = lnres
  lnres_SOM@gsLineColor = "red"
  lnres_SST = lnres
  lnres_SST@gsLineColor = "blue"

  res_scat  = True
  res_scat@gsnFrame = False
  res_scat@gsnDraw = False

  ;res_scat@lgLabelFontHeightF = 0.015

  ;res_scat@pmLegendDisplayMode = "nocreate"
  res_scat@pmLegendWidthF = 0.11
  res_scat@lgLabelFontHeightF = 0.015

  res_scat@xyMarker = 4                ; 1: dot, 4: empty circle, 16: full circle
  ;res_scat@xyMarker = 16                ; 1: dot, 4: empty circle, 16: full circle
  res_scat@gsnXRefLine = 0.
  res_scat@gsnYRefLine = 0.
  res_scat@gsnYRefLineDashPattern = 2
  res_scat@gsnXRefLineDashPattern = res_scat@gsnYRefLineDashPattern

  res_scat@xyMonoMarkLineMode = True
  res_scat@xyMarkLineMode    = "Markers"
  
  res_scat@xyMarkerColors = (/"red","blue"/)
  res_scat@xyLineColors = res_scat@xyMarkerColors

  ;*** VAR1 ***
  res_scat_VAR1 = res_scat
  res_scat_VAR1@tiXAxisString = VAR1_xstring
  res_scat_VAR1@tiYAxisString = VAR1_ystring
  res_scat_VAR1@trXMinF = -12.
  res_scat_VAR1@trXMaxF = 16.
  res_scat_VAR1@trYMinF = -4.
  res_scat_VAR1@trYMaxF = 5.


  ;*** VAR1, DJF ***
  res_scat_VAR1_DJF = res_scat_VAR1
  res_scat_VAR1_DJF@tiMainString = suffix_DJF 
  plot_scat_ensmean_VAR1_DJF = ignorant_xy(wks, VAR1_plot_xx_DJF, VAR1_plot_yy_DJF, labels_VAR1_DJF, res_scat_VAR1_DJF, opt_xy_scat)

  line_SOM_VAR1_DJF = gsn_add_polyline(wks, plot_scat_ensmean_VAR1_DJF, line_xx_ensmean_SOM_VAR1_DJF, line_yy_ensmean_SOM_VAR1_DJF, lnres_SOM)
  line_SST_VAR1_DJF = gsn_add_polyline(wks, plot_scat_ensmean_VAR1_DJF, line_xx_ensmean_SST_VAR1_DJF, line_yy_ensmean_SST_VAR1_DJF, lnres_SST)

  ;*** VAR1, JJA ***
  res_scat_VAR1_JJA = res_scat_VAR1
  res_scat_VAR1_JJA@tiMainString = suffix_JJA 
  plot_scat_ensmean_VAR1_JJA = ignorant_xy(wks, VAR1_plot_xx_JJA, VAR1_plot_yy_JJA, labels_VAR1_JJA, res_scat_VAR1_JJA, opt_xy_scat)

  line_SOM_VAR1_JJA = gsn_add_polyline(wks, plot_scat_ensmean_VAR1_JJA, line_xx_ensmean_SOM_VAR1_JJA, line_yy_ensmean_SOM_VAR1_JJA, lnres_SOM)
  line_SST_VAR1_JJA = gsn_add_polyline(wks, plot_scat_ensmean_VAR1_JJA, line_xx_ensmean_SST_VAR1_JJA, line_yy_ensmean_SST_VAR1_JJA, lnres_SST)

  ;*** VAR2 ***
  ;*** VAR2 ***
  res_scat_VAR2 = res_scat
  res_scat_VAR2@tiXAxisString = VAR2_xstring
  res_scat_VAR2@tiYAxisString = VAR2_ystring
  res_scat_VAR2@trXMinF = -1.5 
  res_scat_VAR2@trXMaxF = 1.5
  res_scat_VAR2@trYMinF = res_scat_VAR1@trXMinF
  res_scat_VAR2@trYMaxF = res_scat_VAR1@trXMaxF

  ;*** VAR2, DJF ***
  res_scat_VAR2_DJF = res_scat_VAR2
  res_scat_VAR2_DJF@tiMainString = suffix_DJF 
  plot_scat_ensmean_VAR2_DJF = ignorant_xy(wks, VAR2_plot_xx_DJF, VAR2_plot_yy_DJF, labels_VAR2_DJF, res_scat_VAR2_DJF, opt_xy_scat)

  line_SOM_VAR2_DJF = gsn_add_polyline(wks, plot_scat_ensmean_VAR2_DJF, line_xx_ensmean_SOM_VAR2_DJF, line_yy_ensmean_SOM_VAR2_DJF, lnres_SOM)
  line_SST_VAR2_DJF = gsn_add_polyline(wks, plot_scat_ensmean_VAR2_DJF, line_xx_ensmean_SST_VAR2_DJF, line_yy_ensmean_SST_VAR2_DJF, lnres_SST)

  ;*** VAR2, JJA ***
  res_scat_VAR2_JJA = res_scat_VAR2
  res_scat_VAR2_JJA@tiMainString = suffix_JJA 
  plot_scat_ensmean_VAR2_JJA = ignorant_xy(wks, VAR2_plot_xx_JJA, VAR2_plot_yy_JJA, labels_VAR2_JJA, res_scat_VAR2_JJA, opt_xy_scat)

  line_SOM_VAR2_JJA = gsn_add_polyline(wks, plot_scat_ensmean_VAR2_JJA, line_xx_ensmean_SOM_VAR2_JJA, line_yy_ensmean_SOM_VAR2_JJA, lnres_SOM)
  line_SST_VAR2_JJA = gsn_add_polyline(wks, plot_scat_ensmean_VAR2_JJA, line_xx_ensmean_SST_VAR2_JJA, line_yy_ensmean_SST_VAR2_JJA, lnres_SST)

;-------------
; panel plot
;-------------
  nplot = 4

  plot_all = new(nplot, "graphic")

  plot_all(0) = plot_scat_ensmean_VAR1_DJF
  plot_all(1) = plot_scat_ensmean_VAR1_JJA
  plot_all(2) = plot_scat_ensmean_VAR2_DJF
  plot_all(3) = plot_scat_ensmean_VAR2_JJA

  resp = True
  resp@gsnMaximize = True
  column1 = 2
  row1 = nplot / column1
  gsn_panel(wks, plot_all, (/row1,column1/), resp)

  ;opt1 = True
  ;trim_fig(wks_type, wks_name, opt1)

end
