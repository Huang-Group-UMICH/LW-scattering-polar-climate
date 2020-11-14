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

  outfile_scattered = "data_plot03-scattered-FLDS_TREFHT_TMQ_FLDS_areamean-SOM_SST-DJF_JJA-Arctic.nc"
  outfile_profile = "data_plot03-profile-T_Q_CLOUD_areamean-SOM_SST-DJF_JJA-Arctic.nc"

  wks_type = "eps"

  plotname_region = "Arctic"

  ;wks_name = "test111_"+plotname_region
  wks_name = "fig03-xy_scatter_profile-VARs-Arctic-DJF_JJA"
  wks = gsn_open_wks(wks_type, wks_name)

  VAR1x_name = "FLDS"
  VAR1y_name = "TREFHT"

  VAR2x_name = "TMQ"
  VAR2y_name = "FLDS"

  VAR1_profile = "T"
  VAR2_profile = "Q"
  VAR3_profile = "CLOUD"

  opt00 = True
    opt00@string_change = "delta"

  xstring_VAR1x = string_change("FLDS (W m~S~-2~N~)",opt00)
  ystring_VAR1y = string_change("SAT (K)",opt00)
  xstring_VAR2x = string_change("TCWV (mm)",opt00)
  ystring_VAR2y = xstring_VAR1x

  xstring_VAR1 = string_change("T (K)",opt00)
  xstring_VAR2 = string_change("Q (g/kg)",opt00)
  xstring_VAR3 = string_change("CLOUD (%)",opt00)

  suffix_DJF = "DJF"
  suffix_JJA = "JJA"

  suffix_SOM = "SOM"
  suffix_SST = "SST"
  suffix_SOMSST = (/suffix_SOM,suffix_SST/)
  suffix_SSTSOM = (/suffix_SST,suffix_SOM/)

;------------------
; open files and read variables
;------------------

;**********************
  out_scattered = addfile(outfile_scattered,"r")

  ;*** VAR1, DJF ***
  VAR1x_SOM_ScatMnoScat_DJF = out_scattered->VAR1x_SOM_ScatMnoScat_DJF
  VAR1y_SOM_ScatMnoScat_DJF = out_scattered->VAR1y_SOM_ScatMnoScat_DJF
  VAR1x_SST_ScatMnoScat_DJF = out_scattered->VAR1x_SST_ScatMnoScat_DJF
  VAR1y_SST_ScatMnoScat_DJF = out_scattered->VAR1y_SST_ScatMnoScat_DJF

  ;*** VAR2, DJF ***
  VAR2x_SOM_ScatMnoScat_DJF = out_scattered->VAR2x_SOM_ScatMnoScat_DJF
  VAR2y_SOM_ScatMnoScat_DJF = out_scattered->VAR2y_SOM_ScatMnoScat_DJF
  VAR2x_SST_ScatMnoScat_DJF = out_scattered->VAR2x_SST_ScatMnoScat_DJF
  VAR2y_SST_ScatMnoScat_DJF = out_scattered->VAR2y_SST_ScatMnoScat_DJF

  ;*** VAR1, JJA ***
  VAR1x_SOM_ScatMnoScat_JJA = out_scattered->VAR1x_SOM_ScatMnoScat_JJA
  VAR1y_SOM_ScatMnoScat_JJA = out_scattered->VAR1y_SOM_ScatMnoScat_JJA
  VAR1x_SST_ScatMnoScat_JJA = out_scattered->VAR1x_SST_ScatMnoScat_JJA
  VAR1y_SST_ScatMnoScat_JJA = out_scattered->VAR1y_SST_ScatMnoScat_JJA

  ;*** VAR2, JJA ***
  VAR2x_SOM_ScatMnoScat_JJA = out_scattered->VAR2x_SOM_ScatMnoScat_JJA
  VAR2y_SOM_ScatMnoScat_JJA = out_scattered->VAR2y_SOM_ScatMnoScat_JJA
  VAR2x_SST_ScatMnoScat_JJA = out_scattered->VAR2x_SST_ScatMnoScat_JJA
  VAR2y_SST_ScatMnoScat_JJA = out_scattered->VAR2y_SST_ScatMnoScat_JJA

;**********************
  out_profile = addfile(outfile_profile,"r")

  VAR1p_ScatMnoScat_SOM_DJF = out_profile->VAR1p_ScatMnoScat_SOM_DJF
  VAR1p_ScatMnoScat_SST_DJF = out_profile->VAR1p_ScatMnoScat_SST_DJF
  VAR2p_ScatMnoScat_SOM_DJF = out_profile->VAR2p_ScatMnoScat_SOM_DJF
  VAR2p_ScatMnoScat_SST_DJF = out_profile->VAR2p_ScatMnoScat_SST_DJF
  VAR3p_ScatMnoScat_SOM_DJF = out_profile->VAR3p_ScatMnoScat_SOM_DJF
  VAR3p_ScatMnoScat_SST_DJF = out_profile->VAR3p_ScatMnoScat_SST_DJF

  VAR1p_ScatMnoScat_SOM_JJA = out_profile->VAR1p_ScatMnoScat_SOM_JJA
  VAR1p_ScatMnoScat_SST_JJA = out_profile->VAR1p_ScatMnoScat_SST_JJA
  VAR2p_ScatMnoScat_SOM_JJA = out_profile->VAR2p_ScatMnoScat_SOM_JJA
  VAR2p_ScatMnoScat_SST_JJA = out_profile->VAR2p_ScatMnoScat_SST_JJA
  VAR3p_ScatMnoScat_SOM_JJA = out_profile->VAR3p_ScatMnoScat_SOM_JJA
  VAR3p_ScatMnoScat_SST_JJA = out_profile->VAR3p_ScatMnoScat_SST_JJA

  pres = VAR3p_ScatMnoScat_SST_JJA&coord
  ;print(pres)

  nn = dimsizes(VAR1p_ScatMnoScat_SST_JJA)
  n1 = nn(0)
  n2 = nn(1)
  nfile = 2
  nmember = n1

  VAR1p_ScatMnoScat_SOMSST_DJF = reshape( (/VAR1p_ScatMnoScat_SOM_DJF, VAR1p_ScatMnoScat_SST_DJF/), (/nfile*n1,n2/) )
  VAR2p_ScatMnoScat_SOMSST_DJF = reshape( (/VAR2p_ScatMnoScat_SOM_DJF, VAR2p_ScatMnoScat_SST_DJF/), (/nfile*n1,n2/) )
  VAR3p_ScatMnoScat_SOMSST_DJF = reshape( (/VAR3p_ScatMnoScat_SOM_DJF, VAR3p_ScatMnoScat_SST_DJF/), (/nfile*n1,n2/) )

  VAR1p_ScatMnoScat_SOMSST_JJA = reshape( (/VAR1p_ScatMnoScat_SOM_JJA, VAR1p_ScatMnoScat_SST_JJA/), (/nfile*n1,n2/) )
  VAR2p_ScatMnoScat_SOMSST_JJA = reshape( (/VAR2p_ScatMnoScat_SOM_JJA, VAR2p_ScatMnoScat_SST_JJA/), (/nfile*n1,n2/) )
  VAR3p_ScatMnoScat_SOMSST_JJA = reshape( (/VAR3p_ScatMnoScat_SOM_JJA, VAR3p_ScatMnoScat_SST_JJA/), (/nfile*n1,n2/) )

;----------
; regline
;----------

  opt_regline = True
    string_input = ""
    opt_regline@printout = False

  string_input = ""

;*******
; VAR1
;*******

  ;*** VAR1, SOM, DJF ***
  regline_SOM_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_SOM_ScatMnoScat_DJF,VAR1y_SOM_ScatMnoScat_DJF,string_input,opt_regline)
  string_SOM_VAR1_ScatMnoScat_DJF = regline_SOM_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_SOM_VAR1_ScatMnoScat_DJF)
  ;print(regline_SOM_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_SOM_VAR1_DJF = regline_SOM_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_SOM_VAR1_DJF = regline_SOM_VAR1_ScatMnoScat_DJF@line_yy

  ;*** VAR1, SST, DJF ***
  regline_SST_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_SST_ScatMnoScat_DJF,VAR1y_SST_ScatMnoScat_DJF,string_input,opt_regline)
  string_SST_VAR1_ScatMnoScat_DJF = regline_SST_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_SST_VAR1_ScatMnoScat_DJF)
  ;print(regline_SST_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_SST_VAR1_DJF = regline_SST_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_SST_VAR1_DJF = regline_SST_VAR1_ScatMnoScat_DJF@line_yy

  VAR1_plot_xx_DJF = (/VAR1x_SOM_ScatMnoScat_DJF, VAR1x_SST_ScatMnoScat_DJF/) 
  VAR1_plot_yy_DJF = (/VAR1y_SOM_ScatMnoScat_DJF, VAR1y_SST_ScatMnoScat_DJF/)

  labels_VAR1_DJF = (/suffix_SOM+": "+regline_SOM_VAR1_ScatMnoScat_DJF@string_fit, \
                      suffix_SST+": "+regline_SST_VAR1_ScatMnoScat_DJF@string_fit/)


  ;*** VAR1, SOM, JJA ***
  regline_SOM_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_SOM_ScatMnoScat_JJA,VAR1y_SOM_ScatMnoScat_JJA,string_input,opt_regline)
  string_SOM_VAR1_ScatMnoScat_JJA = regline_SOM_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_SOM_VAR1_ScatMnoScat_JJA)
  ;print(regline_SOM_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_SOM_VAR1_JJA = regline_SOM_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_SOM_VAR1_JJA = regline_SOM_VAR1_ScatMnoScat_JJA@line_yy

  ;*** VAR1, SST, JJA ***
  regline_SST_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_SST_ScatMnoScat_JJA,VAR1y_SST_ScatMnoScat_JJA,string_input,opt_regline)
  string_SST_VAR1_ScatMnoScat_JJA = regline_SST_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_SST_VAR1_ScatMnoScat_JJA)
  ;print(regline_SST_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_SST_VAR1_JJA = regline_SST_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_SST_VAR1_JJA = regline_SST_VAR1_ScatMnoScat_JJA@line_yy

  VAR1_plot_xx_JJA = (/VAR1x_SOM_ScatMnoScat_JJA, VAR1x_SST_ScatMnoScat_JJA/) 
  VAR1_plot_yy_JJA = (/VAR1y_SOM_ScatMnoScat_JJA, VAR1y_SST_ScatMnoScat_JJA/)

  VAR1_plot_xx_JJA = (/VAR1x_SOM_ScatMnoScat_JJA, VAR1x_SST_ScatMnoScat_JJA/) 
  VAR1_plot_yy_JJA = (/VAR1y_SOM_ScatMnoScat_JJA, VAR1y_SST_ScatMnoScat_JJA/)

  labels_VAR1_JJA = (/suffix_SOM+": "+regline_SOM_VAR1_ScatMnoScat_JJA@string_fit, \
                      suffix_SST+": "+regline_SST_VAR1_ScatMnoScat_JJA@string_fit/)

;*******
; VAR2
;*******

  ;*** VAR2, SOM, DJF ***
  regline_SOM_VAR2_ScatMnoScat_DJF = regline_print(VAR2x_SOM_ScatMnoScat_DJF,VAR2y_SOM_ScatMnoScat_DJF,string_input,opt_regline)
  string_SOM_VAR2_ScatMnoScat_DJF = regline_SOM_VAR2_ScatMnoScat_DJF@string_fit
  ;print(regline_SOM_VAR2_ScatMnoScat_DJF)
  ;print(regline_SOM_VAR2_ScatMnoScat_DJF@string_fit)
  line_xx_SOM_VAR2_DJF = regline_SOM_VAR2_ScatMnoScat_DJF@line_xx
  line_yy_SOM_VAR2_DJF = regline_SOM_VAR2_ScatMnoScat_DJF@line_yy

  ;*** VAR2, SST, DJF ***
  regline_SST_VAR2_ScatMnoScat_DJF = regline_print(VAR2x_SST_ScatMnoScat_DJF,VAR2y_SST_ScatMnoScat_DJF,string_input,opt_regline)
  string_SST_VAR2_ScatMnoScat_DJF = regline_SST_VAR2_ScatMnoScat_DJF@string_fit
  ;print(regline_SST_VAR2_ScatMnoScat_DJF)
  ;print(regline_SST_VAR2_ScatMnoScat_DJF@string_fit)
  line_xx_SST_VAR2_DJF = regline_SST_VAR2_ScatMnoScat_DJF@line_xx
  line_yy_SST_VAR2_DJF = regline_SST_VAR2_ScatMnoScat_DJF@line_yy

  VAR2_plot_xx_DJF = (/VAR2x_SOM_ScatMnoScat_DJF, VAR2x_SST_ScatMnoScat_DJF/) 
  VAR2_plot_yy_DJF = (/VAR2y_SOM_ScatMnoScat_DJF, VAR2y_SST_ScatMnoScat_DJF/)

  labels_VAR2_DJF = (/suffix_SOM+": "+regline_SOM_VAR2_ScatMnoScat_DJF@string_fit, \
                      suffix_SST+": "+regline_SST_VAR2_ScatMnoScat_DJF@string_fit/)

  ;*** VAR2, SOM, JJA ***
  regline_SOM_VAR2_ScatMnoScat_JJA = regline_print(VAR2x_SOM_ScatMnoScat_JJA,VAR2y_SOM_ScatMnoScat_JJA,string_input,opt_regline)
  string_SOM_VAR2_ScatMnoScat_JJA = regline_SOM_VAR2_ScatMnoScat_JJA@string_fit
  ;print(regline_SOM_VAR2_ScatMnoScat_JJA)
  ;print(regline_SOM_VAR2_ScatMnoScat_JJA@string_fit)
  line_xx_SOM_VAR2_JJA = regline_SOM_VAR2_ScatMnoScat_JJA@line_xx
  line_yy_SOM_VAR2_JJA = regline_SOM_VAR2_ScatMnoScat_JJA@line_yy

  ;*** VAR2, SST, JJA ***
  regline_SST_VAR2_ScatMnoScat_JJA = regline_print(VAR2x_SST_ScatMnoScat_JJA,VAR2y_SST_ScatMnoScat_JJA,string_input,opt_regline)
  string_SST_VAR2_ScatMnoScat_JJA = regline_SST_VAR2_ScatMnoScat_JJA@string_fit
  ;print(regline_SST_VAR2_ScatMnoScat_JJA)
  ;print(regline_SST_VAR2_ScatMnoScat_JJA@string_fit)
  line_xx_SST_VAR2_JJA = regline_SST_VAR2_ScatMnoScat_JJA@line_xx
  line_yy_SST_VAR2_JJA = regline_SST_VAR2_ScatMnoScat_JJA@line_yy

  VAR2_plot_xx_JJA = (/VAR2x_SOM_ScatMnoScat_JJA, VAR2x_SST_ScatMnoScat_JJA/) 
  VAR2_plot_yy_JJA = (/VAR2y_SOM_ScatMnoScat_JJA, VAR2y_SST_ScatMnoScat_JJA/)

  VAR2_plot_xx_JJA = (/VAR2x_SOM_ScatMnoScat_JJA, VAR2x_SST_ScatMnoScat_JJA/) 
  VAR2_plot_yy_JJA = (/VAR2y_SOM_ScatMnoScat_JJA, VAR2y_SST_ScatMnoScat_JJA/)

  labels_VAR2_JJA = (/suffix_SOM+": "+regline_SOM_VAR2_ScatMnoScat_JJA@string_fit, \
                      suffix_SST+": "+regline_SST_VAR2_ScatMnoScat_JJA@string_fit/)

;-----------
; share res
;-----------
  res = True
  res@gsnFrame = False
  res@gsnDraw = False
  res@pmLegendDisplayMode = "nocreate"
  
  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = res@tmXBLabelFontHeightF

  res@tiXAxisFontHeightF = res@tmXBLabelFontHeightF
  res@tiYAxisFontHeightF = res@tmXBLabelFontHeightF

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

  res_scat  = res
  ;res_scat@tiMainString = " "

  ;res_scat@pmLegendDisplayMode = "nocreate"
  ;res_scat@pmLegendWidthF = 
  ;res_scat@lgLabelFontHeightF = 0.01

  ;res_scat@xyMarker = 4                ; 1: dot, 4: empty circle, 16: full circle
  res_scat@xyMarker = 16                ; 1: dot, 4: empty circle, 16: full circle
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
  res_scat_VAR1@tiXAxisString = xstring_VAR1x
  res_scat_VAR1@tiYAxisString = ystring_VAR1y
  res_scat_VAR1@trXMinF = -9.
  res_scat_VAR1@trXMaxF = 9.
  res_scat_VAR1@trYMinF = -3.
  res_scat_VAR1@trYMaxF = 3.

  ;*** VAR1, DJF ***
  res_scat_VAR1_DJF = res_scat_VAR1
  ;res_scat_VAR1_DJF@tiMainString = suffix_DJF 
  plot_scat_VAR1_DJF = ignorant_xy(wks, VAR1_plot_xx_DJF, VAR1_plot_yy_DJF, labels_VAR1_DJF, res_scat_VAR1_DJF, opt_xy_scat)

  line_SOM_VAR1_DJF = gsn_add_polyline(wks, plot_scat_VAR1_DJF, line_xx_SOM_VAR1_DJF, line_yy_SOM_VAR1_DJF, lnres_SOM)
  line_SST_VAR1_DJF = gsn_add_polyline(wks, plot_scat_VAR1_DJF, line_xx_SST_VAR1_DJF, line_yy_SST_VAR1_DJF, lnres_SST)

  ;*** VAR1, JJA ***
  res_scat_VAR1_JJA = res_scat_VAR1
  ;res_scat_VAR1_JJA@tiMainString = suffix_JJA 
  plot_scat_VAR1_JJA = ignorant_xy(wks, VAR1_plot_xx_JJA, VAR1_plot_yy_JJA, labels_VAR1_JJA, res_scat_VAR1_JJA, opt_xy_scat)

  line_SOM_VAR1_JJA = gsn_add_polyline(wks, plot_scat_VAR1_JJA, line_xx_SOM_VAR1_JJA, line_yy_SOM_VAR1_JJA, lnres_SOM)
  line_SST_VAR1_JJA = gsn_add_polyline(wks, plot_scat_VAR1_JJA, line_xx_SST_VAR1_JJA, line_yy_SST_VAR1_JJA, lnres_SST)

  ;*** VAR2 ***
  res_scat_VAR2 = res_scat
  res_scat_VAR2@tiXAxisString = xstring_VAR2x
  res_scat_VAR2@tiYAxisString = ystring_VAR2y
  res_scat_VAR2@trXMinF = -0.8
  res_scat_VAR2@trXMaxF = 1.
  res_scat_VAR2@trYMinF = res_scat_VAR1@trXMinF
  res_scat_VAR2@trYMaxF = res_scat_VAR1@trXMaxF

  ;*** VAR2, DJF ***
  res_scat_VAR2_DJF = res_scat_VAR2
  ;res_scat_VAR2_DJF@tiMainString = suffix_DJF
  plot_scat_VAR2_DJF = ignorant_xy(wks, VAR2_plot_xx_DJF, VAR2_plot_yy_DJF, labels_VAR2_DJF, res_scat_VAR2_DJF, opt_xy_scat)

  line_SOM_VAR2_DJF = gsn_add_polyline(wks, plot_scat_VAR2_DJF, line_xx_SOM_VAR2_DJF, line_yy_SOM_VAR2_DJF, lnres_SOM)
  line_SST_VAR2_DJF = gsn_add_polyline(wks, plot_scat_VAR2_DJF, line_xx_SST_VAR2_DJF, line_yy_SST_VAR2_DJF, lnres_SST)

  ;*** VAR2, JJA ***
  res_scat_VAR2_JJA = res_scat_VAR2
  ;res_scat_VAR2_JJA@tiMainString = suffix_JJA 
  plot_scat_VAR2_JJA = ignorant_xy(wks, VAR2_plot_xx_JJA, VAR2_plot_yy_JJA, labels_VAR2_JJA, res_scat_VAR2_JJA, opt_xy_scat)

  line_SOM_VAR2_JJA = gsn_add_polyline(wks, plot_scat_VAR2_JJA, line_xx_SOM_VAR2_JJA, line_yy_SOM_VAR2_JJA, lnres_SOM)
  line_SST_VAR2_JJA = gsn_add_polyline(wks, plot_scat_VAR2_JJA, line_xx_SST_VAR2_JJA, line_yy_SST_VAR2_JJA, lnres_SST)

;===========================================================
;===========================================================
;===========================================================

;-----------------
; plot profiles
;-----------------
  opt_xy_profile = True
    opt_xy_profile@legend_mean = False
    opt_xy_profile@option_group = True     ; if grouping the first dimension
    opt_xy_profile@nmember_per_group = nmember   ; number of member per group 
    opt_xy_profile@XYStyle = "solid"
    ;opt_xy_profile@XYStyle = "fill"

    labels = (/"SOM","SST"/)

  res_profile  = res
  ;res_profile@pmLegendDisplayMode = "Always"

  res_profile@trYReverse = True
  res_profile@gsnXRefLine = 0.
  res_profile@gsnXRefLineDashPattern = res_scat@gsnYRefLineDashPattern
  res_profile@tiYAxisString = "Pressure (hPa)"
  res_profile@tiMainString = " "

  ;*** VAR 1 ***
  res_profile_VAR1 = res_profile
  res_profile_VAR1@tiXAxisString = xstring_VAR1
  res_profile_VAR1@trXMinF = -1.5
  res_profile_VAR1@trXMaxF = 2.

  ;*** VAR1, DJF ***
  res_prof_VAR1_DJF = res_profile_VAR1
  plot_prof_VAR1_DJF = ignorant_xy_ensemble(wks, VAR1p_ScatMnoScat_SOMSST_DJF, pres, labels, res_prof_VAR1_DJF, opt_xy_profile)

  ;*** VAR1, JJA ***
  res_prof_VAR1_JJA = res_profile_VAR1
  plot_prof_VAR1_JJA = ignorant_xy_ensemble(wks, VAR1p_ScatMnoScat_SOMSST_JJA, pres, labels, res_prof_VAR1_JJA, opt_xy_profile)

  ;*** VAR 2 ***
  res_profile_VAR2 = res_profile
  res_profile_VAR2@tiXAxisString = xstring_VAR2
  res_profile_VAR2@trXMinF = -0.09
  res_profile_VAR2@trXMaxF = 0.15

  ;*** VAR2, DJF ***
  res_prof_VAR2_DJF = res_profile_VAR2
  plot_prof_VAR2_DJF = ignorant_xy_ensemble(wks, VAR2p_ScatMnoScat_SOMSST_DJF, pres, labels, res_prof_VAR2_DJF, opt_xy_profile)

  ;*** VAR2, JJA ***
  res_prof_VAR2_JJA = res_profile_VAR2
  plot_prof_VAR2_JJA = ignorant_xy_ensemble(wks, VAR2p_ScatMnoScat_SOMSST_JJA, pres, labels, res_prof_VAR2_JJA, opt_xy_profile)

  ;*** VAR3 ***
  res_profile_VAR3 = res_profile
  res_profile_VAR3@tiXAxisString = xstring_VAR3
  res_profile_VAR3@trXMinF = -1.
  res_profile_VAR3@trXMaxF = 2.5

  ;res_profile_VAR3@trXMinF = -1.5
  ;res_profile_VAR3@trXMaxF = 2.5
  ;res_profile_VAR3@tmXBMode = "Manual"
  ;res_profile_VAR3@tmXBTickStartF = res_profile_VAR3@trXMinF
  ;res_profile_VAR3@tmXBTickEndF = res_profile_VAR3@trXMaxF
  ;res_profile_VAR3@tmXBTickSpacingF = 0.5

  ;*** VAR3, DJF ***
  res_prof_VAR3_DJF = res_profile_VAR3
  plot_prof_VAR3_DJF = ignorant_xy_ensemble(wks, VAR3p_ScatMnoScat_SOMSST_DJF, pres, labels, res_prof_VAR3_DJF, opt_xy_profile)

  ;*** VAR3, JJA ***
  res_prof_VAR3_JJA = res_profile_VAR3
  plot_prof_VAR3_JJA = ignorant_xy_ensemble(wks, VAR3p_ScatMnoScat_SOMSST_JJA, pres, labels, res_prof_VAR3_JJA, opt_xy_profile)

;----------------------------------
; add polygons to block ref line
;----------------------------------
 
  res11 = True
  res11@gsFillColor = 0

  xx1 = (/-1.,-1.,1.,1.,-1./)
  yy1 = (/-3.,-2.,-2.,-3.,-3/)
  pp1a = gsn_add_polygon(wks, plot_scat_VAR1_DJF, xx1, yy1, res11)
  pp1b = gsn_add_polygon(wks, plot_scat_VAR1_JJA, xx1, yy1, res11)

  xx1 = (/-0.1,-0.1,0.1,0.1,-0.1/)
  yy1 = (/-8.,-6.,-6.,-8.,-8/)
  pp2a = gsn_add_polygon(wks, plot_scat_VAR2_DJF, xx1, yy1, res11)
  pp2b = gsn_add_polygon(wks, plot_scat_VAR2_JJA, xx1, yy1, res11)

;--------------
; plot legend
;--------------

    genres                         = True
    genres@XPosPercent             = 3                      ; move to the right
    genres@YPosPercent             = 15                      ; move to the right
    genres@ItemSpacePercent        = 6

    textres                        = True
    textres@lgLabelFontHeights     = 0.022
    textres@lgPerimOn              = False                   ; no perimeter
    textres@lgItemCount            = 2                       ; how many

    lineres                        = True
    lineres@lgLineThicknesses      = 2.5                     ; line thickness
    lineres@LineLengthPercent      = 5                       ; expressed as %, 0->100, length of line
    lineres@lgLineLabelFontHeights = 0.015                   ; font height
    lineres@lgDashIndexes          = (/0,1,0,1/)             ; line patterns
    lineres@lgLineColors           = res_scat@xyMarkerColors

    ;*** scat ***
    textres_scat = textres

    textres_scat@lgLabels               = labels_VAR1_DJF
    dum = simple_legend(wks, plot_scat_VAR1_DJF, genres, lineres, textres_scat)

    textres_scat@lgLabels               = labels_VAR2_DJF
    dum = simple_legend(wks, plot_scat_VAR2_DJF, genres, lineres, textres_scat)

    textres_scat@lgLabels               = labels_VAR1_JJA
    dum = simple_legend(wks, plot_scat_VAR1_JJA, genres, lineres, textres_scat)

    textres_scat@lgLabels               = labels_VAR2_JJA
    dum = simple_legend(wks, plot_scat_VAR2_JJA, genres, lineres, textres_scat)

    ;*** profile ***
    textres_prof = textres
    lineres_prof = lineres
    lineres_prof@lgLineColors           = res_scat@xyMarkerColors(::-1)
    ;lineres_prof@lgLineColors           = res_scat@xyMarkerColors

    genres@XPosPercent             = 75                      ; move to the right
    genres@YPosPercent             = 50                      ; move to the right

    textres_prof@lgLabels               = suffix_SSTSOM
    ;textres_prof@lgLabels               = suffix_SOMSST
    textres_prof@lgLabelFontHeights     = 0.03

    dum1 = simple_legend(wks, plot_prof_VAR1_DJF, genres, lineres_prof, textres_prof)
    dum1 = simple_legend(wks, plot_prof_VAR2_DJF, genres, lineres_prof, textres_prof)
    dum1 = simple_legend(wks, plot_prof_VAR3_DJF, genres, lineres_prof, textres_prof)
    dum1 = simple_legend(wks, plot_prof_VAR1_JJA, genres, lineres_prof, textres_prof)
    dum1 = simple_legend(wks, plot_prof_VAR2_JJA, genres, lineres_prof, textres_prof)
    dum1 = simple_legend(wks, plot_prof_VAR3_JJA, genres, lineres_prof, textres_prof)

;--------------
; plot legend
;--------------
  tres = True
  tres@txFontHeightF = 0.04

  ff_x = 0.85
  ff_y = 0.9
  xx_scat_VAR1 = ff_x*res_scat_VAR1@trXMaxF
  yy_scat_VAR1 = ff_y*res_scat_VAR1@trYMaxF
  ;text1 = gsn_add_text(wks, plot_scat_VAR1_DJF, suffix_DJF, xx_scat_VAR1, yy_scat_VAR1, tres)

;------------- 
; panel plot
;------------- 
  nplot = 10

  plot_all = new(nplot, "graphic")

  plot_all(0) = plot_scat_VAR1_DJF
  plot_all(3) = plot_scat_VAR1_JJA
  plot_all(1) = plot_scat_VAR2_DJF
  plot_all(4) = plot_scat_VAR2_JJA
  plot_all(2) = plot_prof_VAR1_DJF
  plot_all(5) = plot_prof_VAR1_JJA
  plot_all(6) = plot_prof_VAR2_DJF
  plot_all(8) = plot_prof_VAR2_JJA
  plot_all(7) = plot_prof_VAR3_DJF
  plot_all(9) = plot_prof_VAR3_JJA

  resp = True
  resp@gsnPanelFigureStringsFontHeightF= 0.013
  ;resp@gsnPanelFigureStrings= (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)"/)
  ;resp@gsnPanelFigureStrings= (/"(a) DJF","(b) DJF","(c) DJF","(d) JJA","(e) JJA","(f) JJA","(g) DJF","(h) DJF","(i) JJA","(j) JJA"/)
  resp@gsnPanelFigureStrings= (/"(a)                          DJF", \
                               "(b)                           DJF","(c)                           DJF","(d)                           JJA","(e)                           JJA","(f)                           JJA","(g)                           DJF","(h)                           DJF","(i)                           JJA","(j)                           JJA"/)                          
  resp@gsnPanelFigureStringsBackgroundFillColor = "transparent"
  resp@gsnPanelFigureStringsPerimOn = False
  ;resp@gsnPanelFigureStringsJust = "topright"
  resp@gsnPanelFigureStringsJust = "topleft"

  resp@gsnPanelXWhiteSpacePercent = 2 
  resp@gsnPanelYWhiteSpacePercent = 2

  ;column1 = 2
  ;row1 = nplot / column1
  ;gsn_panel(wks, plot_all, (/row1,column1/), resp)

  resp@gsnPanelRowSpec = True
  gsn_panel(wks, plot_all, (/3,3,2,2/), resp)


  opt11 = True
    trim_fig(wks_type, wks_name, opt11)

end
