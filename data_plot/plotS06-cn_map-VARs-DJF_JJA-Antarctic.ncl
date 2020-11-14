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

;------------------
; general settings
;------------------

;---------
; read parameters
;--------

;  VAR1_name = "CLDHGH"
  VAR1_name = "TREFHT"
  VAR2_name = "FLDS"
  VAR3_name = "TMQ"
  VAR4_name = "ICEFRAC"
  VAR_ALL = (/VAR1_name,VAR2_name,VAR3_name,VAR4_name/)

;  ;VAR1_name = "PRECC"
;  VAR1_name = "ICEFRAC"
;  VAR1_name = "CLDLOW"

;  plotname_region = "50S_90S"
;  plotname_region = "Arctic"
  plotname_region = "Antarctic"
;  ;plotname_region = "Southern_Ocean"
;
  map_region = ""
;  map_region = "50S_90S"
;
;  ;choice_landfrac = "land_only"
;  ;choice_landfrac = "non_land"
;  ;choice_landfrac = "all_grids"

;--------------
; set options
;--------------

  opt_cnmap = True
  opt_cnmap@siglvl    = 0.05
  opt_cnmap@sameScale = True
    ;opt_cnmap@sameScale = False
    ;opt_cnmap@percentile = 1
    opt_cnmap@percentile = 5
    ;opt_cnmap@percentile = 10

  wks_type = "eps"  ; or x11,pdf,ncgm
  ;wks_type = "eps"  ; or x11,pdf,ncgm

  if (isatt(opt_cnmap,"sameScale") .and. opt_cnmap@sameScale) then
    opt_cnmap@cnlevels_predefine = True
  end if 

;--------------
; set files
;--------------
  yy_string = "yy06_35"

  plotname_head = "figS08_E2000"

  casename_SOM_Scat = "c10_ensmean-cesm111-E2000_rrtmg_mc6_rtr2"
  casename_SOM_noScat = "c11_ensmean-cesm111-E2000_rrtmg_mc6_rtr2_noScat"
  casename_SST_Scat = "g01_ensmean-cesm111-F2000-mc6_rtr2_Scat"
  casename_SST_noScat = "g02_ensmean-cesm111-F2000-mc6_rtr2_noScat"

  plotname_head = plotname_head ;+"_"+choice_landfrac

  datapath = "../data/"

  filehead_SOM_Scat   = datapath+casename_SOM_Scat+".cam.h0."+yy_string+"-climo_"
  filehead_SOM_noScat   = datapath+casename_SOM_noScat+".cam.h0."+yy_string+"-climo_"
  filehead_SST_Scat   = datapath+casename_SST_Scat+".cam.h0."+yy_string+"-climo_"
  filehead_SST_noScat   = datapath+casename_SST_noScat+".cam.h0."+yy_string+"-climo_"

  ;txstring = txstring_head+", "+yy_string+" ENS, "+VAR1_name ;+", "+choice_landfrac
  txstring = ""

  suffix_Scat = "Scat"
  suffix_noScat = "noScat"
  suffix_SOM = "SOM"
  suffix_SST = "SST"
  suffix_DJF = "DJF"
  suffix_JJA = "JJA"

  ;font_title = "~F22~"
  font_title = "~F22~"

  string_VAR1_ScatMnoScat_SOM_DJF = "~F33~D~F21~SAT, "+suffix_DJF
  string_VAR2_ScatMnoScat_SOM_DJF = "~F33~D~F21~FLDS, "+suffix_DJF
  string_VAR3_ScatMnoScat_SOM_DJF = "~F33~D~F21~TCWV, "+suffix_DJF
  string_VAR4_ScatMnoScat_SOM_DJF = "~F33~D~F21~Sea ice frac, "+suffix_DJF

  string_VAR1_ScatMnoScat_SOM_JJA = "~F33~D~F21~SAT, "+suffix_JJA
  string_VAR2_ScatMnoScat_SOM_JJA = "~F33~D~F21~FLDS, "+suffix_JJA
  string_VAR3_ScatMnoScat_SOM_JJA = "~F33~D~F21~TCWV, "+suffix_JJA
  string_VAR4_ScatMnoScat_SOM_JJA = "~F33~D~F21~Sea ice frac, "+suffix_JJA

  ;suffixes = (/"  Scat - noScat, SOM", "  Scat - noScat, SST"/)

  plotname_type = "cn_map"+"_"+plotname_region
  plotname_var  = VAR1_name+"_diff"
  plotname_step = yy_string+"_DJF_JJA"
 
  ;plotname_var  = VAR1_name

  if ( map_region .ne. "" ) then
    plotname_region1 = map_region
  else
    plotname_region1 = plotname_region
  end if

  wks_name = plotname_head +"-"+ plotname_type +"-"+ plotname_region1 +"-"+ plotname_var +"-"+ plotname_step +"-"
  wks = gsn_open_wks(wks_type, wks_name)

  plotname_step = "DJF"

  infilename_SOM_Scat_DJF_region   = filehead_SOM_Scat+plotname_step+"-"+plotname_region+".nc"
  infilename_SOM_noScat_DJF_region = filehead_SOM_noScat+plotname_step+"-"+plotname_region+".nc"
  infilename_SST_Scat_DJF_region   = filehead_SST_Scat+plotname_step+"-"+plotname_region+".nc"  
  infilename_SST_noScat_DJF_region = filehead_SST_noScat+plotname_step+"-"+plotname_region+".nc"

  plotname_step = "JJA"
  infilename_SOM_Scat_JJA_region   = filehead_SOM_Scat+plotname_step+"-"+plotname_region+".nc"
  infilename_SOM_noScat_JJA_region = filehead_SOM_noScat+plotname_step+"-"+plotname_region+".nc"
  infilename_SST_Scat_JJA_region   = filehead_SST_Scat+plotname_step+"-"+plotname_region+".nc"  
  infilename_SST_noScat_JJA_region = filehead_SST_noScat+plotname_step+"-"+plotname_region+".nc"

;----------------
; read variables
;----------------
  opt111 = True
    opt111@return_type = ""


  VAR1_Scat_SOM_DJF = ignorant_convert_2Ddata(infilename_SOM_Scat_DJF_region,VAR1_name,opt111)
  VAR1_noScat_SOM_DJF = ignorant_convert_2Ddata(infilename_SOM_noScat_DJF_region,VAR1_name,opt111)
  VAR1_Scat_SOM_JJA = ignorant_convert_2Ddata(infilename_SOM_Scat_JJA_region,VAR1_name,opt111)
  VAR1_noScat_SOM_JJA = ignorant_convert_2Ddata(infilename_SOM_noScat_JJA_region,VAR1_name,opt111)

  VAR2_Scat_SOM_DJF = ignorant_convert_2Ddata(infilename_SOM_Scat_DJF_region,VAR2_name,opt111)
  VAR2_noScat_SOM_DJF = ignorant_convert_2Ddata(infilename_SOM_noScat_DJF_region,VAR2_name,opt111)
  VAR2_Scat_SOM_JJA = ignorant_convert_2Ddata(infilename_SOM_Scat_JJA_region,VAR2_name,opt111)
  VAR2_noScat_SOM_JJA = ignorant_convert_2Ddata(infilename_SOM_noScat_JJA_region,VAR2_name,opt111)

  VAR3_Scat_SOM_DJF = ignorant_convert_2Ddata(infilename_SOM_Scat_DJF_region,VAR3_name,opt111)
  VAR3_noScat_SOM_DJF = ignorant_convert_2Ddata(infilename_SOM_noScat_DJF_region,VAR3_name,opt111)
  VAR3_Scat_SOM_JJA = ignorant_convert_2Ddata(infilename_SOM_Scat_JJA_region,VAR3_name,opt111)
  VAR3_noScat_SOM_JJA = ignorant_convert_2Ddata(infilename_SOM_noScat_JJA_region,VAR3_name,opt111)

  VAR4_Scat_SOM_DJF = ignorant_convert_2Ddata(infilename_SOM_Scat_DJF_region,VAR4_name,opt111)
  VAR4_noScat_SOM_DJF = ignorant_convert_2Ddata(infilename_SOM_noScat_DJF_region,VAR4_name,opt111)
  VAR4_Scat_SOM_JJA = ignorant_convert_2Ddata(infilename_SOM_Scat_JJA_region,VAR4_name,opt111)
  VAR4_noScat_SOM_JJA = ignorant_convert_2Ddata(infilename_SOM_noScat_JJA_region,VAR4_name,opt111)

;  ;*** Scat, SOM, DJF ***
;  VARS_Scat_SOM_DJF = ignorant_get_vars_file(infilename_SOM_Scat_DJF_region, VAR_ALL, opt111)
;  delete_VarAtts(VARS_Scat_SOM_DJF, -1)
;
;  VAR1_Scat_SOM_DJF = VARS_Scat_SOM_DJF(0,:,:,:)
;  VAR2_Scat_SOM_DJF = VARS_Scat_SOM_DJF(1,:,:,:)
;  VAR3_Scat_SOM_DJF = VARS_Scat_SOM_DJF(2,:,:,:)
;  VAR4_Scat_SOM_DJF = VARS_Scat_SOM_DJF(3,:,:,:)
;
;  ;*** noScat, SOM, DJF ***
;  VARS_noScat_SOM_DJF = ignorant_get_vars_file(infilename_SOM_noScat_DJF_region, VAR_ALL, opt111)
;  delete_VarAtts(VARS_noScat_SOM_DJF, -1)
;  VAR1_noScat_SOM_DJF = VARS_noScat_SOM_DJF(0,:,:,:)
;  VAR2_noScat_SOM_DJF = VARS_noScat_SOM_DJF(1,:,:,:)
;  VAR3_noScat_SOM_DJF = VARS_noScat_SOM_DJF(2,:,:,:)
;  VAR4_noScat_SOM_DJF = VARS_noScat_SOM_DJF(3,:,:,:)
;
;  ;*** Scat, SOM, JJA ***
;  VARS_Scat_SOM_JJA = ignorant_get_vars_file(infilename_SOM_Scat_JJA_region, VAR_ALL, opt111)
;  delete_VarAtts(VARS_Scat_SOM_JJA, -1)
;  VAR1_Scat_SOM_JJA = VARS_Scat_SOM_JJA(0,:,:,:)
;  VAR2_Scat_SOM_JJA = VARS_Scat_SOM_JJA(1,:,:,:)
;  VAR3_Scat_SOM_JJA = VARS_Scat_SOM_JJA(2,:,:,:)
;  VAR4_Scat_SOM_JJA = VARS_Scat_SOM_JJA(3,:,:,:)
;
;  ;*** noScat, SOM, JJA ***
;  VARS_noScat_SOM_JJA = ignorant_get_vars_file(infilename_SOM_noScat_JJA_region, VAR_ALL, opt111)
;  delete_VarAtts(VARS_Scat_SOM_JJA, -1)
;  VAR1_noScat_SOM_JJA = VARS_noScat_SOM_JJA(0,:,:,:)
;  VAR2_noScat_SOM_JJA = VARS_noScat_SOM_JJA(1,:,:,:)
;  VAR3_noScat_SOM_JJA = VARS_noScat_SOM_JJA(2,:,:,:)
;  VAR4_noScat_SOM_JJA = VARS_noScat_SOM_JJA(3,:,:,:)

;--------
; plot
;--------
  res  = True
  res@gsnFrame = False
  res@gsnDraw = False

  res_sig = res
  ;res_sig@cnFillPatterns = (/0,4/)   ; slash
  res_sig@cnFillPatterns = (/0,17/)   ; dot

  if (plotname_region .eq. "Arctic") then
    opt_cnmap@mapview    = "polar_NH"
    res@mpMinLatF = 66.5
    ;res@mpMinLatF = 63.
  end if

  if (plotname_region .eq. "Antarctic") then
    opt_cnmap@mapview    = "polar_SH"
    res@mpMaxLatF = -66.5
    ;res@mpMinLatF = 63.
  end if

  ;res@lbLabelBarOn = False
  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnLineLabelsOn = False

  res@lbAutoManage = False
  res@lbLabelFontHeightF = 0.025

  res@gsnLeftString = "" 
  res@gsnLeftString = ""
  res@gsnCenterStringFontHeightF = res@lbLabelFontHeightF
  res@gsnRightStringFontHeightF = res@gsnCenterStringFontHeightF

  ;*** VAR1, SOM, DJF ***
  res_VAR1_ScatMnoScat_SOM= res
  if (isatt(opt_cnmap,"sameScale") .and. opt_cnmap@sameScale) then
    ignorant_sameScale( (/VAR1_Scat_SOM_DJF-VAR1_noScat_SOM_DJF, VAR1_Scat_SOM_JJA-VAR1_noScat_SOM_JJA/), res_VAR1_ScatMnoScat_SOM, opt_cnmap)
  end if

  ;*** VAR1, SOM, DJF ***
  res_VAR1_ScatMnoScat_SOM_DJF = res_VAR1_ScatMnoScat_SOM
  res_VAR1_ScatMnoScat_SOM_DJF@tiMainString = string_VAR1_ScatMnoScat_SOM_DJF
  plot_VAR1_ScatMnoScat_SOM_DJF = gsn_csm_contour_map_siglvl(wks, VAR1_Scat_SOM_DJF, VAR1_noScat_SOM_DJF, res_VAR1_ScatMnoScat_SOM_DJF, res_sig,opt_cnmap)

  ;*** VAR1, SOM, JJA ***
  res_VAR1_ScatMnoScat_SOM_JJA = res_VAR1_ScatMnoScat_SOM
  res_VAR1_ScatMnoScat_SOM_JJA@tiMainString = string_VAR1_ScatMnoScat_SOM_JJA 
  plot_VAR1_ScatMnoScat_SOM_JJA = gsn_csm_contour_map_siglvl(wks, VAR1_Scat_SOM_JJA, VAR1_noScat_SOM_JJA, res_VAR1_ScatMnoScat_SOM_JJA, res_sig,opt_cnmap)

  ;*** VAR2, SOM, DJF ***
  res_VAR2_ScatMnoScat_SOM= res
  if (isatt(opt_cnmap,"sameScale") .and. opt_cnmap@sameScale) then
    ignorant_sameScale( (/VAR2_Scat_SOM_DJF-VAR2_noScat_SOM_DJF, VAR2_Scat_SOM_JJA-VAR2_noScat_SOM_JJA/), res_VAR2_ScatMnoScat_SOM, opt_cnmap)
  end if

  ;*** VAR2, SOM, DJF ***
  res_VAR2_ScatMnoScat_SOM_DJF = res_VAR2_ScatMnoScat_SOM
  res_VAR2_ScatMnoScat_SOM_DJF@tiMainString = string_VAR2_ScatMnoScat_SOM_DJF 
  plot_VAR2_ScatMnoScat_SOM_DJF = gsn_csm_contour_map_siglvl(wks, VAR2_Scat_SOM_DJF, VAR2_noScat_SOM_DJF, res_VAR2_ScatMnoScat_SOM_DJF, res_sig,opt_cnmap)

  ;*** VAR2, SOM, JJA ***
  res_VAR2_ScatMnoScat_SOM_JJA = res_VAR2_ScatMnoScat_SOM
  res_VAR2_ScatMnoScat_SOM_JJA@tiMainString = string_VAR2_ScatMnoScat_SOM_JJA 
  plot_VAR2_ScatMnoScat_SOM_JJA = gsn_csm_contour_map_siglvl(wks, VAR2_Scat_SOM_JJA, VAR2_noScat_SOM_JJA, res_VAR2_ScatMnoScat_SOM_JJA, res_sig,opt_cnmap)

  ;*** VAR3, SOM, DJF ***
  res_VAR3_ScatMnoScat_SOM= res
  if (isatt(opt_cnmap,"sameScale") .and. opt_cnmap@sameScale) then
    ignorant_sameScale( (/VAR3_Scat_SOM_DJF-VAR3_noScat_SOM_DJF, VAR3_Scat_SOM_JJA-VAR3_noScat_SOM_JJA/), res_VAR3_ScatMnoScat_SOM, opt_cnmap)
  end if

  ;*** VAR3, SOM, DJF ***
  res_VAR3_ScatMnoScat_SOM_DJF = res_VAR3_ScatMnoScat_SOM
  res_VAR3_ScatMnoScat_SOM_DJF@tiMainString = string_VAR3_ScatMnoScat_SOM_DJF 
  plot_VAR3_ScatMnoScat_SOM_DJF = gsn_csm_contour_map_siglvl(wks, VAR3_Scat_SOM_DJF, VAR3_noScat_SOM_DJF, res_VAR3_ScatMnoScat_SOM_DJF, res_sig,opt_cnmap)

  ;*** VAR3, SOM, JJA ***
  res_VAR3_ScatMnoScat_SOM_JJA = res_VAR3_ScatMnoScat_SOM
  res_VAR3_ScatMnoScat_SOM_JJA@tiMainString = string_VAR3_ScatMnoScat_SOM_JJA 
  plot_VAR3_ScatMnoScat_SOM_JJA = gsn_csm_contour_map_siglvl(wks, VAR3_Scat_SOM_JJA, VAR3_noScat_SOM_JJA, res_VAR3_ScatMnoScat_SOM_JJA, res_sig,opt_cnmap)

  ;*** VAR4, SOM, DJF ***
  res_VAR4_ScatMnoScat_SOM= res
  if (isatt(opt_cnmap,"sameScale") .and. opt_cnmap@sameScale) then
    ignorant_sameScale( (/VAR4_Scat_SOM_DJF-VAR4_noScat_SOM_DJF, VAR4_Scat_SOM_JJA-VAR4_noScat_SOM_JJA/), res_VAR4_ScatMnoScat_SOM, opt_cnmap)
  end if

  ;*** VAR4, SOM, DJF ***
  res_VAR4_ScatMnoScat_SOM_DJF = res_VAR4_ScatMnoScat_SOM
  res_VAR4_ScatMnoScat_SOM_DJF@tiMainString = string_VAR4_ScatMnoScat_SOM_DJF 
  plot_VAR4_ScatMnoScat_SOM_DJF = gsn_csm_contour_map_siglvl(wks, VAR4_Scat_SOM_DJF, VAR4_noScat_SOM_DJF, res_VAR4_ScatMnoScat_SOM_DJF, res_sig,opt_cnmap)

  ;*** VAR4, SOM, JJA ***
  res_VAR4_ScatMnoScat_SOM_JJA = res_VAR4_ScatMnoScat_SOM
  res_VAR4_ScatMnoScat_SOM_JJA@tiMainString = string_VAR4_ScatMnoScat_SOM_JJA 
  plot_VAR4_ScatMnoScat_SOM_JJA = gsn_csm_contour_map_siglvl(wks, VAR4_Scat_SOM_JJA, VAR4_noScat_SOM_JJA, res_VAR4_ScatMnoScat_SOM_JJA, res_sig,opt_cnmap)


;--------
; panel
;--------
  nplot = 8

  plot_all = new(nplot,graphic)

  plot_all(0) = plot_VAR1_ScatMnoScat_SOM_DJF
  plot_all(1) = plot_VAR2_ScatMnoScat_SOM_DJF
  plot_all(2) = plot_VAR3_ScatMnoScat_SOM_DJF
  plot_all(3) = plot_VAR4_ScatMnoScat_SOM_DJF
  plot_all(4) = plot_VAR1_ScatMnoScat_SOM_JJA
  plot_all(5) = plot_VAR2_ScatMnoScat_SOM_JJA
  plot_all(6) = plot_VAR3_ScatMnoScat_SOM_JJA
  plot_all(7) = plot_VAR4_ScatMnoScat_SOM_JJA

  resp = True    ; panel res
  resp@gsnMaximize = True
  resp@txString = txstring
  ;resp@gsnPanelLabelBar = True
  resp@gsnPanelFigureStrings= (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"/) 	; add figure string
  ;resp@gsnPanelFigureStringsBackgroundFillColor = "transparent"
  resp@gsnPanelFigureStringsJust = "bottomleft"
  resp@gsnPanelYWhiteSpacePercent = 5

  row1 = 2
  column1 = nplot / row1
  gsn_panel(wks, plot_all, (/row1,column1/), resp)

  opt1 = True
  trim_fig(wks_type, wks_name, opt1)

end
