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

  plotname_region = "Arctic"

  wks_name = "figS04-scattered-FLDS_TREFHT_areamean-4members-SOM_SST-DJF_JJA-"+plotname_region
  wks = gsn_open_wks(wks_type, wks_name)

  outfile_scattered = "data_plotS04-scattered-FLDS_TREFHT_areamean-SOM_SST-DJF_JJA-"+plotname_region+".nc"

  VAR1x_name = "FLDS"
  VAR1y_name = "TREFHT"

  xstring = "~F33~D~F21~FLDS (W m~S~-2~N~)"
  ystring = "~F33~D~F21~SAT (K)"

  suffix_DJF = "DJF"
  suffix_JJA = "JJA"

  suffix_SOM = "SOM"
  suffix_SST = "SST"

  suffix_en00 = "Member 1"
  suffix_en01 = "Member 2"
  suffix_en02 = "Member 3"
  suffix_en03 = "Member 4"

;------------------
; open files and read variables
;------------------

;**********************
  out_scattered = addfile(outfile_scattered,"r")

  ;*** en00, VAR1, DJF ***
  VAR1x_en00_SOM_ScatMnoScat_DJF = out_scattered->VAR1x_en00_SOM_ScatMnoScat_DJF
  VAR1y_en00_SOM_ScatMnoScat_DJF = out_scattered->VAR1y_en00_SOM_ScatMnoScat_DJF
  VAR1x_en00_SST_ScatMnoScat_DJF = out_scattered->VAR1x_en00_SST_ScatMnoScat_DJF
  VAR1y_en00_SST_ScatMnoScat_DJF = out_scattered->VAR1y_en00_SST_ScatMnoScat_DJF

  ;*** en00, VAR1, JJA ***
  VAR1x_en00_SOM_ScatMnoScat_JJA = out_scattered->VAR1x_en00_SOM_ScatMnoScat_JJA
  VAR1y_en00_SOM_ScatMnoScat_JJA = out_scattered->VAR1y_en00_SOM_ScatMnoScat_JJA
  VAR1x_en00_SST_ScatMnoScat_JJA = out_scattered->VAR1x_en00_SST_ScatMnoScat_JJA
  VAR1y_en00_SST_ScatMnoScat_JJA = out_scattered->VAR1y_en00_SST_ScatMnoScat_JJA

  ;*** en01, VAR1, DJF ***
  VAR1x_en01_SOM_ScatMnoScat_DJF = out_scattered->VAR1x_en01_SOM_ScatMnoScat_DJF
  VAR1y_en01_SOM_ScatMnoScat_DJF = out_scattered->VAR1y_en01_SOM_ScatMnoScat_DJF
  VAR1x_en01_SST_ScatMnoScat_DJF = out_scattered->VAR1x_en01_SST_ScatMnoScat_DJF
  VAR1y_en01_SST_ScatMnoScat_DJF = out_scattered->VAR1y_en01_SST_ScatMnoScat_DJF

  ;*** en01, VAR1, JJA ***
  VAR1x_en01_SOM_ScatMnoScat_JJA = out_scattered->VAR1x_en01_SOM_ScatMnoScat_JJA
  VAR1y_en01_SOM_ScatMnoScat_JJA = out_scattered->VAR1y_en01_SOM_ScatMnoScat_JJA
  VAR1x_en01_SST_ScatMnoScat_JJA = out_scattered->VAR1x_en01_SST_ScatMnoScat_JJA
  VAR1y_en01_SST_ScatMnoScat_JJA = out_scattered->VAR1y_en01_SST_ScatMnoScat_JJA

  ;*** en02, VAR1, DJF ***
  VAR1x_en02_SOM_ScatMnoScat_DJF = out_scattered->VAR1x_en02_SOM_ScatMnoScat_DJF
  VAR1y_en02_SOM_ScatMnoScat_DJF = out_scattered->VAR1y_en02_SOM_ScatMnoScat_DJF
  VAR1x_en02_SST_ScatMnoScat_DJF = out_scattered->VAR1x_en02_SST_ScatMnoScat_DJF
  VAR1y_en02_SST_ScatMnoScat_DJF = out_scattered->VAR1y_en02_SST_ScatMnoScat_DJF

  ;*** en02, VAR1, JJA ***
  VAR1x_en02_SOM_ScatMnoScat_JJA = out_scattered->VAR1x_en02_SOM_ScatMnoScat_JJA
  VAR1y_en02_SOM_ScatMnoScat_JJA = out_scattered->VAR1y_en02_SOM_ScatMnoScat_JJA
  VAR1x_en02_SST_ScatMnoScat_JJA = out_scattered->VAR1x_en02_SST_ScatMnoScat_JJA
  VAR1y_en02_SST_ScatMnoScat_JJA = out_scattered->VAR1y_en02_SST_ScatMnoScat_JJA

  ;*** en03, VAR1, DJF ***
  VAR1x_en03_SOM_ScatMnoScat_DJF = out_scattered->VAR1x_en03_SOM_ScatMnoScat_DJF
  VAR1y_en03_SOM_ScatMnoScat_DJF = out_scattered->VAR1y_en03_SOM_ScatMnoScat_DJF
  VAR1x_en03_SST_ScatMnoScat_DJF = out_scattered->VAR1x_en03_SST_ScatMnoScat_DJF
  VAR1y_en03_SST_ScatMnoScat_DJF = out_scattered->VAR1y_en03_SST_ScatMnoScat_DJF

  ;*** en03, VAR1, JJA ***
  VAR1x_en03_SOM_ScatMnoScat_JJA = out_scattered->VAR1x_en03_SOM_ScatMnoScat_JJA
  VAR1y_en03_SOM_ScatMnoScat_JJA = out_scattered->VAR1y_en03_SOM_ScatMnoScat_JJA
  VAR1x_en03_SST_ScatMnoScat_JJA = out_scattered->VAR1x_en03_SST_ScatMnoScat_JJA
  VAR1y_en03_SST_ScatMnoScat_JJA = out_scattered->VAR1y_en03_SST_ScatMnoScat_JJA

;----------
; regline
;----------

  opt_regline = True
    string_input = ""
    opt_regline@printout = False

  string_input = ""

;****************
; VAR1, en00
;****************

  ;*** VAR1, SOM, DJF ***
  regline_en00_SOM_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_en00_SOM_ScatMnoScat_DJF,VAR1y_en00_SOM_ScatMnoScat_DJF,string_input,opt_regline)
  string_en00_SOM_VAR1_ScatMnoScat_DJF = regline_en00_SOM_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_en00_SOM_VAR1_ScatMnoScat_DJF)
  ;print(regline_en00_SOM_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_en00_SOM_VAR1_DJF = regline_en00_SOM_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_en00_SOM_VAR1_DJF = regline_en00_SOM_VAR1_ScatMnoScat_DJF@line_yy

  ;*** VAR1, SST, DJF ***
  regline_en00_SST_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_en00_SST_ScatMnoScat_DJF,VAR1y_en00_SST_ScatMnoScat_DJF,string_input,opt_regline)
  string_en00_SST_VAR1_ScatMnoScat_DJF = regline_en00_SST_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_en00_SST_VAR1_ScatMnoScat_DJF)
  ;print(regline_en00_SST_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_en00_SST_VAR1_DJF = regline_en00_SST_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_en00_SST_VAR1_DJF = regline_en00_SST_VAR1_ScatMnoScat_DJF@line_yy

  VAR1_en00_plot_xx_DJF = (/VAR1x_en00_SOM_ScatMnoScat_DJF, VAR1x_en00_SST_ScatMnoScat_DJF/) 
  VAR1_en00_plot_yy_DJF = (/VAR1y_en00_SOM_ScatMnoScat_DJF, VAR1y_en00_SST_ScatMnoScat_DJF/)

  labels_en00_VAR1_DJF = (/suffix_SOM+": "+regline_en00_SOM_VAR1_ScatMnoScat_DJF@string_fit, \
                           suffix_SST+": "+regline_en00_SST_VAR1_ScatMnoScat_DJF@string_fit/)

  ;*** VAR1, SOM, JJA ***
  regline_en00_SOM_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_en00_SOM_ScatMnoScat_JJA,VAR1y_en00_SOM_ScatMnoScat_JJA,string_input,opt_regline)
  string_en00_SOM_VAR1_ScatMnoScat_JJA = regline_en00_SOM_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_en00_SOM_VAR1_ScatMnoScat_JJA)
  ;print(regline_en00_SOM_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_en00_SOM_VAR1_JJA = regline_en00_SOM_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_en00_SOM_VAR1_JJA = regline_en00_SOM_VAR1_ScatMnoScat_JJA@line_yy

  ;*** VAR1, SST, JJA ***
  regline_en00_SST_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_en00_SST_ScatMnoScat_JJA,VAR1y_en00_SST_ScatMnoScat_JJA,string_input,opt_regline)
  string_en00_SST_VAR1_ScatMnoScat_JJA = regline_en00_SST_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_en00_SST_VAR1_ScatMnoScat_JJA)
  ;print(regline_en00_SST_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_en00_SST_VAR1_JJA = regline_en00_SST_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_en00_SST_VAR1_JJA = regline_en00_SST_VAR1_ScatMnoScat_JJA@line_yy

  VAR1_en00_plot_xx_JJA = (/VAR1x_en00_SOM_ScatMnoScat_JJA, VAR1x_en00_SST_ScatMnoScat_JJA/) 
  VAR1_en00_plot_yy_JJA = (/VAR1y_en00_SOM_ScatMnoScat_JJA, VAR1y_en00_SST_ScatMnoScat_JJA/)

  labels_en00_VAR1_JJA = (/suffix_SOM+": "+regline_en00_SOM_VAR1_ScatMnoScat_JJA@string_fit, \
                           suffix_SST+": "+regline_en00_SST_VAR1_ScatMnoScat_JJA@string_fit/)

;****************
; VAR1, en01
;****************

  ;*** VAR1, SOM, DJF ***
  regline_en01_SOM_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_en01_SOM_ScatMnoScat_DJF,VAR1y_en01_SOM_ScatMnoScat_DJF,string_input,opt_regline)
  string_en01_SOM_VAR1_ScatMnoScat_DJF = regline_en01_SOM_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_en01_SOM_VAR1_ScatMnoScat_DJF)
  ;print(regline_en01_SOM_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_en01_SOM_VAR1_DJF = regline_en01_SOM_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_en01_SOM_VAR1_DJF = regline_en01_SOM_VAR1_ScatMnoScat_DJF@line_yy

  ;*** VAR1, SST, DJF ***
  regline_en01_SST_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_en01_SST_ScatMnoScat_DJF,VAR1y_en01_SST_ScatMnoScat_DJF,string_input,opt_regline)
  string_en01_SST_VAR1_ScatMnoScat_DJF = regline_en01_SST_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_en01_SST_VAR1_ScatMnoScat_DJF)
  ;print(regline_en01_SST_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_en01_SST_VAR1_DJF = regline_en01_SST_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_en01_SST_VAR1_DJF = regline_en01_SST_VAR1_ScatMnoScat_DJF@line_yy

  VAR1_en01_plot_xx_DJF = (/VAR1x_en01_SOM_ScatMnoScat_DJF, VAR1x_en01_SST_ScatMnoScat_DJF/) 
  VAR1_en01_plot_yy_DJF = (/VAR1y_en01_SOM_ScatMnoScat_DJF, VAR1y_en01_SST_ScatMnoScat_DJF/)

  labels_en01_VAR1_DJF = (/suffix_SOM+": "+regline_en01_SOM_VAR1_ScatMnoScat_DJF@string_fit, \
                           suffix_SST+": "+regline_en01_SST_VAR1_ScatMnoScat_DJF@string_fit/)

  ;*** VAR1, SOM, JJA ***
  regline_en01_SOM_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_en01_SOM_ScatMnoScat_JJA,VAR1y_en01_SOM_ScatMnoScat_JJA,string_input,opt_regline)
  string_en01_SOM_VAR1_ScatMnoScat_JJA = regline_en01_SOM_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_en01_SOM_VAR1_ScatMnoScat_JJA)
  ;print(regline_en01_SOM_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_en01_SOM_VAR1_JJA = regline_en01_SOM_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_en01_SOM_VAR1_JJA = regline_en01_SOM_VAR1_ScatMnoScat_JJA@line_yy

  ;*** VAR1, SST, JJA ***
  regline_en01_SST_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_en01_SST_ScatMnoScat_JJA,VAR1y_en01_SST_ScatMnoScat_JJA,string_input,opt_regline)
  string_en01_SST_VAR1_ScatMnoScat_JJA = regline_en01_SST_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_en01_SST_VAR1_ScatMnoScat_JJA)
  ;print(regline_en01_SST_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_en01_SST_VAR1_JJA = regline_en01_SST_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_en01_SST_VAR1_JJA = regline_en01_SST_VAR1_ScatMnoScat_JJA@line_yy

  VAR1_en01_plot_xx_JJA = (/VAR1x_en01_SOM_ScatMnoScat_JJA, VAR1x_en01_SST_ScatMnoScat_JJA/) 
  VAR1_en01_plot_yy_JJA = (/VAR1y_en01_SOM_ScatMnoScat_JJA, VAR1y_en01_SST_ScatMnoScat_JJA/)

  labels_en01_VAR1_JJA = (/suffix_SOM+": "+regline_en01_SOM_VAR1_ScatMnoScat_JJA@string_fit, \
                           suffix_SST+": "+regline_en01_SST_VAR1_ScatMnoScat_JJA@string_fit/)

;****************
; VAR1, en02
;****************

  ;*** VAR1, SOM, DJF ***
  regline_en02_SOM_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_en02_SOM_ScatMnoScat_DJF,VAR1y_en02_SOM_ScatMnoScat_DJF,string_input,opt_regline)
  string_en02_SOM_VAR1_ScatMnoScat_DJF = regline_en02_SOM_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_en02_SOM_VAR1_ScatMnoScat_DJF)
  ;print(regline_en02_SOM_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_en02_SOM_VAR1_DJF = regline_en02_SOM_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_en02_SOM_VAR1_DJF = regline_en02_SOM_VAR1_ScatMnoScat_DJF@line_yy

  ;*** VAR1, SST, DJF ***
  regline_en02_SST_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_en02_SST_ScatMnoScat_DJF,VAR1y_en02_SST_ScatMnoScat_DJF,string_input,opt_regline)
  string_en02_SST_VAR1_ScatMnoScat_DJF = regline_en02_SST_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_en02_SST_VAR1_ScatMnoScat_DJF)
  ;print(regline_en02_SST_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_en02_SST_VAR1_DJF = regline_en02_SST_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_en02_SST_VAR1_DJF = regline_en02_SST_VAR1_ScatMnoScat_DJF@line_yy

  VAR1_en02_plot_xx_DJF = (/VAR1x_en02_SOM_ScatMnoScat_DJF, VAR1x_en02_SST_ScatMnoScat_DJF/) 
  VAR1_en02_plot_yy_DJF = (/VAR1y_en02_SOM_ScatMnoScat_DJF, VAR1y_en02_SST_ScatMnoScat_DJF/)

  labels_en02_VAR1_DJF = (/suffix_SOM+": "+regline_en02_SOM_VAR1_ScatMnoScat_DJF@string_fit, \
                           suffix_SST+": "+regline_en02_SST_VAR1_ScatMnoScat_DJF@string_fit/)

  ;*** VAR1, SOM, JJA ***
  regline_en02_SOM_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_en02_SOM_ScatMnoScat_JJA,VAR1y_en02_SOM_ScatMnoScat_JJA,string_input,opt_regline)
  string_en02_SOM_VAR1_ScatMnoScat_JJA = regline_en02_SOM_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_en02_SOM_VAR1_ScatMnoScat_JJA)
  ;print(regline_en02_SOM_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_en02_SOM_VAR1_JJA = regline_en02_SOM_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_en02_SOM_VAR1_JJA = regline_en02_SOM_VAR1_ScatMnoScat_JJA@line_yy

  ;*** VAR1, SST, JJA ***
  regline_en02_SST_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_en02_SST_ScatMnoScat_JJA,VAR1y_en02_SST_ScatMnoScat_JJA,string_input,opt_regline)
  string_en02_SST_VAR1_ScatMnoScat_JJA = regline_en02_SST_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_en02_SST_VAR1_ScatMnoScat_JJA)
  ;print(regline_en02_SST_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_en02_SST_VAR1_JJA = regline_en02_SST_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_en02_SST_VAR1_JJA = regline_en02_SST_VAR1_ScatMnoScat_JJA@line_yy

  VAR1_en02_plot_xx_JJA = (/VAR1x_en02_SOM_ScatMnoScat_JJA, VAR1x_en02_SST_ScatMnoScat_JJA/) 
  VAR1_en02_plot_yy_JJA = (/VAR1y_en02_SOM_ScatMnoScat_JJA, VAR1y_en02_SST_ScatMnoScat_JJA/)

  labels_en02_VAR1_JJA = (/suffix_SOM+": "+regline_en02_SOM_VAR1_ScatMnoScat_JJA@string_fit, \
                           suffix_SST+": "+regline_en02_SST_VAR1_ScatMnoScat_JJA@string_fit/)

;****************
; VAR1, en03
;****************

  ;*** VAR1, SOM, DJF ***
  regline_en03_SOM_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_en03_SOM_ScatMnoScat_DJF,VAR1y_en03_SOM_ScatMnoScat_DJF,string_input,opt_regline)
  string_en03_SOM_VAR1_ScatMnoScat_DJF = regline_en03_SOM_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_en03_SOM_VAR1_ScatMnoScat_DJF)
  ;print(regline_en03_SOM_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_en03_SOM_VAR1_DJF = regline_en03_SOM_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_en03_SOM_VAR1_DJF = regline_en03_SOM_VAR1_ScatMnoScat_DJF@line_yy

  ;*** VAR1, SST, DJF ***
  regline_en03_SST_VAR1_ScatMnoScat_DJF = regline_print(VAR1x_en03_SST_ScatMnoScat_DJF,VAR1y_en03_SST_ScatMnoScat_DJF,string_input,opt_regline)
  string_en03_SST_VAR1_ScatMnoScat_DJF = regline_en03_SST_VAR1_ScatMnoScat_DJF@string_fit
  ;print(regline_en03_SST_VAR1_ScatMnoScat_DJF)
  ;print(regline_en03_SST_VAR1_ScatMnoScat_DJF@string_fit)
  line_xx_en03_SST_VAR1_DJF = regline_en03_SST_VAR1_ScatMnoScat_DJF@line_xx
  line_yy_en03_SST_VAR1_DJF = regline_en03_SST_VAR1_ScatMnoScat_DJF@line_yy

  VAR1_en03_plot_xx_DJF = (/VAR1x_en03_SOM_ScatMnoScat_DJF, VAR1x_en03_SST_ScatMnoScat_DJF/) 
  VAR1_en03_plot_yy_DJF = (/VAR1y_en03_SOM_ScatMnoScat_DJF, VAR1y_en03_SST_ScatMnoScat_DJF/)

  labels_en03_VAR1_DJF = (/suffix_SOM+": "+regline_en03_SOM_VAR1_ScatMnoScat_DJF@string_fit, \
                           suffix_SST+": "+regline_en03_SST_VAR1_ScatMnoScat_DJF@string_fit/)

  ;*** VAR1, SOM, JJA ***
  regline_en03_SOM_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_en03_SOM_ScatMnoScat_JJA,VAR1y_en03_SOM_ScatMnoScat_JJA,string_input,opt_regline)
  string_en03_SOM_VAR1_ScatMnoScat_JJA = regline_en03_SOM_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_en03_SOM_VAR1_ScatMnoScat_JJA)
  ;print(regline_en03_SOM_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_en03_SOM_VAR1_JJA = regline_en03_SOM_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_en03_SOM_VAR1_JJA = regline_en03_SOM_VAR1_ScatMnoScat_JJA@line_yy

  ;*** VAR1, SST, JJA ***
  regline_en03_SST_VAR1_ScatMnoScat_JJA = regline_print(VAR1x_en03_SST_ScatMnoScat_JJA,VAR1y_en03_SST_ScatMnoScat_JJA,string_input,opt_regline)
  string_en03_SST_VAR1_ScatMnoScat_JJA = regline_en03_SST_VAR1_ScatMnoScat_JJA@string_fit
  ;print(regline_en03_SST_VAR1_ScatMnoScat_JJA)
  ;print(regline_en03_SST_VAR1_ScatMnoScat_JJA@string_fit)
  line_xx_en03_SST_VAR1_JJA = regline_en03_SST_VAR1_ScatMnoScat_JJA@line_xx
  line_yy_en03_SST_VAR1_JJA = regline_en03_SST_VAR1_ScatMnoScat_JJA@line_yy

  VAR1_en03_plot_xx_JJA = (/VAR1x_en03_SOM_ScatMnoScat_JJA, VAR1x_en03_SST_ScatMnoScat_JJA/) 
  VAR1_en03_plot_yy_JJA = (/VAR1y_en03_SOM_ScatMnoScat_JJA, VAR1y_en03_SST_ScatMnoScat_JJA/)

  labels_en03_VAR1_JJA = (/suffix_SOM+": "+regline_en03_SOM_VAR1_ScatMnoScat_JJA@string_fit, \
                           suffix_SST+": "+regline_en03_SST_VAR1_ScatMnoScat_JJA@string_fit/)


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
  res_scat_VAR1@tiXAxisString = xstring
  res_scat_VAR1@tiYAxisString = ystring
  res_scat_VAR1@trXMinF = -15.
  res_scat_VAR1@trXMaxF = 20.
  res_scat_VAR1@trYMinF = -4.
  res_scat_VAR1@trYMaxF = 6.

  ;*** en00, VAR1, DJF ***
  res_scat_VAR1_DJF = res_scat_VAR1
  res_scat_VAR1_DJF@tiMainString = suffix_en00+", "+suffix_DJF 
  plot_scat_en00_VAR1_DJF = ignorant_xy(wks, VAR1_en00_plot_xx_DJF, VAR1_en00_plot_yy_DJF, labels_en00_VAR1_DJF, res_scat_VAR1_DJF, opt_xy_scat)

  line_en00_SOM_VAR1_DJF = gsn_add_polyline(wks, plot_scat_en00_VAR1_DJF, line_xx_en00_SOM_VAR1_DJF, line_yy_en00_SOM_VAR1_DJF, lnres_SOM)
  line_en00_SST_VAR1_DJF = gsn_add_polyline(wks, plot_scat_en00_VAR1_DJF, line_xx_en00_SST_VAR1_DJF, line_yy_en00_SST_VAR1_DJF, lnres_SST)

  ;*** en00, VAR1, JJA ***
  res_scat_VAR1_JJA = res_scat_VAR1
  res_scat_VAR1_JJA@tiMainString = suffix_en00+", "+suffix_JJA 
  plot_scat_en00_VAR1_JJA = ignorant_xy(wks, VAR1_en00_plot_xx_JJA, VAR1_en00_plot_yy_JJA, labels_en00_VAR1_JJA, res_scat_VAR1_JJA, opt_xy_scat)

  line_en00_SOM_VAR1_JJA = gsn_add_polyline(wks, plot_scat_en00_VAR1_JJA, line_xx_en00_SOM_VAR1_JJA, line_yy_en00_SOM_VAR1_JJA, lnres_SOM)
  line_en00_SST_VAR1_JJA = gsn_add_polyline(wks, plot_scat_en00_VAR1_JJA, line_xx_en00_SST_VAR1_JJA, line_yy_en00_SST_VAR1_JJA, lnres_SST)

  ;*** en01, VAR1, DJF ***
  res_scat_VAR1_DJF = res_scat_VAR1
  res_scat_VAR1_DJF@tiMainString = suffix_en01+", "+suffix_DJF 
  plot_scat_en01_VAR1_DJF = ignorant_xy(wks, VAR1_en01_plot_xx_DJF, VAR1_en01_plot_yy_DJF, labels_en01_VAR1_DJF, res_scat_VAR1_DJF, opt_xy_scat)

  line_en01_SOM_VAR1_DJF = gsn_add_polyline(wks, plot_scat_en01_VAR1_DJF, line_xx_en01_SOM_VAR1_DJF, line_yy_en01_SOM_VAR1_DJF, lnres_SOM)
  line_en01_SST_VAR1_DJF = gsn_add_polyline(wks, plot_scat_en01_VAR1_DJF, line_xx_en01_SST_VAR1_DJF, line_yy_en01_SST_VAR1_DJF, lnres_SST)

  ;*** en01, VAR1, JJA ***
  res_scat_VAR1_JJA = res_scat_VAR1
  res_scat_VAR1_JJA@tiMainString = suffix_en01+", "+suffix_JJA 
  plot_scat_en01_VAR1_JJA = ignorant_xy(wks, VAR1_en01_plot_xx_JJA, VAR1_en01_plot_yy_JJA, labels_en01_VAR1_JJA, res_scat_VAR1_JJA, opt_xy_scat)

  line_en01_SOM_VAR1_JJA = gsn_add_polyline(wks, plot_scat_en01_VAR1_JJA, line_xx_en01_SOM_VAR1_JJA, line_yy_en01_SOM_VAR1_JJA, lnres_SOM)
  line_en01_SST_VAR1_JJA = gsn_add_polyline(wks, plot_scat_en01_VAR1_JJA, line_xx_en01_SST_VAR1_JJA, line_yy_en01_SST_VAR1_JJA, lnres_SST)

  ;*** en02, VAR1, DJF ***
  res_scat_VAR1_DJF = res_scat_VAR1
  res_scat_VAR1_DJF@tiMainString = suffix_en02+", "+suffix_DJF 
  plot_scat_en02_VAR1_DJF = ignorant_xy(wks, VAR1_en02_plot_xx_DJF, VAR1_en02_plot_yy_DJF, labels_en02_VAR1_DJF, res_scat_VAR1_DJF, opt_xy_scat)

  line_en02_SOM_VAR1_DJF = gsn_add_polyline(wks, plot_scat_en02_VAR1_DJF, line_xx_en02_SOM_VAR1_DJF, line_yy_en02_SOM_VAR1_DJF, lnres_SOM)
  line_en02_SST_VAR1_DJF = gsn_add_polyline(wks, plot_scat_en02_VAR1_DJF, line_xx_en02_SST_VAR1_DJF, line_yy_en02_SST_VAR1_DJF, lnres_SST)

  ;*** en02, VAR1, JJA ***
  res_scat_VAR1_JJA = res_scat_VAR1
  res_scat_VAR1_JJA@tiMainString = suffix_en02+", "+suffix_JJA 
  plot_scat_en02_VAR1_JJA = ignorant_xy(wks, VAR1_en02_plot_xx_JJA, VAR1_en02_plot_yy_JJA, labels_en02_VAR1_JJA, res_scat_VAR1_JJA, opt_xy_scat)

  line_en02_SOM_VAR1_JJA = gsn_add_polyline(wks, plot_scat_en02_VAR1_JJA, line_xx_en02_SOM_VAR1_JJA, line_yy_en02_SOM_VAR1_JJA, lnres_SOM)
  line_en02_SST_VAR1_JJA = gsn_add_polyline(wks, plot_scat_en02_VAR1_JJA, line_xx_en02_SST_VAR1_JJA, line_yy_en02_SST_VAR1_JJA, lnres_SST)

  ;*** en03, VAR1, DJF ***
  res_scat_VAR1_DJF = res_scat_VAR1
  res_scat_VAR1_DJF@tiMainString = suffix_en03+", "+suffix_DJF 
  plot_scat_en03_VAR1_DJF = ignorant_xy(wks, VAR1_en03_plot_xx_DJF, VAR1_en03_plot_yy_DJF, labels_en03_VAR1_DJF, res_scat_VAR1_DJF, opt_xy_scat)

  line_en03_SOM_VAR1_DJF = gsn_add_polyline(wks, plot_scat_en03_VAR1_DJF, line_xx_en03_SOM_VAR1_DJF, line_yy_en03_SOM_VAR1_DJF, lnres_SOM)
  line_en03_SST_VAR1_DJF = gsn_add_polyline(wks, plot_scat_en03_VAR1_DJF, line_xx_en03_SST_VAR1_DJF, line_yy_en03_SST_VAR1_DJF, lnres_SST)

  ;*** en03, VAR1, JJA ***
  res_scat_VAR1_JJA = res_scat_VAR1
  res_scat_VAR1_JJA@tiMainString = suffix_en03+", "+suffix_JJA 
  plot_scat_en03_VAR1_JJA = ignorant_xy(wks, VAR1_en03_plot_xx_JJA, VAR1_en03_plot_yy_JJA, labels_en03_VAR1_JJA, res_scat_VAR1_JJA, opt_xy_scat)

  line_en03_SOM_VAR1_JJA = gsn_add_polyline(wks, plot_scat_en03_VAR1_JJA, line_xx_en03_SOM_VAR1_JJA, line_yy_en03_SOM_VAR1_JJA, lnres_SOM)
  line_en03_SST_VAR1_JJA = gsn_add_polyline(wks, plot_scat_en03_VAR1_JJA, line_xx_en03_SST_VAR1_JJA, line_yy_en03_SST_VAR1_JJA, lnres_SST)

;-------------
; panel plot
;-------------
  nplot = 8

  plot_all = new(nplot, "graphic")

  plot_all(0) = plot_scat_en00_VAR1_DJF
  plot_all(1) = plot_scat_en00_VAR1_JJA
  plot_all(2) = plot_scat_en01_VAR1_DJF
  plot_all(3) = plot_scat_en01_VAR1_JJA
  plot_all(4) = plot_scat_en02_VAR1_DJF
  plot_all(5) = plot_scat_en02_VAR1_JJA
  plot_all(6) = plot_scat_en03_VAR1_DJF
  plot_all(7) = plot_scat_en03_VAR1_JJA

  resp = True
  column1 = 2
  row1 = nplot / column1
  gsn_panel(wks, plot_all, (/row1,column1/), resp)

  opt1 = True
  trim_fig(wks_type, wks_name, opt1)

end
