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


  plotname_region = "Arctic"

  VAR1x_name = "FLDS"
  VAR1y_name = "TREFHT"

  VAR2x_name = "TMQ"
  VAR2y_name = "FLDS"

  suffix_DJF = "DJF"
  suffix_JJA = "JJA"

  outfile_scattered = "data_plotS03-scattered-"+VAR1x_name+"_"+VAR1y_name+"_"+VAR2x_name+"_"+VAR2y_name+"_zonal-SOM_SST-DJF_JJA-"+plotname_region+".nc"
  ;print(outfile_scattered)
  ;exit
  create_emptyNC(outfile_scattered)
  out_scattered = addfile(outfile_scattered,"w")

  casename_ensmean = (/ \
                "c10_ensmean-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_", \ 
                "c11_ensmean-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \ 
                "g01_ensmean-cesm111-F2000-mc6_rtr2_Scat.cam.h0.yy06_35-climo_", \
                "g02_ensmean-cesm111-F2000-mc6_rtr2_noScat.cam.h0.yy06_35-climo_"  /)

  datapath = "../data/"


  filename_ensmean_DJF = datapath+casename_ensmean+suffix_DJF+"-"+plotname_region+".nc"
  filename_ensmean_JJA = datapath+casename_ensmean+suffix_JJA+"-"+plotname_region+".nc"

;******************
; program start
;******************

;--------------------------------------
; read variables for scattered plot
;--------------------------------------
  opt_regline = True
    string_input = ""
    opt_regline@printout = False

  opt_areamean = True
    opt_areamean@return_type = "tji2tjOneD_avgi"

  print("read variables for scattered plots...")

  ;*** VAR1x_ensmean, DJF ***
  VAR1x_ensmean_ALL = convert_files_vars(filename_ensmean_DJF, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_ensmean_ALL)
  VAR1x_ensmean_SOM_Scat_DJF = VAR1x_ensmean_ALL(0,0,:)
  VAR1x_ensmean_SOM_noScat_DJF = VAR1x_ensmean_ALL(1,0,:)
  VAR1x_ensmean_SST_Scat_DJF = VAR1x_ensmean_ALL(2,0,:)
  VAR1x_ensmean_SST_noScat_DJF = VAR1x_ensmean_ALL(3,0,:)

  VAR1x_ensmean_SOM_ScatMnoScat_DJF = diff(VAR1x_ensmean_SOM_Scat_DJF,VAR1x_ensmean_SOM_noScat_DJF) 
  VAR1x_ensmean_SST_ScatMnoScat_DJF = diff(VAR1x_ensmean_SST_Scat_DJF,VAR1x_ensmean_SST_noScat_DJF) 

  ;printVarSummary(VAR1x_ensmean_SOM_ScatMnoScat_DJF)

  ;*** VAR1y_ensmean, DJF ***
  VAR1y_ensmean_ALL = convert_files_vars(filename_ensmean_DJF, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_ensmean_ALL)
  VAR1y_ensmean_SOM_Scat_DJF = VAR1y_ensmean_ALL(0,0,:)
  VAR1y_ensmean_SOM_noScat_DJF = VAR1y_ensmean_ALL(1,0,:)
  VAR1y_ensmean_SST_Scat_DJF = VAR1y_ensmean_ALL(2,0,:)
  VAR1y_ensmean_SST_noScat_DJF = VAR1y_ensmean_ALL(3,0,:)

  VAR1y_ensmean_SOM_ScatMnoScat_DJF = diff(VAR1y_ensmean_SOM_Scat_DJF,VAR1y_ensmean_SOM_noScat_DJF) 
  VAR1y_ensmean_SST_ScatMnoScat_DJF = diff(VAR1y_ensmean_SST_Scat_DJF,VAR1y_ensmean_SST_noScat_DJF) 

  ;*** VAR1x_ensmean, JJA ***
  VAR1x_ensmean_ALL = convert_files_vars(filename_ensmean_JJA, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_ensmean_ALL)
  VAR1x_ensmean_SOM_Scat_JJA = VAR1x_ensmean_ALL(0,0,:)
  VAR1x_ensmean_SOM_noScat_JJA = VAR1x_ensmean_ALL(1,0,:)
  VAR1x_ensmean_SST_Scat_JJA = VAR1x_ensmean_ALL(2,0,:)
  VAR1x_ensmean_SST_noScat_JJA = VAR1x_ensmean_ALL(3,0,:)

  VAR1x_ensmean_SOM_ScatMnoScat_JJA = diff(VAR1x_ensmean_SOM_Scat_JJA,VAR1x_ensmean_SOM_noScat_JJA) 
  VAR1x_ensmean_SST_ScatMnoScat_JJA = diff(VAR1x_ensmean_SST_Scat_JJA,VAR1x_ensmean_SST_noScat_JJA) 

  ;*** VAR1y_ensmean, JJA ***
  VAR1y_ensmean_ALL = convert_files_vars(filename_ensmean_JJA, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_ensmean_ALL)
  VAR1y_ensmean_SOM_Scat_JJA = VAR1y_ensmean_ALL(0,0,:)
  VAR1y_ensmean_SOM_noScat_JJA = VAR1y_ensmean_ALL(1,0,:)
  VAR1y_ensmean_SST_Scat_JJA = VAR1y_ensmean_ALL(2,0,:)
  VAR1y_ensmean_SST_noScat_JJA = VAR1y_ensmean_ALL(3,0,:)

  VAR1y_ensmean_SOM_ScatMnoScat_JJA = diff(VAR1y_ensmean_SOM_Scat_JJA,VAR1y_ensmean_SOM_noScat_JJA) 
  VAR1y_ensmean_SST_ScatMnoScat_JJA = diff(VAR1y_ensmean_SST_Scat_JJA,VAR1y_ensmean_SST_noScat_JJA) 

  ;*** VAR2x_ensmean, DJF ***
  VAR2x_ensmean_ALL = convert_files_vars(filename_ensmean_DJF, VAR2x_name, opt_areamean)
  ;printVarSummary(VAR2x_ensmean_ALL)
  VAR2x_ensmean_SOM_Scat_DJF = VAR2x_ensmean_ALL(0,0,:)
  VAR2x_ensmean_SOM_noScat_DJF = VAR2x_ensmean_ALL(1,0,:)
  VAR2x_ensmean_SST_Scat_DJF = VAR2x_ensmean_ALL(2,0,:)
  VAR2x_ensmean_SST_noScat_DJF = VAR2x_ensmean_ALL(3,0,:)

  VAR2x_ensmean_SOM_ScatMnoScat_DJF = diff(VAR2x_ensmean_SOM_Scat_DJF,VAR2x_ensmean_SOM_noScat_DJF) 
  VAR2x_ensmean_SST_ScatMnoScat_DJF = diff(VAR2x_ensmean_SST_Scat_DJF,VAR2x_ensmean_SST_noScat_DJF) 

  ;printVarSummary(VAR2x_ensmean_SOM_ScatMnoScat_DJF)

  ;*** VAR2y_ensmean, DJF ***
  VAR2y_ensmean_ALL = convert_files_vars(filename_ensmean_DJF, VAR2y_name, opt_areamean)
  ;printVarSummary(VAR2y_ensmean_ALL)
  VAR2y_ensmean_SOM_Scat_DJF = VAR2y_ensmean_ALL(0,0,:)
  VAR2y_ensmean_SOM_noScat_DJF = VAR2y_ensmean_ALL(1,0,:)
  VAR2y_ensmean_SST_Scat_DJF = VAR2y_ensmean_ALL(2,0,:)
  VAR2y_ensmean_SST_noScat_DJF = VAR2y_ensmean_ALL(3,0,:)

  VAR2y_ensmean_SOM_ScatMnoScat_DJF = diff(VAR2y_ensmean_SOM_Scat_DJF,VAR2y_ensmean_SOM_noScat_DJF) 
  VAR2y_ensmean_SST_ScatMnoScat_DJF = diff(VAR2y_ensmean_SST_Scat_DJF,VAR2y_ensmean_SST_noScat_DJF) 

  ;*** VAR2x_ensmean, JJA ***
  VAR2x_ensmean_ALL = convert_files_vars(filename_ensmean_JJA, VAR2x_name, opt_areamean)
  ;printVarSummary(VAR2x_ensmean_ALL)
  VAR2x_ensmean_SOM_Scat_JJA = VAR2x_ensmean_ALL(0,0,:)
  VAR2x_ensmean_SOM_noScat_JJA = VAR2x_ensmean_ALL(1,0,:)
  VAR2x_ensmean_SST_Scat_JJA = VAR2x_ensmean_ALL(2,0,:)
  VAR2x_ensmean_SST_noScat_JJA = VAR2x_ensmean_ALL(3,0,:)

  VAR2x_ensmean_SOM_ScatMnoScat_JJA = diff(VAR2x_ensmean_SOM_Scat_JJA,VAR2x_ensmean_SOM_noScat_JJA) 
  VAR2x_ensmean_SST_ScatMnoScat_JJA = diff(VAR2x_ensmean_SST_Scat_JJA,VAR2x_ensmean_SST_noScat_JJA) 

  ;*** VAR2y_ensmean, JJA ***
  VAR2y_ensmean_ALL = convert_files_vars(filename_ensmean_JJA, VAR2y_name, opt_areamean)
  ;printVarSummary(VAR2y_ensmean_ALL)
  VAR2y_ensmean_SOM_Scat_JJA = VAR2y_ensmean_ALL(0,0,:)
  VAR2y_ensmean_SOM_noScat_JJA = VAR2y_ensmean_ALL(1,0,:)
  VAR2y_ensmean_SST_Scat_JJA = VAR2y_ensmean_ALL(2,0,:)
  VAR2y_ensmean_SST_noScat_JJA = VAR2y_ensmean_ALL(3,0,:)

  VAR2y_ensmean_SOM_ScatMnoScat_JJA = diff(VAR2y_ensmean_SOM_Scat_JJA,VAR2y_ensmean_SOM_noScat_JJA) 
  VAR2y_ensmean_SST_ScatMnoScat_JJA = diff(VAR2y_ensmean_SST_Scat_JJA,VAR2y_ensmean_SST_noScat_JJA) 

;-----------
; write out 
;-----------

  print("write variables for scattered plots to "+outfile_scattered)

  ;*** ensmean, VAR1, DJF ***
  out_scattered->VAR1x_ensmean_SOM_ScatMnoScat_DJF = VAR1x_ensmean_SOM_ScatMnoScat_DJF
  out_scattered->VAR1y_ensmean_SOM_ScatMnoScat_DJF = VAR1y_ensmean_SOM_ScatMnoScat_DJF
  out_scattered->VAR1x_ensmean_SST_ScatMnoScat_DJF = VAR1x_ensmean_SST_ScatMnoScat_DJF
  out_scattered->VAR1y_ensmean_SST_ScatMnoScat_DJF = VAR1y_ensmean_SST_ScatMnoScat_DJF

  ;*** ensmean, VAR1, JJA ***
  out_scattered->VAR1x_ensmean_SOM_ScatMnoScat_JJA = VAR1x_ensmean_SOM_ScatMnoScat_JJA
  out_scattered->VAR1y_ensmean_SOM_ScatMnoScat_JJA = VAR1y_ensmean_SOM_ScatMnoScat_JJA
  out_scattered->VAR1x_ensmean_SST_ScatMnoScat_JJA = VAR1x_ensmean_SST_ScatMnoScat_JJA
  out_scattered->VAR1y_ensmean_SST_ScatMnoScat_JJA = VAR1y_ensmean_SST_ScatMnoScat_JJA

  ;*** ensmean, VAR2, DJF ***
  out_scattered->VAR2x_ensmean_SOM_ScatMnoScat_DJF = VAR2x_ensmean_SOM_ScatMnoScat_DJF
  out_scattered->VAR2y_ensmean_SOM_ScatMnoScat_DJF = VAR2y_ensmean_SOM_ScatMnoScat_DJF
  out_scattered->VAR2x_ensmean_SST_ScatMnoScat_DJF = VAR2x_ensmean_SST_ScatMnoScat_DJF
  out_scattered->VAR2y_ensmean_SST_ScatMnoScat_DJF = VAR2y_ensmean_SST_ScatMnoScat_DJF

  ;*** ensmean, VAR2, JJA ***
  out_scattered->VAR2x_ensmean_SOM_ScatMnoScat_JJA = VAR2x_ensmean_SOM_ScatMnoScat_JJA
  out_scattered->VAR2y_ensmean_SOM_ScatMnoScat_JJA = VAR2y_ensmean_SOM_ScatMnoScat_JJA
  out_scattered->VAR2x_ensmean_SST_ScatMnoScat_JJA = VAR2x_ensmean_SST_ScatMnoScat_JJA
  out_scattered->VAR2y_ensmean_SST_ScatMnoScat_JJA = VAR2y_ensmean_SST_ScatMnoScat_JJA

end
