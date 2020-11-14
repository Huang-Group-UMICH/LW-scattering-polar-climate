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

  VAR1x_name = "TMQ"
  VAR1y_name = "FLDS"

  suffix_DJF = "DJF"
  suffix_JJA = "JJA"

  outfile_scattered = "data_plotS05-scattered-"+VAR1x_name+"_"+VAR1y_name+"_areamean-SOM_SST-DJF_JJA-"+plotname_region+".nc"
  ;print(outfile_scattered)
  ;exit
  create_emptyNC(outfile_scattered)
  out_scattered = addfile(outfile_scattered,"w")

  casename_en00 = (/ \
                "c10-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_", \ 
                "c11-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \ 
                "g01-cesm111-F2000-mc6_rtr2_Scat.cam.h0.yy06_35-climo_", \
                "g02-cesm111-F2000-mc6_rtr2_noScat.cam.h0.yy06_35-climo_"  /)

  casename_en01 = (/ \
                "c10_en01-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_", \ 
                "c11_en01-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \ 
                "g01_en01-cesm111-F2000-mc6_rtr2_Scat.cam.h0.yy06_35-climo_", \
                "g02_en01-cesm111-F2000-mc6_rtr2_noScat.cam.h0.yy06_35-climo_"  /)

  casename_en02 = (/ \
                "c10_en02-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_", \ 
                "c11_en02-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \ 
                "g01_en02-cesm111-F2000-mc6_rtr2_Scat.cam.h0.yy06_35-climo_", \
                "g02_en02-cesm111-F2000-mc6_rtr2_noScat.cam.h0.yy06_35-climo_"  /)

  casename_en03 = (/ \
                "c10_en03-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_", \ 
                "c11_en03-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \ 
                "g01_en03-cesm111-F2000-mc6_rtr2_Scat.cam.h0.yy06_35-climo_", \
                "g02_en03-cesm111-F2000-mc6_rtr2_noScat.cam.h0.yy06_35-climo_"  /)


  datapath = "../data/"

  filename_en00_DJF = datapath+casename_en00+suffix_DJF+"-"+plotname_region+".nc"
  filename_en01_DJF = datapath+casename_en01+suffix_DJF+"-"+plotname_region+".nc"
  filename_en02_DJF = datapath+casename_en02+suffix_DJF+"-"+plotname_region+".nc"
  filename_en03_DJF = datapath+casename_en03+suffix_DJF+"-"+plotname_region+".nc"

  filename_en00_JJA = datapath+casename_en00+suffix_JJA+"-"+plotname_region+".nc"
  filename_en01_JJA = datapath+casename_en01+suffix_JJA+"-"+plotname_region+".nc"
  filename_en02_JJA = datapath+casename_en02+suffix_JJA+"-"+plotname_region+".nc"
  filename_en03_JJA = datapath+casename_en03+suffix_JJA+"-"+plotname_region+".nc"


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
    opt_areamean@return_type = "tji2t_avgji"

  print("read variables for scattered plots...")

  ;*** VAR1x_en00, DJF ***
  VAR1x_en00_ALL = convert_files_vars(filename_en00_DJF, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_en00_ALL)
  VAR1x_en00_SOM_Scat_DJF = VAR1x_en00_ALL(0,0,:)
  VAR1x_en00_SOM_noScat_DJF = VAR1x_en00_ALL(1,0,:)
  VAR1x_en00_SST_Scat_DJF = VAR1x_en00_ALL(2,0,:)
  VAR1x_en00_SST_noScat_DJF = VAR1x_en00_ALL(3,0,:)

  VAR1x_en00_SOM_ScatMnoScat_DJF = diff(VAR1x_en00_SOM_Scat_DJF,VAR1x_en00_SOM_noScat_DJF) 
  VAR1x_en00_SST_ScatMnoScat_DJF = diff(VAR1x_en00_SST_Scat_DJF,VAR1x_en00_SST_noScat_DJF) 

  ;*** VAR1y_en00, DJF ***
  VAR1y_en00_ALL = convert_files_vars(filename_en00_DJF, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_en00_ALL)
  VAR1y_en00_SOM_Scat_DJF = VAR1y_en00_ALL(0,0,:)
  VAR1y_en00_SOM_noScat_DJF = VAR1y_en00_ALL(1,0,:)
  VAR1y_en00_SST_Scat_DJF = VAR1y_en00_ALL(2,0,:)
  VAR1y_en00_SST_noScat_DJF = VAR1y_en00_ALL(3,0,:)

  VAR1y_en00_SOM_ScatMnoScat_DJF = diff(VAR1y_en00_SOM_Scat_DJF,VAR1y_en00_SOM_noScat_DJF) 
  VAR1y_en00_SST_ScatMnoScat_DJF = diff(VAR1y_en00_SST_Scat_DJF,VAR1y_en00_SST_noScat_DJF) 

  ;*** VAR1x_en00, JJA ***
  VAR1x_en00_ALL = convert_files_vars(filename_en00_JJA, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_en00_ALL)
  VAR1x_en00_SOM_Scat_JJA = VAR1x_en00_ALL(0,0,:)
  VAR1x_en00_SOM_noScat_JJA = VAR1x_en00_ALL(1,0,:)
  VAR1x_en00_SST_Scat_JJA = VAR1x_en00_ALL(2,0,:)
  VAR1x_en00_SST_noScat_JJA = VAR1x_en00_ALL(3,0,:)

  VAR1x_en00_SOM_ScatMnoScat_JJA = diff(VAR1x_en00_SOM_Scat_JJA,VAR1x_en00_SOM_noScat_JJA) 
  VAR1x_en00_SST_ScatMnoScat_JJA = diff(VAR1x_en00_SST_Scat_JJA,VAR1x_en00_SST_noScat_JJA) 

  ;*** VAR1y_en00, JJA ***
  VAR1y_en00_ALL = convert_files_vars(filename_en00_JJA, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_en00_ALL)
  VAR1y_en00_SOM_Scat_JJA = VAR1y_en00_ALL(0,0,:)
  VAR1y_en00_SOM_noScat_JJA = VAR1y_en00_ALL(1,0,:)
  VAR1y_en00_SST_Scat_JJA = VAR1y_en00_ALL(2,0,:)
  VAR1y_en00_SST_noScat_JJA = VAR1y_en00_ALL(3,0,:)

  VAR1y_en00_SOM_ScatMnoScat_JJA = diff(VAR1y_en00_SOM_Scat_JJA,VAR1y_en00_SOM_noScat_JJA) 
  VAR1y_en00_SST_ScatMnoScat_JJA = diff(VAR1y_en00_SST_Scat_JJA,VAR1y_en00_SST_noScat_JJA) 

  ;*** VAR1x_en01, DJF ***
  VAR1x_en01_ALL = convert_files_vars(filename_en01_DJF, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_en01_ALL)
  VAR1x_en01_SOM_Scat_DJF = VAR1x_en01_ALL(0,0,:)
  VAR1x_en01_SOM_noScat_DJF = VAR1x_en01_ALL(1,0,:)
  VAR1x_en01_SST_Scat_DJF = VAR1x_en01_ALL(2,0,:)
  VAR1x_en01_SST_noScat_DJF = VAR1x_en01_ALL(3,0,:)

  VAR1x_en01_SOM_ScatMnoScat_DJF = diff(VAR1x_en01_SOM_Scat_DJF,VAR1x_en01_SOM_noScat_DJF) 
  VAR1x_en01_SST_ScatMnoScat_DJF = diff(VAR1x_en01_SST_Scat_DJF,VAR1x_en01_SST_noScat_DJF) 

  ;*** VAR1y_en01, DJF ***
  VAR1y_en01_ALL = convert_files_vars(filename_en01_DJF, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_en01_ALL)
  VAR1y_en01_SOM_Scat_DJF = VAR1y_en01_ALL(0,0,:)
  VAR1y_en01_SOM_noScat_DJF = VAR1y_en01_ALL(1,0,:)
  VAR1y_en01_SST_Scat_DJF = VAR1y_en01_ALL(2,0,:)
  VAR1y_en01_SST_noScat_DJF = VAR1y_en01_ALL(3,0,:)

  VAR1y_en01_SOM_ScatMnoScat_DJF = diff(VAR1y_en01_SOM_Scat_DJF,VAR1y_en01_SOM_noScat_DJF) 
  VAR1y_en01_SST_ScatMnoScat_DJF = diff(VAR1y_en01_SST_Scat_DJF,VAR1y_en01_SST_noScat_DJF) 

  ;*** VAR1x_en01, JJA ***
  VAR1x_en01_ALL = convert_files_vars(filename_en01_JJA, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_en01_ALL)
  VAR1x_en01_SOM_Scat_JJA = VAR1x_en01_ALL(0,0,:)
  VAR1x_en01_SOM_noScat_JJA = VAR1x_en01_ALL(1,0,:)
  VAR1x_en01_SST_Scat_JJA = VAR1x_en01_ALL(2,0,:)
  VAR1x_en01_SST_noScat_JJA = VAR1x_en01_ALL(3,0,:)

  VAR1x_en01_SOM_ScatMnoScat_JJA = diff(VAR1x_en01_SOM_Scat_JJA,VAR1x_en01_SOM_noScat_JJA) 
  VAR1x_en01_SST_ScatMnoScat_JJA = diff(VAR1x_en01_SST_Scat_JJA,VAR1x_en01_SST_noScat_JJA) 

  ;*** VAR1y_en01, JJA ***
  VAR1y_en01_ALL = convert_files_vars(filename_en01_JJA, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_en01_ALL)
  VAR1y_en01_SOM_Scat_JJA = VAR1y_en01_ALL(0,0,:)
  VAR1y_en01_SOM_noScat_JJA = VAR1y_en01_ALL(1,0,:)
  VAR1y_en01_SST_Scat_JJA = VAR1y_en01_ALL(2,0,:)
  VAR1y_en01_SST_noScat_JJA = VAR1y_en01_ALL(3,0,:)

  VAR1y_en01_SOM_ScatMnoScat_JJA = diff(VAR1y_en01_SOM_Scat_JJA,VAR1y_en01_SOM_noScat_JJA) 
  VAR1y_en01_SST_ScatMnoScat_JJA = diff(VAR1y_en01_SST_Scat_JJA,VAR1y_en01_SST_noScat_JJA) 

  ;*** VAR1x_en02, DJF ***
  VAR1x_en02_ALL = convert_files_vars(filename_en02_DJF, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_en02_ALL)
  VAR1x_en02_SOM_Scat_DJF = VAR1x_en02_ALL(0,0,:)
  VAR1x_en02_SOM_noScat_DJF = VAR1x_en02_ALL(1,0,:)
  VAR1x_en02_SST_Scat_DJF = VAR1x_en02_ALL(2,0,:)
  VAR1x_en02_SST_noScat_DJF = VAR1x_en02_ALL(3,0,:)

  VAR1x_en02_SOM_ScatMnoScat_DJF = diff(VAR1x_en02_SOM_Scat_DJF,VAR1x_en02_SOM_noScat_DJF) 
  VAR1x_en02_SST_ScatMnoScat_DJF = diff(VAR1x_en02_SST_Scat_DJF,VAR1x_en02_SST_noScat_DJF) 

  ;*** VAR1y_en02, DJF ***
  VAR1y_en02_ALL = convert_files_vars(filename_en02_DJF, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_en02_ALL)
  VAR1y_en02_SOM_Scat_DJF = VAR1y_en02_ALL(0,0,:)
  VAR1y_en02_SOM_noScat_DJF = VAR1y_en02_ALL(1,0,:)
  VAR1y_en02_SST_Scat_DJF = VAR1y_en02_ALL(2,0,:)
  VAR1y_en02_SST_noScat_DJF = VAR1y_en02_ALL(3,0,:)

  VAR1y_en02_SOM_ScatMnoScat_DJF = diff(VAR1y_en02_SOM_Scat_DJF,VAR1y_en02_SOM_noScat_DJF) 
  VAR1y_en02_SST_ScatMnoScat_DJF = diff(VAR1y_en02_SST_Scat_DJF,VAR1y_en02_SST_noScat_DJF) 

  ;*** VAR1x_en02, JJA ***
  VAR1x_en02_ALL = convert_files_vars(filename_en02_JJA, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_en02_ALL)
  VAR1x_en02_SOM_Scat_JJA = VAR1x_en02_ALL(0,0,:)
  VAR1x_en02_SOM_noScat_JJA = VAR1x_en02_ALL(1,0,:)
  VAR1x_en02_SST_Scat_JJA = VAR1x_en02_ALL(2,0,:)
  VAR1x_en02_SST_noScat_JJA = VAR1x_en02_ALL(3,0,:)

  VAR1x_en02_SOM_ScatMnoScat_JJA = diff(VAR1x_en02_SOM_Scat_JJA,VAR1x_en02_SOM_noScat_JJA) 
  VAR1x_en02_SST_ScatMnoScat_JJA = diff(VAR1x_en02_SST_Scat_JJA,VAR1x_en02_SST_noScat_JJA) 

  ;*** VAR1y_en02, JJA ***
  VAR1y_en02_ALL = convert_files_vars(filename_en02_JJA, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_en02_ALL)
  VAR1y_en02_SOM_Scat_JJA = VAR1y_en02_ALL(0,0,:)
  VAR1y_en02_SOM_noScat_JJA = VAR1y_en02_ALL(1,0,:)
  VAR1y_en02_SST_Scat_JJA = VAR1y_en02_ALL(2,0,:)
  VAR1y_en02_SST_noScat_JJA = VAR1y_en02_ALL(3,0,:)

  VAR1y_en02_SOM_ScatMnoScat_JJA = diff(VAR1y_en02_SOM_Scat_JJA,VAR1y_en02_SOM_noScat_JJA) 
  VAR1y_en02_SST_ScatMnoScat_JJA = diff(VAR1y_en02_SST_Scat_JJA,VAR1y_en02_SST_noScat_JJA) 

  ;*** VAR1x_en03, DJF ***
  VAR1x_en03_ALL = convert_files_vars(filename_en03_DJF, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_en03_ALL)
  VAR1x_en03_SOM_Scat_DJF = VAR1x_en03_ALL(0,0,:)
  VAR1x_en03_SOM_noScat_DJF = VAR1x_en03_ALL(1,0,:)
  VAR1x_en03_SST_Scat_DJF = VAR1x_en03_ALL(2,0,:)
  VAR1x_en03_SST_noScat_DJF = VAR1x_en03_ALL(3,0,:)

  VAR1x_en03_SOM_ScatMnoScat_DJF = diff(VAR1x_en03_SOM_Scat_DJF,VAR1x_en03_SOM_noScat_DJF) 
  VAR1x_en03_SST_ScatMnoScat_DJF = diff(VAR1x_en03_SST_Scat_DJF,VAR1x_en03_SST_noScat_DJF) 

  ;*** VAR1y_en03, DJF ***
  VAR1y_en03_ALL = convert_files_vars(filename_en03_DJF, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_en03_ALL)
  VAR1y_en03_SOM_Scat_DJF = VAR1y_en03_ALL(0,0,:)
  VAR1y_en03_SOM_noScat_DJF = VAR1y_en03_ALL(1,0,:)
  VAR1y_en03_SST_Scat_DJF = VAR1y_en03_ALL(2,0,:)
  VAR1y_en03_SST_noScat_DJF = VAR1y_en03_ALL(3,0,:)

  VAR1y_en03_SOM_ScatMnoScat_DJF = diff(VAR1y_en03_SOM_Scat_DJF,VAR1y_en03_SOM_noScat_DJF) 
  VAR1y_en03_SST_ScatMnoScat_DJF = diff(VAR1y_en03_SST_Scat_DJF,VAR1y_en03_SST_noScat_DJF) 

  ;*** VAR1x_en03, JJA ***
  VAR1x_en03_ALL = convert_files_vars(filename_en03_JJA, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_en03_ALL)
  VAR1x_en03_SOM_Scat_JJA = VAR1x_en03_ALL(0,0,:)
  VAR1x_en03_SOM_noScat_JJA = VAR1x_en03_ALL(1,0,:)
  VAR1x_en03_SST_Scat_JJA = VAR1x_en03_ALL(2,0,:)
  VAR1x_en03_SST_noScat_JJA = VAR1x_en03_ALL(3,0,:)

  VAR1x_en03_SOM_ScatMnoScat_JJA = diff(VAR1x_en03_SOM_Scat_JJA,VAR1x_en03_SOM_noScat_JJA) 
  VAR1x_en03_SST_ScatMnoScat_JJA = diff(VAR1x_en03_SST_Scat_JJA,VAR1x_en03_SST_noScat_JJA) 

  ;*** VAR1y_en03, JJA ***
  VAR1y_en03_ALL = convert_files_vars(filename_en03_JJA, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_en03_ALL)
  VAR1y_en03_SOM_Scat_JJA = VAR1y_en03_ALL(0,0,:)
  VAR1y_en03_SOM_noScat_JJA = VAR1y_en03_ALL(1,0,:)
  VAR1y_en03_SST_Scat_JJA = VAR1y_en03_ALL(2,0,:)
  VAR1y_en03_SST_noScat_JJA = VAR1y_en03_ALL(3,0,:)

  VAR1y_en03_SOM_ScatMnoScat_JJA = diff(VAR1y_en03_SOM_Scat_JJA,VAR1y_en03_SOM_noScat_JJA) 
  VAR1y_en03_SST_ScatMnoScat_JJA = diff(VAR1y_en03_SST_Scat_JJA,VAR1y_en03_SST_noScat_JJA) 

;-----------
; write out 
;-----------

  print("write variables for scattered plots to "+outfile_scattered)

  ;*** en00, VAR1, DJF ***
  out_scattered->VAR1x_en00_SOM_ScatMnoScat_DJF = VAR1x_en00_SOM_ScatMnoScat_DJF
  out_scattered->VAR1y_en00_SOM_ScatMnoScat_DJF = VAR1y_en00_SOM_ScatMnoScat_DJF
  out_scattered->VAR1x_en00_SST_ScatMnoScat_DJF = VAR1x_en00_SST_ScatMnoScat_DJF
  out_scattered->VAR1y_en00_SST_ScatMnoScat_DJF = VAR1y_en00_SST_ScatMnoScat_DJF

  ;*** en00, VAR1, JJA ***
  out_scattered->VAR1x_en00_SOM_ScatMnoScat_JJA = VAR1x_en00_SOM_ScatMnoScat_JJA
  out_scattered->VAR1y_en00_SOM_ScatMnoScat_JJA = VAR1y_en00_SOM_ScatMnoScat_JJA
  out_scattered->VAR1x_en00_SST_ScatMnoScat_JJA = VAR1x_en00_SST_ScatMnoScat_JJA
  out_scattered->VAR1y_en00_SST_ScatMnoScat_JJA = VAR1y_en00_SST_ScatMnoScat_JJA

  ;*** en01, VAR1, DJF ***
  out_scattered->VAR1x_en01_SOM_ScatMnoScat_DJF = VAR1x_en01_SOM_ScatMnoScat_DJF
  out_scattered->VAR1y_en01_SOM_ScatMnoScat_DJF = VAR1y_en01_SOM_ScatMnoScat_DJF
  out_scattered->VAR1x_en01_SST_ScatMnoScat_DJF = VAR1x_en01_SST_ScatMnoScat_DJF
  out_scattered->VAR1y_en01_SST_ScatMnoScat_DJF = VAR1y_en01_SST_ScatMnoScat_DJF

  ;*** en01, VAR1, JJA ***
  out_scattered->VAR1x_en01_SOM_ScatMnoScat_JJA = VAR1x_en01_SOM_ScatMnoScat_JJA
  out_scattered->VAR1y_en01_SOM_ScatMnoScat_JJA = VAR1y_en01_SOM_ScatMnoScat_JJA
  out_scattered->VAR1x_en01_SST_ScatMnoScat_JJA = VAR1x_en01_SST_ScatMnoScat_JJA
  out_scattered->VAR1y_en01_SST_ScatMnoScat_JJA = VAR1y_en01_SST_ScatMnoScat_JJA

  ;*** en02, VAR1, DJF ***
  out_scattered->VAR1x_en02_SOM_ScatMnoScat_DJF = VAR1x_en02_SOM_ScatMnoScat_DJF
  out_scattered->VAR1y_en02_SOM_ScatMnoScat_DJF = VAR1y_en02_SOM_ScatMnoScat_DJF
  out_scattered->VAR1x_en02_SST_ScatMnoScat_DJF = VAR1x_en02_SST_ScatMnoScat_DJF
  out_scattered->VAR1y_en02_SST_ScatMnoScat_DJF = VAR1y_en02_SST_ScatMnoScat_DJF

  ;*** en02, VAR1, JJA ***
  out_scattered->VAR1x_en02_SOM_ScatMnoScat_JJA = VAR1x_en02_SOM_ScatMnoScat_JJA
  out_scattered->VAR1y_en02_SOM_ScatMnoScat_JJA = VAR1y_en02_SOM_ScatMnoScat_JJA
  out_scattered->VAR1x_en02_SST_ScatMnoScat_JJA = VAR1x_en02_SST_ScatMnoScat_JJA
  out_scattered->VAR1y_en02_SST_ScatMnoScat_JJA = VAR1y_en02_SST_ScatMnoScat_JJA

  ;*** en03, VAR1, DJF ***
  out_scattered->VAR1x_en03_SOM_ScatMnoScat_DJF = VAR1x_en03_SOM_ScatMnoScat_DJF
  out_scattered->VAR1y_en03_SOM_ScatMnoScat_DJF = VAR1y_en03_SOM_ScatMnoScat_DJF
  out_scattered->VAR1x_en03_SST_ScatMnoScat_DJF = VAR1x_en03_SST_ScatMnoScat_DJF
  out_scattered->VAR1y_en03_SST_ScatMnoScat_DJF = VAR1y_en03_SST_ScatMnoScat_DJF

  ;*** en03, VAR1, JJA ***
  out_scattered->VAR1x_en03_SOM_ScatMnoScat_JJA = VAR1x_en03_SOM_ScatMnoScat_JJA
  out_scattered->VAR1y_en03_SOM_ScatMnoScat_JJA = VAR1y_en03_SOM_ScatMnoScat_JJA
  out_scattered->VAR1x_en03_SST_ScatMnoScat_JJA = VAR1x_en03_SST_ScatMnoScat_JJA
  out_scattered->VAR1y_en03_SST_ScatMnoScat_JJA = VAR1y_en03_SST_ScatMnoScat_JJA

end
