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
  ;plotname_region = "Antarctic"

  VAR1x_name = "FLDS"
  VAR1y_name = "TREFHT"

  VAR2x_name = "TMQ"
  VAR2y_name = "FLDS"

  VAR1_profile = "T"
  VAR2_profile = "Q"
  VAR3_profile = "CLOUD"
  VAR_ALL_profile = (/VAR1_profile,VAR2_profile,VAR3_profile/)

  outfile_scattered = "data_plot03-scattered-"+VAR1x_name+"_"+VAR1y_name+"_"+VAR2x_name+"_"+VAR2y_name+"_areamean-SOM_SST-DJF_JJA-"+plotname_region+".nc"
  create_emptyNC(outfile_scattered)
  out_scattered = addfile(outfile_scattered,"w")

  outfile_profile = "data_plot03-profile-"+VAR1_profile+"_"+VAR2_profile+"_"+VAR3_profile+"_areamean-SOM_SST-DJF_JJA-"+plotname_region+".nc"
  create_emptyNC(outfile_profile)
  out_profile = addfile(outfile_profile,"w")

  suffix_DJF = "DJF"
  suffix_JJA = "JJA"

  casename_ensmean = (/ \
                "c10_ensmean-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_", \ 
                "c11_ensmean-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \ 
                "g01_ensmean-cesm111-F2000-mc6_rtr2_Scat.cam.h0.yy06_35-climo_", \
                "g02_ensmean-cesm111-F2000-mc6_rtr2_noScat.cam.h0.yy06_35-climo_"  /)

  casename_SOM_Scat = (/ \
                "c10-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_", \
                "c10_en01-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_", \
                "c10_en02-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_", \
                "c10_en03-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_"/)

  casename_SOM_noScat = (/ \
                "c11-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \
                "c11_en01-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \
                "c11_en02-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \
                "c11_en03-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_"/)

  casename_SST_Scat = (/ \
                "g01-cesm111-F2000-mc6_rtr2_Scat.cam.h0.yy06_35-climo_",  \
                "g01_en01-cesm111-F2000-mc6_rtr2_Scat.cam.h0.yy06_35-climo_", \
                "g01_en02-cesm111-F2000-mc6_rtr2_Scat.cam.h0.yy06_35-climo_", \
                "g01_en03-cesm111-F2000-mc6_rtr2_Scat.cam.h0.yy06_35-climo_"/)

  casename_SST_noScat = (/ \
                "g02-cesm111-F2000-mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \
                "g02_en01-cesm111-F2000-mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \
                "g02_en02-cesm111-F2000-mc6_rtr2_noScat.cam.h0.yy06_35-climo_", \
                "g02_en03-cesm111-F2000-mc6_rtr2_noScat.cam.h0.yy06_35-climo_"/)

  datapath = "../data/"

  filename_ensmean_DJF = datapath+casename_ensmean+suffix_DJF+"-"+plotname_region+".nc"
  filename_Scat_SOM_DJF = datapath+casename_SOM_Scat+suffix_DJF+"-"+plotname_region+".nc"
  filename_noScat_SOM_DJF = datapath+casename_SOM_noScat+suffix_DJF+"-"+plotname_region+".nc"
  filename_Scat_SST_DJF = datapath+casename_SST_Scat+suffix_DJF+"-"+plotname_region+".nc"
  filename_noScat_SST_DJF = datapath+casename_SST_noScat+suffix_DJF+"-"+plotname_region+".nc"

  filename_ensmean_JJA = datapath+casename_ensmean+suffix_JJA+"-"+plotname_region+".nc"
  filename_Scat_SOM_JJA = datapath+casename_SOM_Scat+suffix_JJA+"-"+plotname_region+".nc"
  filename_noScat_SOM_JJA = datapath+casename_SOM_noScat+suffix_JJA+"-"+plotname_region+".nc"
  filename_Scat_SST_JJA = datapath+casename_SST_Scat+suffix_JJA+"-"+plotname_region+".nc"
  filename_noScat_SST_JJA = datapath+casename_SST_noScat+suffix_JJA+"-"+plotname_region+".nc"

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

  ;*** VAR1x, DJF ***
  VAR1x_ALL = convert_files_vars(filename_ensmean_DJF, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_ALL)
  VAR1x_SOM_Scat_DJF = VAR1x_ALL(0,0,:)
  VAR1x_SOM_noScat_DJF = VAR1x_ALL(1,0,:)
  VAR1x_SST_Scat_DJF = VAR1x_ALL(2,0,:)
  VAR1x_SST_noScat_DJF = VAR1x_ALL(3,0,:)

  VAR1x_SOM_ScatMnoScat_DJF = diff(VAR1x_SOM_Scat_DJF,VAR1x_SOM_noScat_DJF) 
  VAR1x_SST_ScatMnoScat_DJF = diff(VAR1x_SST_Scat_DJF,VAR1x_SST_noScat_DJF) 

  ;*** VAR1y, DJF ***
  VAR1y_ALL = convert_files_vars(filename_ensmean_DJF, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_ALL)
  VAR1y_SOM_Scat_DJF = VAR1y_ALL(0,0,:)
  VAR1y_SOM_noScat_DJF = VAR1y_ALL(1,0,:)
  VAR1y_SST_Scat_DJF = VAR1y_ALL(2,0,:)
  VAR1y_SST_noScat_DJF = VAR1y_ALL(3,0,:)

  VAR1y_SOM_ScatMnoScat_DJF = diff(VAR1y_SOM_Scat_DJF,VAR1y_SOM_noScat_DJF) 
  VAR1y_SST_ScatMnoScat_DJF = diff(VAR1y_SST_Scat_DJF,VAR1y_SST_noScat_DJF) 

  ;*** VAR2x, DJF ***
  VAR2x_ALL = convert_files_vars(filename_ensmean_DJF, VAR2x_name, opt_areamean)
  ;printVarSummary(VAR2x_ALL)
  VAR2x_SOM_Scat_DJF = VAR2x_ALL(0,0,:)
  VAR2x_SOM_noScat_DJF = VAR2x_ALL(1,0,:)
  VAR2x_SST_Scat_DJF = VAR2x_ALL(2,0,:)
  VAR2x_SST_noScat_DJF = VAR2x_ALL(3,0,:)

  VAR2x_SOM_ScatMnoScat_DJF = diff(VAR2x_SOM_Scat_DJF,VAR2x_SOM_noScat_DJF) 
  VAR2x_SST_ScatMnoScat_DJF = diff(VAR2x_SST_Scat_DJF,VAR2x_SST_noScat_DJF) 

  ;*** VAR2y, DJF ***
  VAR2y_ALL = convert_files_vars(filename_ensmean_DJF, VAR2y_name, opt_areamean)
  ;printVarSummary(VAR2y_ALL)
  VAR2y_SOM_Scat_DJF = VAR2y_ALL(0,0,:)
  VAR2y_SOM_noScat_DJF = VAR2y_ALL(1,0,:)
  VAR2y_SST_Scat_DJF = VAR2y_ALL(2,0,:)
  VAR2y_SST_noScat_DJF = VAR2y_ALL(3,0,:)

  VAR2y_SOM_ScatMnoScat_DJF = diff(VAR2y_SOM_Scat_DJF,VAR2y_SOM_noScat_DJF) 
  VAR2y_SST_ScatMnoScat_DJF = diff(VAR2y_SST_Scat_DJF,VAR2y_SST_noScat_DJF) 

  ;*** VAR1x, JJA ***
  VAR1x_ALL = convert_files_vars(filename_ensmean_JJA, VAR1x_name, opt_areamean)
  ;printVarSummary(VAR1x_ALL)
  VAR1x_SOM_Scat_JJA = VAR1x_ALL(0,0,:)
  VAR1x_SOM_noScat_JJA = VAR1x_ALL(1,0,:)
  VAR1x_SST_Scat_JJA = VAR1x_ALL(2,0,:)
  VAR1x_SST_noScat_JJA = VAR1x_ALL(3,0,:)

  VAR1x_SOM_ScatMnoScat_JJA = diff(VAR1x_SOM_Scat_JJA,VAR1x_SOM_noScat_JJA) 
  VAR1x_SST_ScatMnoScat_JJA = diff(VAR1x_SST_Scat_JJA,VAR1x_SST_noScat_JJA) 

  ;*** VAR1y, JJA ***
  VAR1y_ALL = convert_files_vars(filename_ensmean_JJA, VAR1y_name, opt_areamean)
  ;printVarSummary(VAR1y_ALL)
  VAR1y_SOM_Scat_JJA = VAR1y_ALL(0,0,:)
  VAR1y_SOM_noScat_JJA = VAR1y_ALL(1,0,:)
  VAR1y_SST_Scat_JJA = VAR1y_ALL(2,0,:)
  VAR1y_SST_noScat_JJA = VAR1y_ALL(3,0,:)

  VAR1y_SOM_ScatMnoScat_JJA = diff(VAR1y_SOM_Scat_JJA,VAR1y_SOM_noScat_JJA) 
  VAR1y_SST_ScatMnoScat_JJA = diff(VAR1y_SST_Scat_JJA,VAR1y_SST_noScat_JJA) 

  ;*** VAR2x, JJA ***
  VAR2x_ALL = convert_files_vars(filename_ensmean_JJA, VAR2x_name, opt_areamean)
  ;printVarSummary(VAR2x_ALL)
  VAR2x_SOM_Scat_JJA = VAR2x_ALL(0,0,:)
  VAR2x_SOM_noScat_JJA = VAR2x_ALL(1,0,:)
  VAR2x_SST_Scat_JJA = VAR2x_ALL(2,0,:)
  VAR2x_SST_noScat_JJA = VAR2x_ALL(3,0,:)

  VAR2x_SOM_ScatMnoScat_JJA = diff(VAR2x_SOM_Scat_JJA,VAR2x_SOM_noScat_JJA) 
  VAR2x_SST_ScatMnoScat_JJA = diff(VAR2x_SST_Scat_JJA,VAR2x_SST_noScat_JJA) 

  ;*** VAR2y, JJA ***
  VAR2y_ALL = convert_files_vars(filename_ensmean_JJA, VAR2y_name, opt_areamean)
  ;printVarSummary(VAR2y_ALL)
  VAR2y_SOM_Scat_JJA = VAR2y_ALL(0,0,:)
  VAR2y_SOM_noScat_JJA = VAR2y_ALL(1,0,:)
  VAR2y_SST_Scat_JJA = VAR2y_ALL(2,0,:)
  VAR2y_SST_noScat_JJA = VAR2y_ALL(3,0,:)

  VAR2y_SOM_ScatMnoScat_JJA = diff(VAR2y_SOM_Scat_JJA,VAR2y_SOM_noScat_JJA) 
  VAR2y_SST_ScatMnoScat_JJA = diff(VAR2y_SST_Scat_JJA,VAR2y_SST_noScat_JJA) 

;-----------
; write out 
;-----------

  print("write variables for scattered plots to "+outfile_scattered)

  ;*** VAR1, DJF ***
  out_scattered->VAR1x_SOM_ScatMnoScat_DJF = VAR1x_SOM_ScatMnoScat_DJF
  out_scattered->VAR1y_SOM_ScatMnoScat_DJF = VAR1y_SOM_ScatMnoScat_DJF
  out_scattered->VAR1x_SST_ScatMnoScat_DJF = VAR1x_SST_ScatMnoScat_DJF
  out_scattered->VAR1y_SST_ScatMnoScat_DJF = VAR1y_SST_ScatMnoScat_DJF

  ;*** VAR2, DJF ***
  out_scattered->VAR2x_SOM_ScatMnoScat_DJF = VAR2x_SOM_ScatMnoScat_DJF
  out_scattered->VAR2y_SOM_ScatMnoScat_DJF = VAR2y_SOM_ScatMnoScat_DJF
  out_scattered->VAR2x_SST_ScatMnoScat_DJF = VAR2x_SST_ScatMnoScat_DJF
  out_scattered->VAR2y_SST_ScatMnoScat_DJF = VAR2y_SST_ScatMnoScat_DJF

  ;*** VAR1, JJA ***
  out_scattered->VAR1x_SOM_ScatMnoScat_JJA = VAR1x_SOM_ScatMnoScat_JJA
  out_scattered->VAR1y_SOM_ScatMnoScat_JJA = VAR1y_SOM_ScatMnoScat_JJA
  out_scattered->VAR1x_SST_ScatMnoScat_JJA = VAR1x_SST_ScatMnoScat_JJA
  out_scattered->VAR1y_SST_ScatMnoScat_JJA = VAR1y_SST_ScatMnoScat_JJA

  ;*** VAR2, JJA ***
  out_scattered->VAR2x_SOM_ScatMnoScat_JJA = VAR2x_SOM_ScatMnoScat_JJA
  out_scattered->VAR2y_SOM_ScatMnoScat_JJA = VAR2y_SOM_ScatMnoScat_JJA
  out_scattered->VAR2x_SST_ScatMnoScat_JJA = VAR2x_SST_ScatMnoScat_JJA
  out_scattered->VAR2y_SST_ScatMnoScat_JJA = VAR2y_SST_ScatMnoScat_JJA
;
;------------------------------
; read variables for profile
;------------------------------
  opt_profile = True
    opt_profile@return_type = "ftlji2fp_p_avgtji"
    opt_profile@write2file = True

  print("read variables for profile plots...")

  ;*** Scat, SOM, DJF ***
  VARS_profile_Scat_SOM_DJF = convert_files_vars(filename_Scat_SOM_DJF, VAR_ALL_profile, opt_profile)
  ;printVarSummary(VARS_profile_Scat_SOM_DJF)
  VAR1p_Scat_SOM_DJF = VARS_profile_Scat_SOM_DJF(:,0,:)
  VAR2p_Scat_SOM_DJF = VARS_profile_Scat_SOM_DJF(:,1,:)
  VAR3p_Scat_SOM_DJF = VARS_profile_Scat_SOM_DJF(:,2,:)

  ;*** noScat, SOM, DJF ***
  VARS_profile_noScat_SOM_DJF = convert_files_vars(filename_noScat_SOM_DJF, VAR_ALL_profile, opt_profile)
  ;printVarSummary(VARS_profile_noScat_SOM_DJF)
  VAR1p_noScat_SOM_DJF = VARS_profile_noScat_SOM_DJF(:,0,:)
  VAR2p_noScat_SOM_DJF = VARS_profile_noScat_SOM_DJF(:,1,:)
  VAR3p_noScat_SOM_DJF = VARS_profile_noScat_SOM_DJF(:,2,:)

  VAR1p_ScatMnoScat_SOM_DJF = diff(VAR1p_Scat_SOM_DJF, VAR1p_noScat_SOM_DJF)
  VAR2p_ScatMnoScat_SOM_DJF = diff(VAR2p_Scat_SOM_DJF, VAR2p_noScat_SOM_DJF)
  VAR3p_ScatMnoScat_SOM_DJF = diff(VAR3p_Scat_SOM_DJF, VAR3p_noScat_SOM_DJF)
  ;printVarSummary(VAR1p_ScatMnoScat_SOM_DJF)

  ;*** Scat, SST, DJF ***
  VARS_profile_Scat_SST_DJF = convert_files_vars(filename_Scat_SST_DJF, VAR_ALL_profile, opt_profile)
  ;printVarSummary(VARS_profile_Scat_SST_DJF)
  VAR1p_Scat_SST_DJF = VARS_profile_Scat_SST_DJF(:,0,:)
  VAR2p_Scat_SST_DJF = VARS_profile_Scat_SST_DJF(:,1,:)
  VAR3p_Scat_SST_DJF = VARS_profile_Scat_SST_DJF(:,2,:)

  ;*** noScat, SST, DJF ***
  VARS_profile_noScat_SST_DJF = convert_files_vars(filename_noScat_SST_DJF, VAR_ALL_profile, opt_profile)
  ;printVarSummary(VARS_profile_noScat_SST_DJF)
  VAR1p_noScat_SST_DJF = VARS_profile_noScat_SST_DJF(:,0,:)
  VAR2p_noScat_SST_DJF = VARS_profile_noScat_SST_DJF(:,1,:)
  VAR3p_noScat_SST_DJF = VARS_profile_noScat_SST_DJF(:,2,:)

  VAR1p_ScatMnoScat_SST_DJF = diff(VAR1p_Scat_SST_DJF, VAR1p_noScat_SST_DJF)
  VAR2p_ScatMnoScat_SST_DJF = diff(VAR2p_Scat_SST_DJF, VAR2p_noScat_SST_DJF)
  VAR3p_ScatMnoScat_SST_DJF = diff(VAR3p_Scat_SST_DJF, VAR3p_noScat_SST_DJF)

  ;*** Scat, SOM, JJA ***
  VARS_profile_Scat_SOM_JJA = convert_files_vars(filename_Scat_SOM_JJA, VAR_ALL_profile, opt_profile)
  ;printVarSummary(VARS_profile_Scat_SOM_JJA)
  VAR1p_Scat_SOM_JJA = VARS_profile_Scat_SOM_JJA(:,0,:)
  VAR2p_Scat_SOM_JJA = VARS_profile_Scat_SOM_JJA(:,1,:)
  VAR3p_Scat_SOM_JJA = VARS_profile_Scat_SOM_JJA(:,2,:)

  ;*** noScat, SOM, JJA ***
  VARS_profile_noScat_SOM_JJA = convert_files_vars(filename_noScat_SOM_JJA, VAR_ALL_profile, opt_profile)
  ;printVarSummary(VARS_profile_noScat_SOM_JJA)
  VAR1p_noScat_SOM_JJA = VARS_profile_noScat_SOM_JJA(:,0,:)
  VAR2p_noScat_SOM_JJA = VARS_profile_noScat_SOM_JJA(:,1,:)
  VAR3p_noScat_SOM_JJA = VARS_profile_noScat_SOM_JJA(:,2,:)

  VAR1p_ScatMnoScat_SOM_JJA = diff(VAR1p_Scat_SOM_JJA, VAR1p_noScat_SOM_JJA)
  VAR2p_ScatMnoScat_SOM_JJA = diff(VAR2p_Scat_SOM_JJA, VAR2p_noScat_SOM_JJA)
  VAR3p_ScatMnoScat_SOM_JJA = diff(VAR3p_Scat_SOM_JJA, VAR3p_noScat_SOM_JJA)

  ;*** Scat, SST, JJA ***
  VARS_profile_Scat_SST_JJA = convert_files_vars(filename_Scat_SST_JJA, VAR_ALL_profile, opt_profile)
  ;printVarSummary(VARS_profile_Scat_SST_JJA)
  VAR1p_Scat_SST_JJA = VARS_profile_Scat_SST_JJA(:,0,:)
  VAR2p_Scat_SST_JJA = VARS_profile_Scat_SST_JJA(:,1,:)
  VAR3p_Scat_SST_JJA = VARS_profile_Scat_SST_JJA(:,2,:)

  ;*** noScat, SST, JJA ***
  VARS_profile_noScat_SST_JJA = convert_files_vars(filename_noScat_SST_JJA, VAR_ALL_profile, opt_profile)
  ;printVarSummary(VARS_profile_noScat_SST_JJA)
  VAR1p_noScat_SST_JJA = VARS_profile_noScat_SST_JJA(:,0,:)
  VAR2p_noScat_SST_JJA = VARS_profile_noScat_SST_JJA(:,1,:)
  VAR3p_noScat_SST_JJA = VARS_profile_noScat_SST_JJA(:,2,:)

  VAR1p_ScatMnoScat_SST_JJA = diff(VAR1p_Scat_SST_JJA, VAR1p_noScat_SST_JJA)
  VAR2p_ScatMnoScat_SST_JJA = diff(VAR2p_Scat_SST_JJA, VAR2p_noScat_SST_JJA)
  VAR3p_ScatMnoScat_SST_JJA = diff(VAR3p_Scat_SST_JJA, VAR3p_noScat_SST_JJA)

  ;printVarSummary(VAR1p_ScatMnoScat_SST_JJA)

  nn = dimsizes(VAR1p_ScatMnoScat_SST_JJA)
  n1 = nn(0)
  n2 = nn(1)
  nfile = 2

  VAR1p_ScatMnoScat_SOMSST_DJF = reshape( (/VAR1p_ScatMnoScat_SOM_DJF, VAR1p_ScatMnoScat_SST_DJF/), (/nfile*n1,n2/) )
  VAR2p_ScatMnoScat_SOMSST_DJF = reshape( (/VAR2p_ScatMnoScat_SOM_DJF, VAR2p_ScatMnoScat_SST_DJF/), (/nfile*n1,n2/) )
  VAR3p_ScatMnoScat_SOMSST_DJF = reshape( (/VAR3p_ScatMnoScat_SOM_DJF, VAR3p_ScatMnoScat_SST_DJF/), (/nfile*n1,n2/) )

  VAR1p_ScatMnoScat_SOMSST_JJA = reshape( (/VAR1p_ScatMnoScat_SOM_JJA, VAR1p_ScatMnoScat_SST_JJA/), (/nfile*n1,n2/) )
  VAR2p_ScatMnoScat_SOMSST_JJA = reshape( (/VAR2p_ScatMnoScat_SOM_JJA, VAR2p_ScatMnoScat_SST_JJA/), (/nfile*n1,n2/) )
  VAR3p_ScatMnoScat_SOMSST_JJA = reshape( (/VAR3p_ScatMnoScat_SOM_JJA, VAR3p_ScatMnoScat_SST_JJA/), (/nfile*n1,n2/) )

  ;copy_VarAtts(VAR1p_ScatMnoScat_SOM_DJF, VAR1p_ScatMnoScat_SOMSST_DJF)
  ;printVarSummary(VAR1p_ScatMnoScat_SOM_DJF)
  ;printVarSummary(VAR1p_ScatMnoScat_SOMSST_DJF)

  print("write variables for profile plots to "+outfile_profile)

  out_profile->VAR1p_ScatMnoScat_SOM_DJF = VAR1p_ScatMnoScat_SOM_DJF
  out_profile->VAR1p_ScatMnoScat_SST_DJF = VAR1p_ScatMnoScat_SST_DJF
  out_profile->VAR2p_ScatMnoScat_SOM_DJF = VAR2p_ScatMnoScat_SOM_DJF
  out_profile->VAR2p_ScatMnoScat_SST_DJF = VAR2p_ScatMnoScat_SST_DJF
  out_profile->VAR3p_ScatMnoScat_SOM_DJF = VAR3p_ScatMnoScat_SOM_DJF
  out_profile->VAR3p_ScatMnoScat_SST_DJF = VAR3p_ScatMnoScat_SST_DJF

  out_profile->VAR1p_ScatMnoScat_SOM_JJA = VAR1p_ScatMnoScat_SOM_JJA
  out_profile->VAR1p_ScatMnoScat_SST_JJA = VAR1p_ScatMnoScat_SST_JJA
  out_profile->VAR2p_ScatMnoScat_SOM_JJA = VAR2p_ScatMnoScat_SOM_JJA
  out_profile->VAR2p_ScatMnoScat_SST_JJA = VAR2p_ScatMnoScat_SST_JJA
  out_profile->VAR3p_ScatMnoScat_SOM_JJA = VAR3p_ScatMnoScat_SOM_JJA
  out_profile->VAR3p_ScatMnoScat_SST_JJA = VAR3p_ScatMnoScat_SST_JJA

end
