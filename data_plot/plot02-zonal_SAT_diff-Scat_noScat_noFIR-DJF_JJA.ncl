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
; read start
;--------
  yy_string = "yy06_35"

  VAR1_name = "TREFHT"
  VAR2_name = "FLDS"

  ;casename_SOM_Scat = "c10_en01-cesm111-E2000_rrtmg_mc6_rtr2"
  ;casename_SOM_noScat = "c11_en01-cesm111-E2000_rrtmg_mc6_rtr2_noScat"
  ;casename_SST_Scat = "g01_en01-cesm111-F2000-mc6_rtr2_Scat"
  ;casename_SST_noScat = "g02_en01-cesm111-F2000-mc6_rtr2_noScat"

  txstring = ""

  wks_name_ScatAnoScat = "fig02-zonal_SAT_diff-Scat_noScat_noFIR-DJF_JJA"
  wks_type = "eps"  ; or x11,pdf,ncgm

  txstring_head = "CESM/SOM & AMIP, DJF & JJA, Scat-noScat"
  plotname_head = "xxfig_E2000_F2000"

  ;choice_landfrac = "land_only"
  ;choice_landfrac = "non_land"
  choice_landfrac = "all_grids"

  datapath = "./data_plot02-"

  casename_SOM_Scat = (/"c10_en01-cesm111-E2000_rrtmg_mc6_rtr2", \
                        "c10_en02-cesm111-E2000_rrtmg_mc6_rtr2", \
                        "c10_en03-cesm111-E2000_rrtmg_mc6_rtr2", \
                        "c10-cesm111-E2000_rrtmg_mc6_rtr2"/)

  casename_SOM_noFIR = (/"q01_en01-cesm111-E2000_rrtmg_mc6_rtr2_ScatFIRoff", \
                         "q01_en02-cesm111-E2000_rrtmg_mc6_rtr2_ScatFIRoff", \
                         "q01_en03-cesm111-E2000_rrtmg_mc6_rtr2_ScatFIRoff", \
                         "q01-cesm111-E2000_rrtmg_mc6_rtr2_ScatFIRoff" /)

  casename_SOM_noScat = (/"c11_en01-cesm111-E2000_rrtmg_mc6_rtr2_noScat", \
                        "c11_en02-cesm111-E2000_rrtmg_mc6_rtr2_noScat", \
                        "c11_en03-cesm111-E2000_rrtmg_mc6_rtr2_noScat", \
                        "c11-cesm111-E2000_rrtmg_mc6_rtr2_noScat"/)

  casename_SST_Scat = (/"g01_en01-cesm111-F2000-mc6_rtr2_Scat", \
                        "g01_en02-cesm111-F2000-mc6_rtr2_Scat", \
                        "g01_en03-cesm111-F2000-mc6_rtr2_Scat", \
                        "g01-cesm111-F2000-mc6_rtr2_Scat"/)

  casename_SST_noScat = (/"g02_en01-cesm111-F2000-mc6_rtr2_noScat", \
                        "g02_en02-cesm111-F2000-mc6_rtr2_noScat", \
                        "g02_en03-cesm111-F2000-mc6_rtr2_noScat", \
                        "g02-cesm111-F2000-mc6_rtr2_noScat"/)

  nmember = dimsizes(casename_SOM_Scat)
  
  plotname_step1 = "DJF"
  plotname_step2 = "JJA"

;---------
; read end
;--------
  plotname_head = plotname_head ;+"_"+choice_landfrac


  filehead_SOM_Scat   = datapath+casename_SOM_Scat+".cam.h0."+yy_string+"-climo_"
  filehead_SOM_noScat   = datapath+casename_SOM_noScat+".cam.h0."+yy_string+"-climo_"
  filehead_SST_Scat   = datapath+casename_SST_Scat+".cam.h0."+yy_string+"-climo_"
  filehead_SST_noScat   = datapath+casename_SST_noScat+".cam.h0."+yy_string+"-climo_"
  filehead_SOM_noFIR   = datapath+casename_SOM_noFIR+".cam.h0."+yy_string+"-climo_"

  ;txstring = txstring_head+", "+yy_string+", "+plotname_step1+"_"+plotname_step2+", "+choice_landfrac

  plotname_type = "xy_zonal"
  plotname_file = "SOM_SST"
  ;plotname_step = yy_string+"_DJF_JJA"
  plotname_step = yy_string+"_"+plotname_step1+"_"+plotname_step2
 
  plotname_var  = VAR1_name+"_"+VAR2_name+"_"+choice_landfrac


  suffixes_ScatAnoScat = (/"  Scat - noScat, SOM", "  Scat - noScat, SST", "  noFIR - noScat, SOM"/)
  plotname_var  = VAR1_name+"_"+VAR2_name+"_diff"
  ;wks_name_ScatAnoScat = plotname_head +"-"+ plotname_type +"-"+ plotname_file +"-"+ plotname_var +"-"+ plotname_step +"-"
  wks_ScatAnoScat = gsn_open_wks(wks_type, wks_name_ScatAnoScat)

  suffixes_noScat = (/"  noScat, SOM", "  noScat, SST"/)
  plotname_var  = VAR1_name+"_"+VAR2_name+"_noScat"
  wks_name_noScat = plotname_head +"-"+ plotname_type +"-"+ plotname_file +"-"+ plotname_var +"-"+ plotname_step +"-"
  ;wks_noScat = gsn_open_wks(wks_type, wks_name_noScat)

  opt_ttest = True
  siglvl    = 0.05

  ;suffixes = (/"Scat_SOM - noScat_SOM","noFIR_SOM - noScat_SOM", "Scat_SST - noScat_SST"/)
  ;suffixes = (/"  Scat - noScat, SOM","  noFIR - noScat, SOM", "  Scat - noScat, SST"/)

  plotname_file = "global-TREFHT"
  ;plotname_step = "DJF"
  infilename_SOM_Scat_DJF_global   = filehead_SOM_Scat+plotname_step1+"-"+plotname_file+".nc"
  infilename_SOM_noScat_DJF_global = filehead_SOM_noScat+plotname_step1+"-"+plotname_file+".nc"
  infilename_SOM_noFIR_DJF_global = filehead_SOM_noFIR+plotname_step1+"-"+plotname_file+".nc"
  infilename_SST_Scat_DJF_global   = filehead_SST_Scat+plotname_step1+"-"+plotname_file+".nc"  
  infilename_SST_noScat_DJF_global = filehead_SST_noScat+plotname_step1+"-"+plotname_file+".nc"

  plotname_file = "global-TREFHT"
  ;plotname_step = "JJA"
  infilename_SOM_Scat_JJA_global   = filehead_SOM_Scat+plotname_step2+"-"+plotname_file+".nc"
  infilename_SOM_noScat_JJA_global = filehead_SOM_noScat+plotname_step2+"-"+plotname_file+".nc"
  infilename_SOM_noFIR_JJA_global = filehead_SOM_noFIR+plotname_step2+"-"+plotname_file+".nc"
  infilename_SST_Scat_JJA_global   = filehead_SST_Scat+plotname_step2+"-"+plotname_file+".nc"  
  infilename_SST_noScat_JJA_global = filehead_SST_noScat+plotname_step2+"-"+plotname_file+".nc"

  opt_1Ddata_zonal = True
  opt_1Ddata_zonal@return_type = "ftji2j_fj_avgti"      ; input: CESM (files,time,lat,lon), return: var_xx=(lat), var_yy=(files,lat), averaged over time and lon
  opt_1Ddata_zonal@resolution = "cesm_1.9x2.5_"+plotname_file
  ;opt_1Ddata_zonal@select_landfrac = True
  ;opt_1Ddata_zonal@landfrac = choice_landfrac

;-----------------
; read variables
;-----------------

  ;*** DJF ***
  VAR1_SOM_Scat_DJF0 = convert_files_vars(infilename_SOM_Scat_DJF_global, VAR1_name, opt_1Ddata_zonal)
  VAR1_SOM_Scat_DJF0!2 = "latitude"
  VAR1_fj_SOM_Scat_DJF = VAR1_SOM_Scat_DJF0(:,0,:)
  ;printVarSummary(VAR1_fj_SOM_Scat_DJF)

  VAR1_SOM_noScat_DJF0 = convert_files_vars(infilename_SOM_noScat_DJF_global, VAR1_name, opt_1Ddata_zonal)
  VAR1_SOM_noScat_DJF0!2 = "latitude"
  VAR1_fj_SOM_noScat_DJF = VAR1_SOM_noScat_DJF0(:,0,:)
  ;printVarSummary(VAR1_fj_SOM_noScat_DJF)

  VAR1_SOM_noFIR_DJF0 = convert_files_vars(infilename_SOM_noFIR_DJF_global, VAR1_name, opt_1Ddata_zonal)
  VAR1_SOM_noFIR_DJF0!2 = "latitude"
  VAR1_fj_SOM_noFIR_DJF = VAR1_SOM_noFIR_DJF0(:,0,:)
  ;printVarSummary(VAR1_fj_SOM_noFIR_DJF)

  VAR1_SST_Scat_DJF0 = convert_files_vars(infilename_SST_Scat_DJF_global, VAR1_name, opt_1Ddata_zonal)
  VAR1_SST_Scat_DJF0!2 = "latitude"
  VAR1_fj_SST_Scat_DJF = VAR1_SST_Scat_DJF0(:,0,:)
  ;printVarSummary(VAR1_fj_SST_Scat_DJF)

  VAR1_SST_noScat_DJF0 = convert_files_vars(infilename_SST_noScat_DJF_global, VAR1_name, opt_1Ddata_zonal)
  VAR1_SST_noScat_DJF0!2 = "latitude"
  VAR1_fj_SST_noScat_DJF = VAR1_SST_noScat_DJF0(:,0,:)
  ;printVarSummary(VAR1_fj_SST_noScat_DJF)
  
  ;*** JJA ***
  VAR1_SOM_Scat_JJA0 = convert_files_vars(infilename_SOM_Scat_JJA_global, VAR1_name, opt_1Ddata_zonal)
  VAR1_SOM_Scat_JJA0!2 = "latitude"
  VAR1_fj_SOM_Scat_JJA = VAR1_SOM_Scat_JJA0(:,0,:)
  ;printVarSummary(VAR1_fj_SOM_Scat_JJA)

  VAR1_SOM_noScat_JJA0 = convert_files_vars(infilename_SOM_noScat_JJA_global, VAR1_name, opt_1Ddata_zonal)
  VAR1_SOM_noScat_JJA0!2 = "latitude"
  VAR1_fj_SOM_noScat_JJA = VAR1_SOM_noScat_JJA0(:,0,:)
  ;printVarSummary(VAR1_fj_SOM_noScat_JJA)

  VAR1_SOM_noFIR_JJA0 = convert_files_vars(infilename_SOM_noFIR_JJA_global, VAR1_name, opt_1Ddata_zonal)
  VAR1_SOM_noFIR_JJA0!2 = "latitude"
  VAR1_fj_SOM_noFIR_JJA = VAR1_SOM_noFIR_JJA0(:,0,:)
  ;printVarSummary(VAR1_fj_SOM_noFIR_JJA)

  VAR1_SST_Scat_JJA0 = convert_files_vars(infilename_SST_Scat_JJA_global, VAR1_name, opt_1Ddata_zonal)
  VAR1_SST_Scat_JJA0!2 = "latitude"
  VAR1_fj_SST_Scat_JJA = VAR1_SST_Scat_JJA0(:,0,:)
  ;printVarSummary(VAR1_fj_SST_Scat_JJA)

  VAR1_SST_noScat_JJA0 = convert_files_vars(infilename_SST_noScat_JJA_global, VAR1_name, opt_1Ddata_zonal)
  VAR1_SST_noScat_JJA0!2 = "latitude"
  VAR1_fj_SST_noScat_JJA = VAR1_SST_noScat_JJA0(:,0,:)
  ;printVarSummary(VAR1_fj_SST_noScat_JJA)

  VAR1_fj_SOM_ScatMnoScat_DJF = diff(VAR1_fj_SOM_Scat_DJF,VAR1_fj_SOM_noScat_DJF)
  VAR1_fj_SOM_ScatMnoScat_JJA = diff(VAR1_fj_SOM_Scat_JJA,VAR1_fj_SOM_noScat_JJA)
  VAR1_fj_SST_ScatMnoScat_DJF = diff(VAR1_fj_SST_Scat_DJF,VAR1_fj_SST_noScat_DJF)
  VAR1_fj_SST_ScatMnoScat_JJA = diff(VAR1_fj_SST_Scat_JJA,VAR1_fj_SST_noScat_JJA)
  VAR1_fj_SOM_noFIRMnoScat_DJF = diff(VAR1_fj_SOM_noFIR_DJF,VAR1_fj_SOM_noScat_DJF)
  VAR1_fj_SOM_noFIRMnoScat_JJA = diff(VAR1_fj_SOM_noFIR_JJA,VAR1_fj_SOM_noScat_JJA)

  nn = dimsizes(VAR1_fj_SOM_noScat_JJA)
  n1 = nn(0)
  n2 = nn(1)

  VAR1_fj_ScatAnoScat_DJF = reshape( (/VAR1_fj_SOM_ScatMnoScat_DJF,VAR1_fj_SST_ScatMnoScat_DJF,VAR1_fj_SOM_noFIRMnoScat_DJF/), (/3*n1,n2/))
  VAR1_fj_noScat_DJF = reshape( (/VAR1_fj_SOM_noScat_DJF,VAR1_fj_SST_noScat_DJF/), (/2*n1,n2/))

  ;printVarSummary(VAR1_fj_ScatAnoScat_DJF)
  ;printVarSummary(VAR1_fj_noScat_DJF)

  VAR1_fj_ScatAnoScat_JJA = reshape( (/VAR1_fj_SOM_ScatMnoScat_JJA,VAR1_fj_SST_ScatMnoScat_JJA,VAR1_fj_SOM_noFIRMnoScat_JJA/), (/3*n1,n2/))
  VAR1_fj_noScat_JJA = reshape( (/VAR1_fj_SOM_noScat_JJA,VAR1_fj_SST_noScat_JJA/), (/2*n1,n2/))

  ;res00 = True
  ;opt00 = True
  ;printVarSummary(VAR1_fj_ScatMnoScat_DJF)
  ;print(res00)
  ;exit 

  lat_fj_Scat_DJF = VAR1_SOM_Scat_DJF0&latitude
  lat_fj_Scat_JJA = VAR1_SOM_Scat_JJA0&latitude

  lat_fj_Scat_DJF@units = "degrees_north"
  lat_fj_Scat_JJA@units = "degrees_north"

;--------
; plot
;--------
  res  = True
  res@gsnFrame = False
  res@gsnDraw = False
  res@pmLegendHeightF = 0.15
  ;res@pmLegendHeightF = 0.22
  ;res@pmLegendHeightF = 0.24

  res@xyLineColors = (/"red","blue","goldenrod"/)
  res@pmLegendDisplayMode = "NoCreate"

  res@xyMonoDashPattern = True
  res@xyDashPattern = 0

  res@tiMainString = " "

  res_zonal = res
  res_zonal@gsnYRefLine = 0.
  res_zonal@tiXAxisString = "Latitude"

  res_zonal@trYMinF = -1.5
  res_zonal@trYMaxF = 2.5

  res_zonal@tmYLMode = "Manual"
  res_zonal@tmYLTickStartF = res_zonal@trYMinF
  res_zonal@tmYLTickEndF = res_zonal@trYMaxF
  res_zonal@tmYLTickSpacingF = 0.5

  resp = True    ; panel res
  resp@txString = txstring


  polyres = True
  polyres@gsMarkerIndex = 16   ; solid dot
  polyres@gsMarkerSizeF = 5. 

  opt_xy = True                         ; option for ignorant_xy
  opt_xy@legend_mean = False   ; turn on/off mean values on legends
  opt_xy@option_group = True     ; if grouping the first dimension
  opt_xy@nmember_per_group = nmember   ; number of member per group 
  ;opt_xy@XYFill = True
  ;opt_xy@XYFill = False
  opt_xy@XYStyle = "fill"
  ;opt_xy@XYStyle = "dash"

  res_VAR1 = res_zonal
  change_var_atts(VAR1_name, VAR1_fj_ScatAnoScat_DJF, res_VAR1, True)
  ;res_VAR1@gsnRightString = VAR1_fj_SST_noScat_JJA@units
  ;res_VAR1@gsnLeftString = "Surf air temp diff"
  ;res_VAR1@tiYAxisString = res_VAR1@gsnLeftString+" ("+res_VAR1@gsnRightString+")"
  res_VAR1@tiYAxisString = "~F33~D~F21~SAT (K)"
  res_VAR1@gsnRightString = ""
  res_VAR1@gsnLeftString = ""
 ;print(res_VAR1)

  ;res_VAR1@tiMainString = VAR1_name+", "+plotname_step1
  res_VAR1@tiMainString = plotname_step1
  ;res_VAR1@tiMainString = "~F22~DJF"
  plot_VAR1_ScatAnoScat_DJF = ignorant_xy_ensemble(wks_ScatAnoScat, lat_fj_Scat_DJF , VAR1_fj_ScatAnoScat_DJF, suffixes_ScatAnoScat, res_VAR1, opt_xy)
  ;plot_VAR1_noScat_DJF = ignorant_xy_ensemble(wks_noScat, lat_fj_Scat_DJF , VAR1_fj_noScat_DJF, suffixes_noScat, res_VAR1, opt_xy)

  ;res_VAR1@tiMainString = VAR1_name+", "+plotname_step2
  res_VAR1@tiMainString = plotname_step2
  ;res_VAR1@tiMainString = "~F22~JJA"
  plot_VAR1_ScatAnoScat_JJA = ignorant_xy_ensemble(wks_ScatAnoScat, lat_fj_Scat_JJA , VAR1_fj_ScatAnoScat_JJA, suffixes_ScatAnoScat, res_VAR1, opt_xy)
  ;plot_VAR1_noScat_JJA = ignorant_xy_ensemble(wks_noScat, lat_fj_Scat_JJA , VAR1_fj_noScat_JJA, suffixes_noScat, res_VAR1, opt_xy)

;  ;----------------------------
;  ; overlay significan levels
;  ;----------------------------
;  plot_polymarker = new(10000, graphic)  ; id for each p-values (polymarker)
;  polyres1 = polyres
;  iip = 0  ; id number of each polymarker
;
;  ;*** zonal VAR1, DJF & JJA ***
;  npp1 = dimsizes(VAR1_ftj_noScat_DJF(:,0,0))
;  npp2 = dimsizes(VAR1_ftj_noScat_DJF(0,0,:))
;  pvalues_VAR1_ScatAnoScat_DJF = VAR1_ftj_noScat_DJF(:,0,:)
;  pvalues_VAR1_ScatAnoScat_DJF = 1.
;  pvalues_VAR1_ScatAnoScat_JJA = VAR1_ftj_noScat_JJA(:,0,:)
;  pvalues_VAR1_ScatAnoScat_JJA = 1.
;
;      do i=0,npp1-1
;        pvalues_VAR1_ScatAnoScat_DJF(i,:) = ignorant_ttest(VAR1_ftj_Scat_DJF(i,:,:),VAR1_ftj_noScat_DJF(i,:,:),opt_ttest)
;        pvalues_VAR1_ScatAnoScat_JJA(i,:) = ignorant_ttest(VAR1_ftj_Scat_JJA(i,:,:),VAR1_ftj_noScat_JJA(i,:,:),opt_ttest)
;      end do
;
;      pval_xx_DJF = lat_fj_Scat_DJF
;      pval_yy_DJF = VAR1_fj_ScatAnoScat_DJF
;      pval_xx_JJA = lat_fj_Scat_JJA
;      pval_yy_JJA = VAR1_fj_ScatAnoScat_JJA
;
;          do i=0,npp1-1
;            polyres1@gsMarkerColor = res@xyLineColors(i)
;
;          do j=0,npp2-1
;            if (.not.ismissing(pvalues_VAR1_ScatAnoScat_DJF(i,j)) .and. pvalues_VAR1_ScatAnoScat_DJF(i,j) .le. siglvl) then
;              plot_polymarker(iip) = gsn_add_polymarker(wks, plot_VAR1_ScatAnoScat_DJF, \
;                                            pval_xx_DJF(j),pval_yy_DJF(i,j),polyres1)
;              iip = iip+1
;            end if
;            if (.not.ismissing(pvalues_VAR1_ScatAnoScat_JJA(i,j)) .and. pvalues_VAR1_ScatAnoScat_JJA(i,j) .le. siglvl) then
;              plot_polymarker(iip) = gsn_add_polymarker(wks, plot_VAR1_ScatAnoScat_JJA, \
;                                            pval_xx_JJA(j),pval_yy_JJA(i,j),polyres1)
;              iip = iip+1
;            end if
;          end do
;          end do
; 
; ;*** zonal VAR2, DJF & JJA ***
;  delete(pval_xx_DJF)
;  delete(pval_yy_DJF)
;  delete(pval_xx_JJA)
;  delete(pval_yy_JJA)
;
;  npp1 = dimsizes(VAR2_ftj_noScat_DJF(:,0,0))
;  npp2 = dimsizes(VAR2_ftj_noScat_DJF(0,0,:))
;  pvalues_VAR2_ScatAnoScat_DJF = VAR2_ftj_noScat_DJF(:,0,:)
;  pvalues_VAR2_ScatAnoScat_DJF = 1.
;  pvalues_VAR2_ScatAnoScat_JJA = VAR2_ftj_noScat_JJA(:,0,:)
;  pvalues_VAR2_ScatAnoScat_JJA = 1.
;
;      do i=0,npp1-1
;        pvalues_VAR2_ScatAnoScat_DJF(i,:) = ignorant_ttest(VAR2_ftj_Scat_DJF(i,:,:),VAR2_ftj_noScat_DJF(i,:,:),opt_ttest)
;        pvalues_VAR2_ScatAnoScat_JJA(i,:) = ignorant_ttest(VAR2_ftj_Scat_JJA(i,:,:),VAR2_ftj_noScat_JJA(i,:,:),opt_ttest)
;      end do
;
;      pval_xx_DJF = lat_fj_Scat_DJF
;      pval_yy_DJF = VAR2_fj_ScatAnoScat_DJF
;      pval_xx_JJA = lat_fj_Scat_JJA
;      pval_yy_JJA = VAR2_fj_ScatAnoScat_JJA
;
;          do i=0,npp1-1
;            polyres1@gsMarkerColor = res@xyLineColors(i)
;
;          do j=0,npp2-1
;            if (.not.ismissing(pvalues_VAR2_ScatAnoScat_DJF(i,j)) .and. pvalues_VAR2_ScatAnoScat_DJF(i,j) .le. siglvl) then
;              plot_polymarker(iip) = gsn_add_polymarker(wks, plot_VAR2_ScatAnoScat_DJF, \
;                                            pval_xx_DJF(j),pval_yy_DJF(i,j),polyres1)
;              iip = iip+1
;            end if
;            if (.not.ismissing(pvalues_VAR2_ScatAnoScat_JJA(i,j)) .and. pvalues_VAR2_ScatAnoScat_JJA(i,j) .le. siglvl) then
;              plot_polymarker(iip) = gsn_add_polymarker(wks, plot_VAR2_ScatAnoScat_JJA, \
;                                            pval_xx_JJA(j),pval_yy_JJA(i,j),polyres1)
;              iip = iip+1
;            end if
;          end do
;          end do

;--------------
; draw legend
;--------------
;
  lgres = True
  lgres@lgMonoLineColor = False
  lgres@lgLineColors = res@xyLineColors(::-1)
  lgres@lgMonoDashIndex = res@xyMonoDashPattern
  lgres@lgDashIndex = res@xyDashPattern
  ;lgres@lgDashIndexes = res@xyDashPatterns
  lgres@lgMonoLineThickness = True
  lgres@lgLineThicknessF = 5. 

  lgres@lgAutoManage = False
  lgres@lgLabelFontHeightF = 0.2

  lgres@vpWidthF = 0.2
  lgres@vpHeightF = 0.1

  lgres@lgPerimOn = False
  ;lgres@lgPerimOn = True

  label_names_ScatAnoScat = suffixes_ScatAnoScat(::-1)
  label_names_noScat = suffixes_noScat(::-1)
  ;drawNDCGrid(wks_ScatAnoScat)

  lgx1 = 0.18
  lgx2 = 0.7
  lgy1 = 0.62
  lgy2 = lgy1

  lg_xx_ScatAnoScat = (/lgx1,lgx2/)
  lg_yy_ScatAnoScat = (/lgy1,lgy2/)

  nn1 = dimsizes(lg_xx_ScatAnoScat)

  do i=0,nn1-1
    lg_x = lg_xx_ScatAnoScat(i)
    lg_y = lg_yy_ScatAnoScat(i)
    gsn_legend_ndc(wks_ScatAnoScat,dimsizes(label_names_ScatAnoScat),label_names_ScatAnoScat,lg_x,lg_y,lgres)

    lg_x = lg_xx_ScatAnoScat(i)
    lg_y = lg_yy_ScatAnoScat(i)
    ;gsn_legend_ndc(wks_noScat,dimsizes(label_names_noScat),label_names_noScat,lg_x,lg_y,lgres)
  end do

;--------
; panel
;--------
  plot_all = new(2,graphic)
  row1 = 1
  column1 = 2

  plot_all(0) = plot_VAR1_ScatAnoScat_DJF
  plot_all(1) = plot_VAR1_ScatAnoScat_JJA

  gsn_panel(wks_ScatAnoScat, plot_all, (/row1,column1/), resp)

  opt_trim = True
  trim_fig(wks_type, wks_name_ScatAnoScat, opt_trim)

  ;plot_all(0) = plot_VAR1_noScat_DJF
  ;plot_all(1) = plot_VAR1_noScat_JJA
  ;gsn_panel(wks_noScat, plot_all, (/row1,column1/), resp)

  ;opt1 = True
  ;trim_fig(wks_type, wks_name, opt1)

end
