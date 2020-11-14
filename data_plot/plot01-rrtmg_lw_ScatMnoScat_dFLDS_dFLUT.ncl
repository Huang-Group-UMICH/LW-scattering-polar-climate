load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "/home/yihsuan//script/ncl//ignorantNCL.ncl"
load "./ignorantNCL.ncl"
;load "/home/yihsuan//script/ncl//ignorantNCL-colormaps.ncl"
;load "zz-read_parameters.ncl"

begin

;------------------
; general settings
;------------------

  ;*** read from command line ***
  datapath = "./"
  filehead = "data_plot01-bbxx3_all_cld_atm_surf-atm_"
  ;iciwp = "iciwp_0.1"
  ;iciwp = "iciwp_0.02"
  iciwp = "iciwp_0.05"

  wks_name = "fig01-rrtmg_lw_ScatMnoScat_dFLDS_dFLUT-sas_HGH_LOW"+"-"+iciwp
  wks_type = "eps"  ; or x11,pdf,ncgm
  ;wks_type = "png"  ; or x11,pdf,ncgm

  ;wks_name = "test111-rrtmg_lw_ScatMnoScat_dFLDS_dFLUT_all_cld_surf_atm-"+atm+"-"+iciwp
  ;wks_name = "ggtest111-"+filehead+"-"+iciwp

  if (filehead .eq. "bbxx3_all_cld_atm_surf-atm_") then
    title = "Scat - noScat, Sub-Arctic summer, CLDase + SURFe + ATMae"
  end if

  if (filehead .eq. "bbxx1_CLDase_SURFe_ATMx-atm_") then
    title = "Scat - noScat, Sub-Arctic summer, CLDase + SURFe + ATMx"
  end if

  if (filehead .eq. "bbxx2_CLDase_SURFx_ATMx-atm_") then
    title = "Scat - noScat, Sub-Arctic summer, CLDase + SURFx + ATMx"
  end if

  ;filehed = "bbb1-atm_"
  ;filehead = "bbb2_CPKuo-atm_"

  title = "Modified RRTMG offline calculation, Scat - noScat, Sub-Arctic summer"

  atm = "sas"

  filenames = (/ \
                filehead+atm+"-kcld_42-"+iciwp+"-rei_q", \
                filehead+atm+"-kcld_46-"+iciwp+"-rei_q", \
                filehead+atm+"-kcld_48-"+iciwp+"-rei_q" /)

  ;*** iwp=100 g/m2, rei=5., 10., 20., 50., 100.***
  tau_Scat_iwp01 = (/13.61,7.80,4.02,1.57,0.80/)

  if (iciwp.eq."iciwp_0.2") then
    tau_Scat   = tau_Scat_iwp01 * 2.
    string_iwp = "IWP= 200 g/m2"
  end if

  if (iciwp.eq."iciwp_0.1") then
    tau_Scat   = tau_Scat_iwp01
    string_iwp = "IWP= 100 g/m2"
  end if

  if (iciwp.eq."iciwp_0.05") then
    tau_Scat   = tau_Scat_iwp01/2.
    string_iwp = "IWP= 50 g m~S~-2"
  end if

  if (iciwp.eq."iciwp_0.02") then
    tau_Scat   = tau_Scat_iwp01/5.
    string_iwp = "IWP= 20 g m~S~-2"
  end if

  if (iciwp.eq."iciwp_0.01") then
    tau_Scat   = tau_Scat_iwp01/10.
    string_iwp = "IWP= 10 g/m2"
  end if

  if (iciwp.eq."iciwp_0.002") then
    tau_Scat   = tau_Scat_iwp01 / 50.
    string_iwp = "IWP= 2 g/m2"
  end if
  tau_noScat = tau_Scat


  ;title = title+title_case

  tistring_high_dFLUT = "OLR, high cloud"
  tistring_low_dFLUT = "OLR, low cloud"
  tistring_high_dFLDS = "FLDS, high cloud"
  tistring_low_dFLDS = "FLDS, low cloud"

  opt00 = True
    opt00@string_change = "delta"
  tistring_high_dFLUT = string_change(tistring_high_dFLUT,opt00)
  tistring_low_dFLUT = string_change(tistring_low_dFLUT,opt00)
  tistring_high_dFLDS = string_change(tistring_high_dFLDS,opt00)
  tistring_low_dFLDS = string_change(tistring_low_dFLDS,opt00)

  opt00 = True

  nrei = 5  ; # of rei
  ntmq = 5  ; # of tmq
  ncol = 7  ; # of column
  nrow = 6  ; # of row for each case

;-------------
; read files
;-------------

  filenames = datapath+filenames
  nfile = dimsizes(filenames)

  FLDS_Scat = new( (/nfile,nrei,ntmq/) , "float")
  FLDS_Scat!0 = "pcld"
  FLDS_Scat!1 = "rei"
  FLDS_Scat!2 = "tmq"

  FLUT_Scat = FLDS_Scat

  FLDS_noScat = FLDS_Scat
  FLUT_noScat = FLUT_Scat

  ;FLUT_Scat = FLDS_Scat
  ;FLUT_noScat = FLDS_Scat

  pcld = FLDS_Scat(:,0,0)
  npcld = dimsizes(pcld)
  ;rei = FLDS_Scat(0,:,0)
  ;tmq = FLDS_Scat(0,0,:)

  do ff=0,nfile-1
  ;do ff=0,0
    file1 = filenames(ff)

    ;***********
    ; read Scat
    ;***********
    file1_Scat = file1+"-Scat-.txt"
    data1_Scat = asciiread(file1_Scat,(/nrei*ntmq,ncol/),"float")
    ;print(data1_Scat)
    ;printVarSummary(data1_Scat)

    rei_Scat = data1_Scat(::nrow,1)
    ;print(rei_Scat)

    ;tau_Scat = data1_Scat(::nrow,4)
    ;print(rei_Scat)

    tmq_Scat = data1_Scat(0:ntmq-1,3)
    ;print(tmq_Scat)

    pcld(ff) = (/data1_Scat(0,0)/)
    ;print(pcld_Scat)

    ;FLDS_Scat&pcld(ff) = pcld_Scat
    ;printVarSummary(FLDS_Scat)
    FLDS_Scat&rei = rei_Scat
    FLDS_Scat&tmq = tmq_Scat
    FLDS_Scat_1d = data1_Scat(:,5)
    FLDS_Scat_2d = reshape(FLDS_Scat_1d,(/nrei,ntmq/))
    FLDS_Scat(ff,:,:) = (/FLDS_Scat_2d/)

    FLUT_Scat_1d = data1_Scat(:,6)
    FLUT_Scat_2d = reshape(FLUT_Scat_1d,(/nrei,ntmq/))
    FLUT_Scat(ff,:,:) = (/FLUT_Scat_2d/)

    delete(data1_Scat)
    delete(FLDS_Scat_1d)
    delete(FLDS_Scat_2d)
    delete(FLUT_Scat_1d)
    delete(FLUT_Scat_2d)

    ;***********
    ; read noScat
    ;***********
    file1_noScat = file1+"-noScat-.txt"
    data1_noScat = asciiread(file1_noScat,(/nrei*ntmq,ncol/),"float")
    ;print(data1_noScat)

    rei_noScat = data1_noScat(::nrow,1)
    ;print(rei_noScat)

    ;tau_noScat = data1_noScat(::nrow,4)
    ;print(rei_noScat)

    tmq_noScat = data1_noScat(0:ntmq-1,3)
    ;print(tmq_Scat)

    pcld(ff) = (/data1_noScat(0,0)/)
    ;print(pcld_noScat)

    ;FLDS_noScat&pcld(ff) = pcld_noScat
    FLDS_noScat&rei = rei_noScat
    FLDS_noScat&tmq = tmq_noScat
    FLDS_noScat_1d = data1_noScat(:,5)
    FLDS_noScat_2d = reshape(FLDS_noScat_1d,(/nrei,ntmq/))
    FLDS_noScat(ff,:,:) = (/FLDS_noScat_2d/)

    FLUT_noScat_1d = data1_noScat(:,6)
    FLUT_noScat_2d = reshape(FLUT_noScat_1d,(/nrei,ntmq/))
    FLUT_noScat(ff,:,:) = (/FLUT_noScat_2d/)

    delete(data1_noScat)
    delete(FLDS_noScat_1d)
    delete(FLDS_noScat_2d)
    delete(FLUT_noScat_1d)
    delete(FLUT_noScat_2d)

    ;***********
    ; print
    ;***********

    do i=0,nrei-1
    do j=0,ntmq-1
      print("file="+file1_Scat+", rei="+FLDS_Scat&rei(i)+", tmq="+FLDS_Scat&tmq(j)+ \
            ", FLDS="+FLDS_Scat(ff,i,j)+", FLUT="+FLUT_Scat(ff,i,j))
      print("file="+file1_noScat+", rei="+FLDS_noScat&rei(i)+", tmq="+FLDS_noScat&tmq(j)+ \
            ", FLDS="+FLDS_noScat(ff,i,j)+", FLUT="+FLUT_noScat(ff,i,j))
      tt1 = FLDS_Scat(ff,i,j)-FLDS_noScat(ff,i,j)
      tt2 = FLUT_Scat(ff,i,j)-FLUT_noScat(ff,i,j)
      print("cloud: "+pcld(ff)+", Scat - noScat, dFLDS= "+tt1+", dFLUT="+tt2)
      print("--------------------------------")
    end do
    end do

  end do ; end of files

  FLDS_Scat&pcld = (/pcld/)
  FLDS_noScat&pcld = (/pcld/)
  ;printVarSummary(FLDS_Scat)
  ;printVarSummary(FLDS_noScat)

  FLUT_Scat&pcld = (/pcld/)
  FLUT_noScat&pcld = (/pcld/)
  ;printVarSummary(FLUT_Scat)

  ;printVarSummary(FLUT_noScat)
  ;print(FLUT_Scat(0,0,0)+" , "+FLUT_noScat(0,0,0))
  ;print(FLUT_Scat+" , "+FLDS_Scat)

;--------------------
; change coordinate
;--------------------
  rei00 = FLDS_Scat&rei
  dei00 = multiply(rei00, 2.)

  tmq00 = FLDS_Scat&tmq

  rei01 = fspan(0., tofloat(nrei-1), nrei)
  tmq01 = fspan(0., tofloat(ntmq-1), ntmq)

  FLDS_Scat&rei = rei01
  FLDS_Scat&tmq = tmq01
  FLUT_Scat&rei = rei01
  FLUT_Scat&tmq = tmq01

  FLDS_noScat&rei = rei01
  FLDS_noScat&tmq = tmq01
  FLUT_noScat&rei = rei01
  FLUT_noScat&tmq = tmq01

  FLDS_ScatMnoScat = diff(FLDS_Scat, FLDS_noScat)
  FLUT_ScatMnoScat = diff(FLUT_Scat, FLUT_noScat)

;printVarSummary(FLDS_Scat)
;print(rei00+" , "+tau_Scat)

;-------
; plot
;------

  ;*** share res ***
  res = True
  res@gsnFrame = False
  res@gsnDraw = False

  res@trXMinF = -0.5
  res@trXMaxF = res@trXMinF+tofloat(ntmq) 
  res@trYMinF = -0.5
  res@trYMaxF = res@trYMinF+tofloat(nrei) 

  res@cnFillOn = True
  res@cnFillMode = "RasterFill"

  res@cnLinesOn = False
  res@cnLineLabelsOn = False
  res@cnInfoLabelOn = False

  res@lbOrientation = "vertical"

  res@tmXBMode = "Explicit"
  res@tmXBValues = tmq01
  string_xx = "~C~("+sprintf("%3.2f",(/0.25,0.5,1,1.25,1.5/))+"x)"
  res@tmXBLabels = sprintf(" %4.2f" , tmq00)+string_xx

  res@tmYLMode = "Explicit"
  res@tmYLValues = rei01
  ;yllabels = sprintf("%3.0f" , rei00) + "~C~(~F8~t~F21~="+sprintf("%4.2f" , tau_Scat)+")"
  ;yllabels = sprintf("%3.0f" , rei00) + "~C~(~F8~t~F21~"+"~B~820-980~N~="+sprintf("%4.2f" , tau_Scat)+")"
  yllabels = sprintf("%3.0f" , dei00) + "~C~(~F8~t~F21~"+"="+sprintf("%4.2f" , tau_Scat)+")"
  res@tmYLLabels = yllabels

  res@tiXAxisString = "Total column water vapor (mm)"
  res@tiYAxisString = "Ice effective diameter (microns)"
  res@gsnRightString = "W m~S~-2"
  ;print(res)

  resp = True
  ;resp@txString = title
  resp@gsnPanelFigureStrings= (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)"/)
  resp@gsnPanelFigureStringsBackgroundFillColor = "transparent"
  resp@gsnPanelFigureStringsJust = "topright"
  resp@gsnPanelYWhiteSpacePercent = 3 
 
  ;*** FLDS ***
  res_FLDS = res
  res_FLDS@cnFillPalette = "WhiteYellowOrangeRed"
  ;res_FLDS@gsnLeftString = "Surf Dw LW flux diff"

  res_FLDS@cnLevelSelectionMode = "ManualLevels"
  res_FLDS@cnMinLevelValF = 0.5
  res_FLDS@cnMaxLevelValF = 4.
  res_FLDS@cnLevelSpacingF = 0.25

  res_FLDS@lbLabelStride = 2

  ;tt1_FLDS = reshape(FLDS_ScatMnoScat,(/npcld*nrei,ntmq/))
  ;opt_cncolors_FLDS = True
  ;cnColors_PosNeg(tt1_FLDS,res_FLDS,opt_cncolors_FLDS)
  ;print(res_FLDS)

  ;*** FLUT ***
  res_FLUT = res
  ;res_FLUT@gsnLeftString = "Outgoing LW radiation diff"

  ;res_FLUT@cnFillPalette = "CBR_wet"
  ;res_FLUT@cnFillPalette = "MPL_GnBu"
  ;res_FLUT@cnFillPalette = "AMWG"
  ;res_FLUT@gsnSpreadColorEnd = 1
  ;res_FLUT@gsnSpreadColorStart = -1

  cmap_name = "MPL_GnBu"
  cmap_values = read_colormap_file(cmap_name)
  cmap_values = cmap_values(::-1,:) ; reverse the color map 
  ;res_FLUT@cnFillColors = cmap_values
  res_FLUT@cnFillPalette = cmap_values

  res_FLUT@cnLevelSelectionMode = "ManualLevels"
  res_FLUT@cnMinLevelValF = -12.
  res_FLUT@cnMaxLevelValF = -1.
  res_FLUT@cnLevelSpacingF = 0.5

  ;tt1_FLUT = reshape(FLUT_ScatMnoScat,(/npcld*nrei,ntmq/))
  ;opt_cncolors_FLUT = True
  ;opt_cncolors_FLUT@max_levels = 30
  ;cnColors_PosNeg(tt1_FLUT,res_FLUT,opt_cncolors_FLUT)
  ;print(res_FLUT)

;-------
; plot
;-------
  ;plot = new(npcld*2, "graphic")
  ;nplot = dimsizes(plot)
  wks = gsn_open_wks(wks_type, wks_name)

  ;printVarSummary(FLDS_ScatMnoScat)
  ;print(FLDS_ScatMnoScat)

  ;http://eodg.atm.ox.ac.uk/RFM/atm/
  ; High cloud: 7-8 km
  ; Low cloud:  1-2 km
  tiMainString_FLUT = (/"High cloud", "dOLR, Middle cloud","Low cloud"/)
  tiMainString_FLDS = (/"High cloud", "dFLDS, Middle cloud","Low cloud"/)
  ;tiMainString_FLUT = (/"dOLR, High cloud (386 hPa)", "dOLR, Middle cloud","dOLR, Low cloud (844 hPa)"/)
  ;tiMainString_FLDS = (/"dFLDS, High cloud (386 hPa)", "dFLDS, Middle cloud","dFLDS, Low cloud (844 hPa)"/)

  ;tiMainString_FLUT = tiMainString_FLUT+", "+string_iwp
  ;tiMainString_FLDS = tiMainString_FLDS+", "+string_iwp
  ;print(tiMainString_FLUT)
    
    pp = 0
    res_FLUT@tiMainString = tistring_high_dFLUT
    plot1 = gsn_csm_contour(wks, FLUT_ScatMnoScat(pp,:,:), res_FLUT)

    pp = 2
    res_FLUT@tiMainString = tistring_low_dFLUT
    plot2 = gsn_csm_contour(wks, FLUT_ScatMnoScat(pp,:,:), res_FLUT)

    pp1 = 0
    res_FLDS@tiMainString = tistring_high_dFLDS
    plot3 = gsn_csm_contour(wks, FLDS_ScatMnoScat(pp1,:,:), res_FLDS)

    pp1 = 2
    res_FLDS@tiMainString = tistring_low_dFLDS
    plot4 = gsn_csm_contour(wks, FLDS_ScatMnoScat(pp1,:,:), res_FLDS)

    pp=0
    ;print(min(FLUT_ScatMnoScat(pp,:,:)))
    ;print(max(FLUT_ScatMnoScat(pp,:,:)))

    ;print(min(FLDS_ScatMnoScat(pp,:,:)))
    ;print(max(FLDS_ScatMnoScat(pp,:,:)))

;  do pp=0,npcld-1
;    res_FLUT@tiMainString = tiMainString_FLUT(pp)
;    plot(pp) = gsn_csm_contour(wks, FLUT_ScatMnoScat(pp,:,:), res_FLUT)
;  end do
;
;  do pp=npcld,nplot-1
;    pp1 = pp-npcld
;    res_FLDS@tiMainString = tiMainString_FLDS(pp1)
;    plot(pp) = gsn_csm_contour(wks, FLDS_ScatMnoScat(pp1,:,:), res_FLDS)
;  end do




  lnres = True
  lnres@gsLineColor = (/1./)
  lnres@gsLineThicknessF = 1.
  lnres@gsLineDashPattern = 2

  ;saw_pw = 4.55
  ;sas_pw = 24.52
  
  y1 = -1.
  y2 = 5.

  x1_sas = 2.
  x2_sas = x1_sas

  x1_saw = -0.13
  x2_saw = x1_saw

  ;do i=0,3
  ;  ln1 = gsn_add_polyline(wks,plot1, (/x1_sas,x2_sas/), (/y1,y2/), lnres)
  ;end do

  plot = new(4, "graphic")
  plot(0) = plot1
  plot(1) = plot2
  plot(2) = plot3
  plot(3) = plot4

  lns = new(1000, "graphic")
  do i=0,3
    lns(j) = gsn_add_polyline(wks,plot(i), (/x1_sas,x2_sas/), (/y1,y2/), lnres)
    lns(j+1) = gsn_add_polyline(wks,plot(i), (/x1_saw,x2_saw/), (/y1,y2/), lnres)
    j = j+2
  end do

  row1    = 2  ; row number of panel plot
  column1 = 2  ; column number of panel plot
  gsn_panel(wks,plot,(/row1,column1/),resp)

  opt11 = True
  trim_fig(wks_type, wks_name, opt11)

end
