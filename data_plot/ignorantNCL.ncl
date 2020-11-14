;*****************************************************
;  NCL program developed by ignorant
;  
;  History:
;           before 2010/10/25  ignorant_contour,ignorant_nc_contour
;           2010/10/25         ignorant_HMLcloud, first release
;           2010/11/07         ignorant_LWP, first release
;           2010/12/27         ignorant_write_nc, first release
;                              test OK in float varible, character & string varible is not test
;           2010/12/31         CHANGE *ignorant_write_nc*, [opt] statement
;                              old: opt=0,1,2 ; new: opt=-1,-2,N
;                              for writing value in file more flexibility, see discription of it
;
;           2011/01/04         CHANGE *ignorant_write_nc*, [opt] statement
;                              old: opt=-1,-2,N ; new: opt=-1,-2,-3,N
;                              opt=-3 allow automatically create non-pre-existing varible and dimensions
;
;           2011/08/29         ADD function "ignorant_SelectColorMap"
;
;           2011/09/30         ADD function "ignorant_merge_colormaps"
; 
;           2011/10/31         ADD function "ignorant_xy"
;           2011/12/02         ADD function "ignorant_combine_colormaps"
;           2011/12/07         CHANGE function name "ignorant_SelectColorMap" 
;                              to "ignorant_select_colormap"
;           2013/04/16         ADD function "ignorant_get_var_files"
;           2013/11/08         note that function "ignorant_select_colormap" MUST call 
;                              before call "gsn_open_wks", otherwise it would show error message:
;                              "fatal:_NclCreateHLUObjOp: Parent value is not an hlu object or is undefined"
;                              in RCEC nuwa machine
;           2013/12/10         ADD a choice in function "ignorant_select_colormap", let user can select all colors
;                              in given color map.
;           2014/04/03         CHANGE function "ignorant_select_colormap", add extra input argument to include/exclude
;                              foreground and backgroung colors
;           2014/04/24         ADD funtion "gsn_csm_contour_map_avg"
;           2014/04/25         ADD function "diff"
;           2014/09/25         ADD function "multiply"
;           2014/09/29         ADD function "divide"
;           2016/02/20         ADD function "lat_wgt_avg"
;           2016/02/21         REWRITE function "lat_wgt_avg", to avoid missing values in data_in
;           2016/02/20         ADD function "plus"
;           2016/03/20         ADD function "binmap"
;           2016/03/29         ADD attribute "pdf" in function binmap
;           2016/05/10         CHANGE function binmap sorting code to speed up
;           2016/10/07         ADD function "taylor_diagram"
;           2016/10/21         MODIFY cosine latitude weighted method in function "gsn_csm_contour_map_avg"
;                              MODIFY function ignorant_select_colormap by using read_colormap_file
;           2016/11/03         FIX a bug at "ignorant_select_colormap" by adding fore/background colors
;           2017/02/17         REWRITE function "lat_wgt_avg" by applying NCL function
;           2017/08/30         CHANGE adding mean values in "ignorant_xy" legends, by setting opt@legend_mean = True,False and/or opt@MeanFormat = "%5.2f"
;           2017/09/27         ADD function "unit_convert"
;           2017/10/18         ADD function "sort_for_xyplot"
;           2018/02/20         ADD procedure "cnColors_PosNeg"
;           2018/02/25         DEBUG procedure "cnColors_PosNeg", when num_positive and num_negative equals one, fspan will fail.
;                              ADD return opt in function "cnColors_PosNeg"
;           2018/03/19         MODIFY "ignorant_get_var_files" to get coordinate variables
;           2018/05/03         ADD "ignorant_get_vars_file"
;           2018/05/04         ADD "ignorant_convert_2Ddata"
;           2018/05/10         ADD opt_2Ddata@return_type = "lat_time-lon_avg" into "ignorant_convert_2Ddata"
;           2018/05/11         MODIFY "ignorant_convert_2Ddata"
;           2018/05/17         ADD "pressure_lat-timelon_avg","time_lat_lon-plev" into ignorant_convert_2Ddata 
;           2018/05/19         ADD ms2knots in "unit_convert"
;                              MODIFY "ignorant_convert_2Ddata"
;           2018/05/30         ADD W/m22mm/day in "unit_convert"
;           2018/05/31         ADD mm/day2m/s in "unit_convert"
;                              MODIFY "unit_convert", not do if input units differs the unit attribute
;           2018/06/07         ADD function gsn_csm_contour_map_siglvl
;           2018/06/08         ADD procedure change_var_atts
;           2018/06/14         ADD set cn fill colors in the procedure "change_var_atts"
;           2018/06/16         MODIFY longname setup in the case "cnColors_PosNeg"
;           2018/07/11         MODIFY ignorant_convert_2Ddata to make lon conversion automatically
;           2018/07/15         MODIFY reture_type name of ignorant_convert_2Ddata
;           2018/07/26         DEBUG in "cnColors_PosNeg", moving CopyAtt otherwise the previous cmnmap will override the correct one
;                              CHANGE dry_wet colors in "cnColors_PosNeg"
;           2018/08/02         DEBUG gsn_csm_contour_map_siglvl, user can decide whether use cnColors_PosNeg 
;           2018/08/04         ADD opt_chgatt@chg_cncolors in "change_var_atts"
;           2018/08/04         ADD a new procedure "ignorant_convert_1Ddata"
;           2018/08/05         MODIFY "ignorant_convert_1Ddata"
;           2018/08/06         ADD "tlji2pn_crossji" into ignorant_convert_2Ddata
;           2018/08/06         MIDIFY "change_var_atts"
;           2018/08/09         MODIFY option "tlji2pn_crossji" in the ignorant_convert_2Ddata
;           2018/08/11         DEBUG cnColors_PosNeg, it doesn't read input colormaps 
;           2018/08/14         DEBUG "tlji2pn_crossji"
;           2018/08/24         DEBUG several functions/procedures
;           2018/08/25         ADD a function, "ignorant_diagnose_var"
;           2018/08/26         DEBUG cnColors_PosNeg for input colormaps
;           2018/08/27         ADD "PW" in ignorant_diagnose_var
;                              MOVE change_var_atts forward 
;           2018/09/10         ADD new varialbes in ignorant_diagnose_var
;           2018/09/13         ADD a function, "dimsizes_compare"
;           2018/09/19         ADD opt@mapview in gsn_csm_contour_map_avg & gsn_csm_contour_map_siglvl
;           2018/09/21         ADD ignorant_diagnose_var into ignorant_convert_2Ddata
;                              ADD ### end of *func ### at the end of each of them
;           2018/10/01         ADD new colormap in cnColors and change_att
;                              ADD mean in gsn_csm_contour_map_siglvl
;                              ADD Meridional_energy_transport in ignorant_diagnose_var
;           2018/10/02         MODIFY get_vars_file & get_var_files to include diagnose_var
;           2018/10/05         ADD option "white_middle" in cnColors_PosNeg
;           2018/10/07         ADD "ftji2j_fj_avgti" & "vtji2j_vj_avgti" in ignorant_convert_1Ddata
;           2018/10/08         ADD new variable in ignorant_diagnose_var
;           2018/10/17         ADD "tji2ijOneD_avgt" in ignorant_convert_1Ddata
;           2018/10/19         MODIFY slightly
;           2018/10/21         DEBUG "tji2ijOneD_avgt" in ignorant_convert_1Ddata
;           2018/10/23         DEBUG cnColor_PosNeg
;           2018/10/23         ADD varname conversion in ignorant_diag_var
;           2018/10/24         ADD "ignorant_convert_varname"
;                              MODIFY varname conversion in ignorant_diag_var
;                              ADD x_coordinate "MPACE_3hr" in ignorant_convert_1Ddata
;           2018/10/25         Slightly modify
;           2018/10/27         Slightly modify
;           2018/10/28         add @cmap in cnColors_PosNeg
;                              Slightly modify
;           2018/10/29         Slightly modify
;           2018/11/14         ADD ignorant_ttest
;           2018/11/14         ADD opt@marker = "p_values" in ignorant_xy
;           2018/11/18         ADD cases in ignorant_convert_1Ddata
;           2018/11/26         Slightly modify
;           2018/11/29         Slightly modify
;           2018/12/05         Modify cnColors_PosNeg to use weaker color when number of positive/negative is less than 3
;           2018/12/22         ADD "monthly_series_tlji2ypm_avgji" into ignorant_convert_2Ddata
;           2019/01/02         Debug cnColors_PosNeg if all levels are > or < 0.
;           2019/01/10         ADD procedure ignorant_sameLB
;           2019/01/26         Update
;           2019/01/31         Update
;           2019/03/03         ADD "qv_zint" into ignorant_diagnose_var
;           2019/03/06         Modify ignorant_convert_2Ddata, "monthly_series_tji2yjm_avgi" to allow monthly series is not 12x times
;                              Modify gsn_csm_contour with siglvl to allow var_any & var_base dimension can be different
;           2019/03/09         Update
;           2019/04/10         ADD procedures, "select_landfrac" & "trim_fig"
;           2019/04/12         ADD functions, "ensemble_compute" & "ignorant_xy_ensemble"
;           2019/04/28         ADD function "convert_files_vars"
;           2019/04/29         ADD function "ignorant_get_vars_files"
;           2019/05/21         ADD procedure ignorant_sameScale
;           2019/05/24         ADD function diff_group
;           2019/05/26         Modify ignorant_sameScale
;           2019/05/31         Debug lat_wgt_avg, weighting factor is WRONG
;           2019/06/09         ADD procedure check_filedims
;           2019/08/15         Update
;                              ADD function get_VarLatLon
;                              include get_VarLatLon in gsn_csm_contour_map_avg & gsn_csm_contour_map_siglvl
;           2019/09/18         Update ignorant_xy_ensemble, make fill color more transparent
;           2019/09/24         Update ignorant_diagnose_var, change name of variable needed interpolated
;           2019/10/02         Add function "regline_print"
;           2019/10/03         Modify gsn_csm_contour_map_siglvl to avoid significance levels plot out to the map
;           2019/10/09         ADD procedure "overlay_siglvl_dots"
;           2019/10/19         ADD function "string_change" and "create_emptyNC"
;*****************************************************

;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;  Content of ignorant NCL program
;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;
; load "/home/chishin/ncarg/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "/home/chishin/ncarg/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;
;procedure overlay_siglvl_dots(wks[1]:graphic, plot[1]:graphic, xx:numeric, yy:numeric, pp:numeric, res[1]:logical, opt_plot[1]:logical)
;  add dots that pass p-value on a XY plot
;    wks : work station
;    plot: a XY plot
;    xx  : X values, a 1-D or 2D array
;    yy  : Y values, a 1-D or 2D array
;    pp  : p-values, a 1-D or 2D array
;    res : NCL res for this plot, e.g gs*
;    opt_plot: option
;      opt_plot@siglvl = 0.5  ; significance level. p-value less than this will be plotted
;      opt_plot@Draw = True   ; draw and frame the plot
;
; function regline_print(x[*]:numeric, y[*]:numeric, string_input[1]:string, opt[1]:logical)
;   print results from NCL linear regression function, regline_stats
;     x: an 1-D array of x values
;     y: an 1-D array of y values
;     string_input: string that will be printed on screen
;     opt: unused
;
; procedure check_filedims(infilenames[*]:string) 
;   check whether input files have the same dimension names and sizes
;   infilenames: file names, e.g. (/"file1","file2/)
;
; function get_VarLatLon(data_input:numeric, opt:logical)
;   retrieve variable in a certain lat/lon regions
;     return: a variable in a certain lat/lon regions
;     data_input: input data with 2 to 6 dimensions. The last two dimension must be (lat,lon)
;     opt       : option
;       opt@mpMaxLatF = 10.   ; upper boundary of lat
;       opt@mpMinLatF = 30.   ; lower boundary of lat. If max<min, max and min will exchange with each other
;       opt@mpMaxLonF = -30.  ; upper boundary of lon
;       opt@mpMinLonF = 50.   ; lower boundary of lon. If max<min, max and min will exchange with each other  
;                                                      If either LonF<0, LonFlip will be applied (change lon to -180~180)
; 
; function convert_files_vars(filenames[*]:string, varnames[*]:string, opt[1]:logical)
;   use ignorant_convert_1D(2D)data to create a variable
;   filenames: netCDF file names, e.g. (/"file1","file2/)
;   varnames : variable names, e.g. (/"T","Q"/)
;   opt      ; logical
;     opt@return_type = choice from 1D or 2D data, e.g. "ftji2t_ft_avgji"
;
; procedure ignorant_sameScale(var_work:numeric, nclres[1]:logical, opt_sameLB[1]:logical)
;   get cnLevels from var_work
;     var_work   : variable name
;     nclres     : NCL res
;     opt_sameLB : option
;       opt_sameLB@OptionScale   = "XScale", "YScale", "LBScale" (default), 
;                                  modify (trX), (trY), and (cnLevels) respectively
;       opt_sameLB@percentile    = 0,1,10  ; percentile, default is 5%
;       opt_sameLB@factor        = 0.9 (default). ; set factor to adjust min/max value
;       opt_sameLB@max_levels    = 15  (default). ; set contour levels
;       opt_sameLW@cnColors_PosNeg = True         ; use cnColors_PosNeg
;
; procedure ignorant_sameLB(infilenames[*]:string, varvar[1]:string, nclres[1]:logical, opt_sameLB[1]:logical)
;   get cnLevels from input data
;     infilenames: input file names
;     varvar     : variable name
;     nclres     : NCL res
;     opt_sameLB : option
;       opt_sameLB@factor        = 0.9 (default). set factor to adjust min/max value
;       opt_sameLB@max_levels    = 15  (default). set contour levels
;       opt_sameLB@diff2base     = True.        . find cn values for difference with respect to filename_base
;       opt_sameLB@filename_base = ""           . name of base file. If filename_base is set, diff2base will be set to True
;
; function ignorant_ttest(x:numeric, y:numeric, opt[1]:logical)
;   do Student's t-test except the first dimension and return p-values, prob. 
;   if x=(time,lat,lon) and y=(time,lat,lon), prob = (lat,lon)
;   x   : a 1D, 2D, or 3D array
;   y   : dimension is as the same as x
;   opt : option
;     opt@iflag = True ; set iflag = True in NCL built-in function, ttest
;
; function ignorant_diagnose_var (infilename[1]:string, varname[1]:string, opt[1]:logical)
;   diagnose a variable from existing varialbes in a file. These diagnosed varialbes can be written out to the input file.
;     infilename: input file name
;     varname   : varialbe name. If existing, return directly
;     opt       : option
;       opt@writeout = False/True  ; write diagnosed variables to the input file
;
; procedure change_var_atts (varname[1]:string, varin, res[1]:logical, opt_chgatt[1]:logical)
;   change units and resources of the input variable
;     varname	: a string that decides which kind of changes will be applied
;     varin  	: input variable
;     res    	: a NCL res
;     opt_chgatt: option. Turn on/off this procedure
;       opt_chgatt@chg_cncolors = True/False ; turn on/off changing contour colors
;
; function gsn_csm_contour_map_siglvl(wks:graphic, varvar_any[*][*][*]:numeric, varvar_base[*][*][*]:numeric, res_var[1]:logical, res_sig[1]:logical, opt[1]:logical)
;   plot the difference between varvar_any and varvar_base over map and overlay the significance level from the t-test results
;     wks		: work station
;     varvar_any	: a three dimension variable (time,lat,lon)
;     varvar_base	: a three dimension variable (time,lat,lon)
;     res_var		: NCL res for plotting varvar_any minus varvar_base
;     res_sig		: NCL res for plotting significance level
;     opt 		: options
;	opt@siglvl = 0.01  		; set significance level (the default is 0.05)
;       opt@cnColors_PosNeg = False	; turn off cnColors_PosNeg
;       opt@mapview    = "polar_NH"; use polar stereographic map centered at the North pole
;                      = "polar_SH"; use polar stereographic map centered at the South pole
;
; procedure ignorant_convert_1Ddata (infilenames[*]:string, varvars[*]:string, res[1]:logical, opt_1Ddata[1]:logical)
;   convert input variable to 1-D for plotting a XY plot
;     Return	: var_xx = opt@var_xx
;		: var_yy = opt@var_yy
;		:   if opt@var_xx_number = True, return var_xx = 0.,1.,2.,...N
;		:   if opt@var_yy_number = True, return var_yy = 0.,1.,2.,...N
;
;     infilename: input file name (CESM for now)
;     varvar    : variables in infilename
;     opt_1Ddata: option
;       opt_1Ddata@return_type = "vtji2t_vt_avgji"	; input: CESM (variables,time,lat,lon), return: var_xx=(time), var_yy=(varible,time), averaged over lat and lon
;       opt_1Ddata@return_type = "vtlji2vp_p_avgtji"	; input: CESM (variables,lev/ilev,time,lat,lon), return: var_xx=(variable,pressure), var_yy=(pressure), averaged over time, lat, and lon
;       opt_1Ddata@return_type = "vtji2j_vj_avgti"      ; input: CESM (variables,time,lat,lon), return: var_xx=(lat), var_yy=(variable,lat), averaged over time and lon
;       opt_1Ddata@return_type = "ftlji2fp_p_avgtji"	; input: CESM (files,time,lev/ilev,lat,lon), return: var_xx=(files,pressure), var_yy=(pressure), averaged over time, lat, and lon
;       opt_1Ddata@return_type = "ftji2t_ft_avgji"	; input: CESM (files,time,lat,lon), return: var_xx=(time), var_yy=(files,time), averaged over lat and lon
;       opt_1Ddata@return_type = "ftj2j_fj_avgt"        ; input: CESM (files,time,lat), return: var_xx=(lat), var_yy=(files,lat), averaged over time
;       opt_1Ddata@return_type = "ftji2j_fj_avgti"      ; input: CESM (files,time,lat,lon), return: var_xx=(lat), var_yy=(files,lat), averaged over time and lon
;       opt_1Ddata@return_type = "tji2ijOneD_avgt"	; input: (time,lat,lon), return: var_xx=1D(lat,lon), yy=1D(lat,lon), averaged over time
;
;       opt_1Ddata@@plev = (/250., 300., 350., 400., 450., 500./)	; pressure levels for interpolation
;
;       opt_1Ddata@@update_res = False	; turn off updating res
;
;       opt_1Ddata@x_coordinate = "12months_number" 		; set x-coordinate to 1,2,3,...,12
;       opt_1Ddata@x_coordinate = "12months_1characters" 	; set x-coordiante to J,F,M,...,D
;       opt_1Ddata@x_coordinate = "12months_3characters" 	; set x-coordinate to Jan,Feb,Mar,...,Dec
;
; function ignorant_convert_2Ddata (infilename[1]:string, varvar[1]:string, opt_2Ddata:logical)
;   convert input dataset to 2-D dataset for plotting a contour
;     infilename: input file name (CESM for now)
;     varvar    : variable in infilename
;     opt_2Ddata: option
;
;      if applicable,    opt_2Ddata@time_step   = 1	; set return time step
;
;        opt_2Ddata@return_type = "tlji2pj_avgti"	; input: CESM var(time,lev/ilev,lat,lon). return var(pressure,lat) averaged over time & lon 
;
;        opt_2Ddata@return_type = "tlji2pt_avgji"	; input: CESM var(time,lev/ilev,lat,lon). return var(pressure,time) averaged over lat & lon 
;           opt_2Ddata@plev        = (/500.,880.1000./)	;   pressure levels (hPa) for interpolation
;
;        opt_2Ddata@return_type = "tlji2tji_intp"	; input: CESM var(time,lev/ilev,lat,lon). return var(time,lat,lon) at the given pressure level
;           opt_2Ddata@plev        = 300.		;   interpolate the var at this pressure level (hPa)
;
;        opt_2Ddata@return_type = "tji2jt_avgi"         ; input: var(time,lat,lon). return var(latitdue,time) averaged over lon
;
;        opt_2Ddata@return_type = "tlji2pn_crossji"           ; input: CESM var(time,lev/ilev,lat,lon). return var(time,pressure,npoints), vertical cross section between two points.
          ;opt_2Ddata@linint2_points_cyclic = False            ; For function linint2_points. This should be set to True only if you have global data, but your longitude values don't quite wrap all the way around the globe. For example, if your longitude values go from, say, -179.75 to 179.75, or 0.5 to 359.5, then you would set this to True.
          ;opt_2Ddata@linint2_points_cyclic = True            
          ;opt_2Ddata@plev     = (/250., 300., 350., 400., 450., 500./)       ; pressure levels for interpolation

          ;opt_2Ddata@leftlat  =  30.         ; lat/lon of the left point,  lat=(-90,90), lon=(-180,180)
          ;opt_2Ddata@leftlon  =  -70.
          ;opt_2Ddata@rightlat =  -30.        ; lat/lon of the right point, lat=(-90,90), lon=(-180,180)
          ;opt_2Ddata@rightlon =  12.

          ;opt_2Ddata@time_step = "avg"                ; = 0,N -> plot at a specific step. = "avg" plot the time-averaged 
          ;opt_2Ddata@npts      =  100                 ; number of points of resulting transection
          ;opt_2Ddata@plot_map  = "transec_map"        ; plot 2-d contour with the transect on a map
          ;opt_2Ddata@plot_map  = "transec_plev_cnmap" ; plot 2-d contour with transect and contour at a pressure level on a map
          ;opt_2Ddata@plev_crossji = 800.            ;   set the pressure level for plotting for "transec_plev_cnmap"
;
; function ignorant_convert_varname(varname_in[1]:string, opt[1]:logical)
;   convert variable name among different datasets
;   varname_in: input variable name
;   opt       : option
;     opt@convert_dataset = "in2out", e.g. "cesm2scm"
;
; procedure cnColors_PosNeg (varvar[*][*]:numeric, res[1]:logical, opt_cncolors[1]:logical)
;   assign colors to positive and negative values in NCL resources for plotting contour
;   return: res@cnFillColors
;     varvar       : 2-D array for plotting
;     res          : NCL resources
;     opt_cncolors : a logical varialbe to control availabe features in this procedure
;        opt_cncolors@max_levels = 15  ; level for using function nice_mnmxintvl
;      
;        opt_cncolors@cnlevels_predefine = False    ; True  - contour levels are given in res
;                                                   ; False - contour levels will be computed in this procedure
;      
;        opt_cncolors@colomap_positive = "WhiteYellowOrangeRed"  ; colormap for plotting positive colors
;        opt_cncolors@colomap_positive_index_start = 2           ; the least positive color
;        opt_cncolors@colomap_positive_index_end = 999           ; the most  positive color
;      
;        opt_cncolors@colomap_negative = "MPL_GnBu"              ; colormap for plotting negative colors
;        opt_cncolors@colomap_negative_index_start = 127         ; the most  negative color
;        opt_cncolors@colomap_negative_index_end = 2             ; the least negative color

; function binmap2(xdata[*]:numeric, ydata:numeric, nbin:integer, opt:logical)
;          sort ydata into xdata bin based on user defined min,max,interval values
;          this is an updated version of binmap2
;            xdata: an input 1-D data
;            ydata: an input 1D or 2D data, the rightmost dimension must be the same as xdata
;            nbin : number of bins for xdata
;            opt  : logical variable
;              opt@x_max = 100.  ; max value for bin x, if not set, max(xdata) will be used
;              opt@x_min =   0.  ; min value for bin x, if not set, min(xdata) will be used
;
;          return: a 3D array, [variable | 1] x [index | 5] x [bin | 4]
;                     dim variable: the same as the first dimension of ydata
;                     dim index   : Index 0: number in each bin; 
;                                   Index 1: Sum in each bin; 
;                                   Index 2: Average value in each bin; 
;                                   Index 3: standard deviation in each bin; 
;                                   Index 4: PDF of each bin
;                     dim bin     : number of bins
;
; function sort_for_xyplot(xdata[*]:numeric, ydata[*]:numeric, zdata[*]:numeric, nbin:integer)
;          sort input zdata into several bins, and then return a data that can use for plotting a XY plot
;            xdata: an input 1-D data that would be plot on X-axis
;            ydata: an input 1-D data that would be plot on Y-axis
;            zdata: an input 1-D data that would be sorted into several bins
;            nbin : number of bins for zdata
;
;          return: a 3-D data plotvars (2,nbin,dimension of X,Y)
;                   plotvars@index = "0: x_plot_var, 1: y_plot_var"
;                   plotvars@nbin_labels = labels
;                   plotvars@nbin_mid    = layers
;                   plotvars@nbin_edge   = levels
;
; function unit_convert(var_from:numeric, unit_from:string, unit_to:string)
;          convert the unit of an input variable, var_from
;          var_from : any numerical variable
;          unit_from: unit of var_from. Can set "" if @units is an attribute of var_from
;          unit_to  : new unit
;
; function taylor_diagram (wks:graphic ,RATIO[*][*]:numeric, CC[*][*]:numeric \
;                                     ,rOpts:logical)
;          wks  : work station
;          RATIO: ratio of the normalized standard deviation respect to the reference variable  
;          CC   : correlation coefficient
;          rOpts: NCL resources for taylor diagram
;
;          ref: http://www.ncl.ucar.edu/Applications/taylor.shtml
;
; function binmap(xx:numeric,yy:numeric,xx_min[1]:numeric,xx_max[1]:numeric,xx_int[1]:numeric)
;          sort yy into xx bin based on user defined min,max,interval values
;          Return: zz(3,binsize), Index 0: number in each bin; Index 1: Average value in each bin; Index 2: standard deviation in each bin
;          xx    : multidimension array
;          yy    : multidimension array
;          xx_min: minimum value of bin
;          xx_max: maximum value of bin
;          xx_int: bin interal
;
; function lat_wgt_avg(data_in:numeric)
;          compute latitude-weighted mean fields
;          Return  : longitude and latitude-weighted mean fields
;          data_in : numeric array, acceptable dimension orders are (lat,lon),(time,lat,lon),(time,?,lat,lon), (time,?,?,lat,lon).
;
; function plus(a:numeric, b:numeric)
;          compute a plus b, return variable's dimension,attributes,coordinated is the same as a.
;          Return  : a plus b
;          a       : numeric array
;          b       : numeric number 
;
; function divide(a:numeric, b:numeric)
;          compute a divided by b value, return variable's dimension,attributes,coordinated is the same as a.
;          Return  : a divide b
;          a       : numeric array
;          b       : numeric number 
;
; function multiply(a:numeric, b:numeric[1])
;          compute a multiply b value, return variable's dimension,attributes,coordinated is the same as a.
;          Return  : a multiply b
;          a       : numeric array
;          b       : numeric number 
;
; function diff(a:numeric, b:numeric)
;          compute a minus b value, return variable's dimension,attributes,coordinated is the same as a.
;          Return  : a minus b
;          a       : numeric array
;          b       : numeric array, dimension is must the same as a
;
; function diff_group(var_in:numeric, option_group[1]:logical)
;          compute a
;          Return  : separate var_in into several group and do difference with the base variable
;                    e.g. var_in = var[10], nfile_per_group=5, base=var[4 & 9]
;          var_in  : numeric array
;          option_group : option
;            option_group@nfile_per_group = 3 ; file# per group
;
; function dimsizes_compare (array1, array2)
;          compare the dimension size of two input arrays
;          Return  : True/False
;          array1  : an array
;          array2  : an array
;
; function gsn_csm_contour_map_avg(wks[1]:graphic, data_input[*][*]:numeric, \
;                                  nclres[1]:logical, opt[1]:logical)
;          using NCL funtion gsn_csm_contour_map and ADD 2D-data average at @gsnCenterString 
;          there are also user-defined resources, opt
;
;          Return    : graphic, contour at map plot
;          wks       : work station
;          data_input: input 2-D data array
;          nclres    : any resources at NCL
;          opt       : user-define options
;            opt@LatWeighted = False    ; do average considering latitude transformation effects
;                                       ; if turn on, data_input coordinate variable's MUST in correct format
;            opt@MeanFormat = "%7.3f"   ; set mean value format
;            opt@MeanPrefix = "mean: "  ; set prefix of mean value
;            opt@mapview    = "polar_NH"; use polar stereographic map centered at the North pole
;                           = "polar_SH"; use polar stereographic map centered at the South pole
;
; function ignorant_get_vars_files(infilenames[*]:string, varvars[*]:string, opt)
;          get all varibles from different files, if [varvars] has dimension as (time,lat,lon)
;          the return dimension could be (files,variable,time,lat,lon)
;
;          Return         : same type as [read_var_name], but add two extra dimensions [files,variable]
;          infilenames    : string, ONE path of the read-in file
;          read_var_name  : string, read-in variables with the SAME dimension
;          opt            : option for ignorant_diagnose_var
;
; function ignorant_get_var_files(infile_names:string, read_var_name[1]:string)
;          get same varible from different files, if [read_var_name] has dimension as (time,lat,lon)
;          the return dimension could be (file,time,lat,lon)
;
;          Return         : same type as [read_var_name], but add a extra dimension [file]
;          infile_names   : one-dimension string, full path of read-in files
;          read_var_name  : string, ONE varible name in infile_names
;
; function ignorant_get_vars_file(infilenames[0]:string, varvars[*]:string)
;          get all varibles from a file, if [varvars] has dimension as (time,lat,lon)
;          the return dimension could be (variable,time,lat,lon)
;
;          Return         : same type as [read_var_name], but add a extra dimension [variable]
;          infilename s   : string, ONE path of the read-in file
;          read_var_name  : string, read-in variables with the SAME dimension
;
; function ignorant_combine_colormaps( cmap1[*][3]:numeric, cmap2[*][3]:numeric, \
;                                      choice:integer  )
;          combine 2 RGB array to one RGB array.
;          Can use for combine different color array.
;           
;          Return      : numeric, combined color map array
;          cmap1,cmap2 : First RGB array, in dimension [*]{3}.
;          choice      : .eq. 1 -> combine total cmap1,cmap2
;                        else   -> don't combine cmap2 fore/background color
;
; function ignorant_xy(wks:graphic, xx:numeric, yy:numeric , \
;                      labels:string, nclres:logical, opt:logical)
;          ignorant defined XY-plot function. gsn_csm_xy is the basic.
;
;          Return : graphic, XY plot.
;          wks    : work station
;          xx     : x-value
;          yy     : y-value
;          labels : each line title
;          nclres : any resources in NCL
;          opt    : user resources, not define
;            opt@legend_mean = True or False   ; turn on/off mean values on legends
;            opt@MeanFormat = "%6.2f"          ; format of the mean values
;
;            plot significance on a XY plot.  DO NOT SUPPORT panel plots
;              opt@marker = "p_values"
;              opt@prob   = p-values, 1D or 2D arrays. SHOULD be the same as either xx or yy
;              opt@siglvl = 0.05 ; significance level. prob<siglvl will be plotted. Defalut is 0.05
;
; function ignorant_merge_colormaps(icmap)
;          merge NCL predifined color map.
;          Return : merged-colormap RGB value array.
;          icmap  : an array contain NCL colormaps name, ex:(/"Cat12", "gsltod", "BrownBlue12"/)
;
; function ignorant_select_colormap(icmap_name:string, ocmap_index[*]:numeric, \
;                                   option_foreback:logical)
;          select some color in NCL-named colormap, and write out to an array, which can used to
;          define a new colormap.
;          icmap_name      : NCL-named colormap, ex: "gui_default"
;          ocmap_index     : an index array of choosed color in icmap_name, ex: (/3,6,8/)
;                            if ocmap_index = 999 then return total colors in given colormap
;          option_foreback : option that add fore- and background colors in return variables
;                            .True.  : add fore- and background colors
;                            .False. : not add fore- and background colors
;
; procedure ignorant_write_nc(filename:string,outvarname:string,invar,opt:integer)
;           Write value to an existing netCDF file with existing varible
;           filename: file which write in, note should be a support file type of NCL
;           outvarname : varible name in *filename* file, must be existing
;           invar : value, dimension MUST be equal to *outvarname* in *filename*, exclude unlimited dimension
;           opt : option, INTEGER, -1   -- outfile have unlimited time, overwrite from max time
;                                  -2   -- Don't write in unlimited time varible
;                                  -3   -- *outvarname* not exist in *filename*
;                                          create new varible and new dimsensions  (lazy method, not recommended)
;                                N.ge.0 -- outfile have unlimited time, overwrite in N time
;
; function ignorant_LWP(qc[*][*][*],p[*][*][*])
;          Calculate liquid water path using NCL "prcwater_dp" function
;          qc : cloud/rain/ice water mixing ratio, dimension is nz*ny*nx, unit in kg/kg
;           p : pressure, dimension is same as qc, unit in Pa
;
;          Return: a two dimension array ny*nx, mean every grid column liquid water path
;
;
; function ignorant_HMLcloud(P[*][*][*]:numeric,Tk[*][*][*]:numeric,qv[*][*][*]:numeric,qall[*][*][*]:numeric)
;          Calculate every grid box High(300-550mb),Middle(550-800mb),Low(800-1000mb) cloud amount using maximum overlap assumption
;          algorithm is same as Xu and Krueger(1991), Mon. Wea. Rev. 119, 342-367
;             P: pressure, dimension is nz*ny*nx , unit in hPa
;            Tk: temperature, dimension is same as P, unit in K
;            qv: water vapor mixing ratio, dimension is same as P, no unit (ex: kg/kg)
;          qall: condenced water mixing ratio, dimension is same as P, no unit (ex: kg/kg)
;
;          Return: High,Middle,Low level cloud amount, 3*ny*nx dimension array, index 0,1,2 corresponding to L,M,H cloud
;
;
; function ignorant_contour(wks[1]:graphic,data[*][*]:numeric,opt[1]:logical)
;          Default Behavier: Filled Color Contour with Label Bar
;
;
; function ignorant_nc_contour( wks[1]:graphic, nc_file:file, data[*][*]:numeric, opt:logical, nclres:logical )
;          Default Behavier: Draw Filled Color Contour on Map
;
;          opt@WRFfile = True        ; Set True if nc_file is WRF output
;
;          opt@REF_LAT = 20.023      ; Use in Map Projection
;          opt@REF_LON = 117.2858    ; KNOWNI,KNOWNJ correspond lon,lat value
;          opt@KNOWNI = 0            ;
;          opt@KNOWNJ = 0            ;
;
;          opt@XLeft   = XL          ; Set Map Area
;          opt@XRight  = XR          ; default is full range
;          opt@YTop    = YT          ;
;          opt@YBottom = YB          ;
;
;          nclres : Any resources in NCL
; 
; procedure select_landfrac(var, opt)
;           Using a land mask from a file, the user can set missing values to the part of input arrasy that does not fulfill the criteria given by the user.
;
;           var: a 3-D or 4-D array. The rightmost dimension must be (lat,lon) and consistent with the land mask file.
;           opt@resolution  : set land mask file. e.g. opt@resolution = "cesm_1.9x2.5"
;           opt@landfrac    : criteria, "land_only" (landfrac=1), "ocean_only" (landfrac=0), "land_0.5" (0.5<=landfrac<=1.)
; 
; procedure trim_fig(wks_type, wks_name, opt1)
;   trim the white space of a figure using "convert" command
;   
;   wks_type : the type of the figure, e.g. png
;   wks_name : the name of the figure
;     ==> The name of figure: wks_name+"."+wks_type
;
;   opt1@fig_name : directly give the the name of figure and omit wks_type & wks_name 
;
; function ensemble_compute(ens[*][*]:numeric, opt:logical)
;   compute ensemble mean, max, and min of the leftmost dimension
;
;     ens   : a two-dimension array (member, npoints). The leftmost dimension will be used to compute ensemble.
;     Return: a three-dimension array (group, value, npoints)
;
;     opt: options
;       opt@return_type = "mxmnag"  ; return max/min/avg values of the ensemble. Otherwise, return (group,member+mean,npoints)
;       opt@option_group = True     ; if grouping the leftmost dimension
;       opt@nmember_per_group = 3   ; number of member per group
;
; function ignorant_xy_ensemble (wks:graphic, xx:numeric, yy:numeric , \
;                                labels:string, nclres:logical, opt:logical)
;          plot ensemble mean & range using ignorant_xy
;
;          Return : graphic, XY plot.
;          wks    : work station
;          xx     : x-value. Either xx and yy MUST be a two-dimensional array.
;          yy     : y-value
;          labels : no effect here
;          nclres : any resources in NCL
;            nclres@gsnFrame = True    ; plot each figure
;            nclres@xyLineColors =     ; set line colors
;            nclres@gsnXYFillColors =  ; set fill colors
;          opt    : user resources, not define
;            opt_ens@ensemble = "x","y"   ; set which array is the ensemble. If not, this function will set the 2-D array to ensemble array
;            opt@return_type = "mxmnag"  ; (return group,max/min/avg,npoints)
;            opt@option_group = True     ; if grouping the first dimension
;            opt@nmember_per_group = 3   ; number of member per group
;            opt@XYStyle = "fill", "solid", "dash" (default is solid)
;
;!!!!!!!!!!!!!!!!!!!!!!
;  Program Code
;!!!!!!!!!!!!!!!!!!!!!!

;#############################################
;undef("")
;function ()
;begin
;
;end
;### end of ###
;#############################################

;#############################################
undef("string_change")
function string_change(string_in:string, opt:logical)
begin

  ; get changes of string_in 
  if (isatt(opt, "string_change")) then
    change_options = str_split(opt@string_change,"_")
  else
    change_options = ""
  end if

  nn = dimsizes(change_options)

  string_add = ""
  do i=0,nn-1
    change1 = change_options(i)

    ;*** put Delta before string_in ***
    if (change1.eq."delta") then
      string_add = "~F33~D~F21~"+string_add
    end if

    ;*** bold font of string_in ***
    if (change1.eq."bold") then
      string_add = "~F22~"+string_add
    end if
  end do
 
  ; return 
  string_out = string_add+string_in
  return(string_out)

end
;### end of string_change ###
;#############################################

;#############################################
undef("create_emptyNC")
procedure create_emptyNC(outfilename[*]:string)
begin

  file_empty = "~/script/ncl/empty.nc"
  file_emptytime = "~/script/ncl/emptytime.nc"
  file_emptyTime = "~/script/ncl/emptyTime.nc"

  file_empty1 = file_empty

  wrkdir=systemfunc("pwd")
  outfile1 = wrkdir+"/"+outfilename
  
  system("cp -i "+file_empty+" "+outfile1)
  system("echo create ["+outfile1+"]")

end
;### end of ###
;#############################################

;#############################################
undef("string_1Dto0D")
function string_1Dto0D(string_in[*]:string)
begin
  nn = dimsizes(string_in)
  string_out = ""
  do i=0,nn-1
    string_out = string_out + string_in(i) + ", "
  end do
  return(string_out)
end
;### end of string_1Dto0D ###
;#############################################

;#############################################
undef("get_VarLatLon")
function get_VarLatLon(data_input:numeric, opt:logical)
begin

  ;********************************
  ; get lat/lon from data_input
  ;********************************
  ;  lat = data_input&lat 
  ;  lon = data_input&lon 
  if (iscoord(data_input, "lat")) then
    lat = data_input&lat 
  else if (iscoord(data_input, "latitude")) then
    lat = data_input&latitude
  end if
  end if

  if (iscoord(data_input, "lon")) then
    lon = data_input&lon 
  else if (iscoord(data_input, "longitude")) then
    lon = data_input&longitude
  end if
  end if

  ;***********************
  ; get lat boundaries
  ;***********************
  ;*** find upper lat ***
  if (isatt(opt,"mpMaxLatF")) then
    lat_max = opt@mpMaxLatF
  else
    lat_max = max(lat)
  end if

  ;*** find lower lat ***
  if (isatt(opt,"mpMinLatF")) then
    lat_min = opt@mpMinLatF
  else
    lat_min = min(lat)
  end if

  if (lat_max.lt.lat_min) then
    tt = lat_max
    lat_max = lat_min
    lat_min = tt
  end if

  ;***********************
  ; get lon boundaries
  ;***********************

  ;*** if the lon is across 0., then flip the lon ***
  if (isatt(opt,"mpMaxLonF") .and. isatt(opt,"mpMinLonF")) then
    if (opt@mpMaxLonF.lt.0. .or. opt@mpMinLonF.lt.0) then
      data_input = lonFlip(data_input)
    end if
  end if

  ;*** find upper lon ***
  if (isatt(opt,"mpMaxLonF")) then
    lon_max = opt@mpMaxLonF
  else
    lon_max = max(lon)
  end if

  ;*** find lower lon ***
  if (isatt(opt,"mpMinLonF")) then
    lon_min = opt@mpMinLonF
  else
    lon_min = min(lon)
  end if

  if (lon_max.lt.lon_min) then
    tt = lon_max
    lon_max = lon_min
    lon_min = tt
  end if

  ;***********************
  ; return data_output
  ;***********************
  if      (dimsizes(dimsizes(data_input)).eq.2) then
    data_output = data_input({lat_min:lat_max},{lon_min:lon_max})
  else if (dimsizes(dimsizes(data_input)).eq.3) then
    data_output = data_input(:,{lat_min:lat_max},{lon_min:lon_max})
  else if (dimsizes(dimsizes(data_input)).eq.4) then
    data_output = data_input(:,:,{lat_min:lat_max},{lon_min:lon_max})
  else if (dimsizes(dimsizes(data_input)).eq.5) then
    data_output = data_input(:,:,:,{lat_min:lat_max},{lon_min:lon_max})
  else if (dimsizes(dimsizes(data_input)).eq.6) then
    data_output = data_input(:,:,:,:,{lat_min:lat_max},{lon_min:lon_max})
  end if
  end if
  end if
  end if
  end if

  return(data_output)

  ;printVarSummary(data_input)
  ;printVarSummary(data_output)
  ;print("lat min/max: "+lat_min+","+lat_max+"; lon min/max: "+lon_min+", "+lon_max)


end
;### end of get_VarLatLon ###
;#############################################

;#############################################
undef("check_filedims")
procedure check_filedims(infilenames[*]:string)
begin

  nfile = dimsizes(infilenames)

  ftemp_name_base = infilenames(nfile-1)
  ftemp1_base = addfile(ftemp_name_base,"r")
  fdims_base = getfiledimsizes(ftemp1_base)
  fdims_name_base = getvardims(ftemp1_base)
  ndim_base = dimsizes(fdims_base)
  ;print(fdims_base)  ; dimension size
  ;print(fdims_name_base) ; dimension name
  ;printVarSummary(ftemp1_base)

  do ff=0,nfile-2
    ftemp_name_any = infilenames(ff)
    ftemp1_any = addfile(ftemp_name_any,"r")
    fdims_any = getfiledimsizes(ftemp1_any)
    fdims_name_any = getvardims(ftemp1_any)

    if (all(fdims_base.eq.fdims_any) .and. all(fdims_name_base.eq.fdims_name_any)) then
      aa=0

    else

      do dd=0,ndim_base-1
        name1 = fdims_name_base(dd)
        ;name1 = "111"
        idx_any = ind(fdims_name_any.eq.name1)
        idx_base = ind(fdims_name_base.eq.name1)

        if (ismissing(idx_any)) then
          print("ERROR: dimension ["+name1+"] is not in ["+ftemp_name_any+"]")
          exit
        else
          dim1_base = fdims_base(idx_base)
          dim1_any  = fdims_any(idx_any)
          ;print(dim1_base+" , "+dim1_any)
          if (dim1_base.ne.dim1_any) then
            print("ERROR: dimension ["+name1+"] is not the same")
            ;print(ftemp_name_base+": "+name1+", "+dim1_base)
            print("  "+name1+" = "+dim1_base+", file: "+ftemp_name_base)
            print("  "+name1+" = "+dim1_any+", file: "+ftemp_name_any)
            exit
          end if
        end if
        ;print(idx_any)

      end do ; end do of dimensions
    end if ; end if of dimensions are not the same

    delete(ftemp_name_any)
    delete(ftemp1_any)
    delete(fdims_any)
    delete(fdims_name_any)
  end do ; end do of files

end
;### end of check_filedims ###
;#############################################


;#############################################
undef("ensemble_compute")
function ensemble_compute(ens[*][*]:numeric, opt:logical)
begin
  ;opt = True
  ;opt@return_type = "mxmnag"  ; (return group,max/min/avg,npoints)
  ;opt@option_group = True     ; if grouping the first dimension
  ;opt@nmember_per_group = 3   ; number of member per group

  ;----------------
  ; get dimensions
  ;----------------
  ntotal  = dimsizes(ens(:,0))
  npoints = dimsizes(ens(0,:))

  ;----------------
  ; if grouping
  ;----------------
  if (       isatt(opt,"option_group") \
       .and. isatt(opt,"nmember_per_group")   .and. isint(opt@nmember_per_group) \
       .and. opt@nmember_per_group.ne.ntotal  .and. mod(ntotal,opt@nmember_per_group).eq.0 ) then
    nmember = opt@nmember_per_group
    ngroup  = ntotal/nmember

    ;--------------------------------------------------------
    ; create the return array as (group, mx/mn/ag, npoints)
    ;--------------------------------------------------------
    ens_mxmnag_temp = reshape(ens, (/ngroup,nmember,npoints/))
    ;print(ens_mxmnag)

    ens_mxmnag_avg = dim_avg_n_Wrap(ens_mxmnag_temp,1)
    ens_mxmnag_max = dim_max_n_Wrap(ens_mxmnag_temp,1)
    ens_mxmnag_min = dim_min_n_Wrap(ens_mxmnag_temp,1)
    ;print("max= "+ens_mxmnag_max+", min= "+ens_mxmnag_min+", avg= "+ens_mxmnag_avg)

    ens_mxmnag_temp0 = (/ens_mxmnag_max,ens_mxmnag_min,ens_mxmnag_avg/)
    ens_mxmnag_temp0!1 = "group"
    ens_mxmnag_temp0!0 = "value"
    ens_mxmnag_temp0!2 = "points"

    ens_mxmnag = ens_mxmnag_temp0(group|:, value|:, points|:)
    ens_mxmnag@value_index = "0:ensemble max, 1:ensembel min, 2:ensemble mean"

    ;----------------------------------------------------------
    ; create the return array as (group, member+mean, npoints)
    ;----------------------------------------------------------
    ens_mem = new( (/ngroup,nmember+1,npoints/), typeof(ens) )

    ens_temp = reshape(ens, (/ngroup,nmember,npoints/))
    ens_mem(:,0:nmember-1,:) = (/ens_temp/)
    ens_mem(:,nmember,:)     = dim_avg_n_Wrap(ens_temp,1)

    ens_mem!0 = "group"
    ens_mem!1 = "value"
    ens_mem!2 = "points"
    ens_mem@value_index = "0 to n-1:ensemble value, n:ensemble mean"
  
  ;----------------
  ; not grouping 
  ;----------------
  else
    print("WARNING (ensemble_compute): grouping fails. Compute ensembel mean using all input")

    ;--------------------------------------------------------
    ; create the return array as (group, mx/mn/ag, npoints)
    ;--------------------------------------------------------
    ens_avg = dim_avg_n_Wrap(ens,0)
    ens_max = dim_max_n_Wrap(ens,0)
    ens_min = dim_min_n_Wrap(ens,0)
    ;print("max= "+ens_max+", min= "+ens_min+", avg= "+ens_avg)

    ens_mxmnag = new((/1,3,npoints/),typeof(ens))

    ens_mxmnag(0,0,:) = (/ens_max/)
    ens_mxmnag(0,1,:) = (/ens_min/)
    ens_mxmnag(0,2,:) = (/ens_avg/)

    ens_mxmnag!0 = "group"
    ens_mxmnag!1 = "value"
    ens_mxmnag@value_index = "0:ensemble max, 1:ensembel min, 2:ensemble mean"

    ;----------------------------------------------------------
    ; create the return array as (group, member+mean, npoints)
    ;----------------------------------------------------------
    nmember = dimsizes(ens(:,0))

    ens_mem = new((/1,nmember+1,npoints/),typeof(ens))

    ens_mem(0,0:nmember-1,:) = ens(:,:)
    ens_mem(0,nmember,:)     = dim_avg_n_Wrap(ens,0)

    ens_mem!0 = "group"
    ens_mem!1 = "value"
    ens_mem!2 = "points"
    ens_mem@value_index = "0 to n-1:ensemble value, n:ensemble mean"

  end if

  ;----------------
  ; return
  ;----------------
  if (isatt(opt,"return_type") .and. opt@return_type.eq."mxmnag") then
    ens_return = ens_mxmnag 
  else
    ens_return = ens_mem
  end if

  copy_VarAtts(ens, ens_return)
  return(ens_return)
  ;printVarSummary(ens)
  ;printVarSummary(ens_return)
  ;print(ens_return)

end
;### end of ensemble_compute ###
;#############################################


;#############################################
undef("select_landfrac")
procedure select_landfrac(varvar_tmp, opt[1]:logical)
begin

;-------------------------------
; determine land fraction file 
;-------------------------------
  if (isatt(opt,"resolution") .and. opt@resolution .eq. "cesm_1.9x2.5_global") then
    ;filename_land = "/glade/u/home/yihsuan/data/landfrac/data_landfrac-cesm_1.9x2.5_global.nc"
    filename_land = "~/script/data/landfrac/data_landfrac-cesm_1.9x2.5_global.nc"
    varname_land  = "LANDFRAC"
  end if

  if (isatt(opt,"resolution") .and. opt@resolution .eq. "cesm_1.9x2.5_Arctic") then
    filename_land = "~/script/data/landfrac/data_landfrac-cesm_1.9x2.5_Arctic.nc"
    varname_land  = "LANDFRAC"
  end if

  if (isatt(opt,"resolution") .and. opt@resolution .eq. "cesm_1.9x2.5_Antarctic") then
    filename_land = "~/script/data/landfrac/data_landfrac-cesm_1.9x2.5_Antarctic.nc"
    varname_land  = "LANDFRAC"
  end if

  file_land = addfile(filename_land, "r")
  landfrac = file_land->$varname_land$(0,:,:)

  option_landfrac = True

;--------------------------------------------------
; determine land fraction upper and lower boundary
;   keep lower <= any points <= upper
;--------------------------------------------------
  if (isatt(opt,"landfrac") .and. opt@landfrac .eq. "land_only") then
    landfrac_upper = 1. ; only keep land=1 points
    landfrac_lower = 1.

  else if (isatt(opt,"landfrac") .and. opt@landfrac .eq. "non_land") then
    landfrac_upper =  0. ; only keep land=0 points, i.e. ocean & lakes
    landfrac_lower =  0.

  else if (isatt(opt,"landfrac") .and. opt@landfrac .eq. "land_0.5") then
    landfrac_upper = 1.  ; keep 0.5<landfrac<1 points
    landfrac_lower = 0.5
  else
    option_landfrac = False
  end if
  end if
  end if

;------------------
; change values
;------------------
  dim_vartmp = dimsizes(dimsizes(varvar_tmp))

  if (option_landfrac .and. dim_vartmp.ge.2) then
    ; get index
    landfrac_1D      = ndtooned(landfrac)
    dsizes_landfrac  = dimsizes(landfrac)
    landfrac_indices = ind_resolve(ind(landfrac_1D.le.landfrac_upper.and.landfrac_1D.ge.landfrac_lower),dsizes_landfrac)

    varvar_tmp0 = varvar_tmp
    varvar_tmp = -999.99
    varvar_tmp@_FillValue = -999.99
  
    npts = dimsizes(landfrac_indices(:,0))
    do i=0,npts-1
      lat1 = landfrac_indices(i,0)
      lon1 = landfrac_indices(i,1)
      ;print(i+" , "+lat1+" , "+lon1+" , "+ varvar_tmp(0,lat1,lon1))
      if (dim_vartmp.eq.3) then
        varvar_tmp(:,lat1,lon1) = varvar_tmp0(:,lat1,lon1)
      else if (dim_vartmp.eq.2) then
        varvar_tmp(lat1,lon1) = varvar_tmp0(lat1,lon1)
      else if (dim_vartmp.eq.4) then
        varvar_tmp(:,:,lat1,lon1) = varvar_tmp0(:,:,lat1,lon1)
      end if
      end if
      end if
    end do

    varvar_tmp@land_selection = opt@landfrac
    print("func select_landfrac: "+opt@landfrac)
  end if ; end if of option_landfrac

end  ; end procedure of select_landfrac
;#############################################

;#############################################
undef("ignorant_ttest")
function ignorant_ttest(x:numeric, y:numeric, opt[1]:logical)
begin

;---------------
; initialzation
;---------------
  ; iflag
  ; Set to False if the two original samples are assumed to have the same population variance. Set to True if the two original samples are assumed to have different population variances. The latter applies the Welsh degree-of-freedom modification (Welsh's t-test).

  ; tval_opt
  ;Set to True if the Student t-values are to be returned in addition to the statistical probabilities. Set to False if only the probabilities are desired.

  if (isatt(opt,"iflag")) then
    iflag = opt@iflag
  else
    iflag = False
  end if

  ; tval_opt only support False at present
  tval_opt = False

  ;if (isatt(opt,"tval_opt")) then
  ;  tval_opt = opt@tval_opt
  ;else
  ;  tval_opt = False
  ; end if

;-------------
; do t-test
;-------------
  nn = dimsizes(dimsizes(x))
  prob = dim_avg_n_Wrap(x,0)
  prob = 0.
  prob@long_name = "p-value from Student's t-test"
  prob@iflag     = iflag
  prob@tval_opt  = tval_opt

  if (nn.eq.1) then
    xtmp = x
    ytmp = y

    ;*** compute time-average and variance ***
    xAve = dim_avg_Wrap (xtmp)              ; calculate means at each grid point
    yAve = dim_avg_Wrap (ytmp)
    xVar = dim_variance_Wrap (xtmp)         ; calculate variances
    yVar = dim_variance_Wrap (ytmp)

    ;*** get size ***
    xN = dimsizes(xtmp)
    yN = dimsizes(ytmp)

    prob = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, tval_opt)
  end if

  if (nn.eq.2) then
    qq = dimsizes(prob)
    ii = qq(0)

    do i=0,ii-1
      xtmp = x(:,i)
      ytmp = y(:,i)

      ;*** compute time-average and variance ***
      xAve = dim_avg_Wrap (xtmp)              ; calculate means at each grid point
      yAve = dim_avg_Wrap (ytmp)
      xVar = dim_variance_Wrap (xtmp)         ; calculate variances
      yVar = dim_variance_Wrap (ytmp)

      ;*** get size ***
      xN = dimsizes(xtmp)
      yN = dimsizes(ytmp)

      prob(i) = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, tval_opt)
    end do
  end if

  if (nn.eq.3) then
    qq = dimsizes(prob)
    ii = qq(0)
    jj = qq(1)

    do i=0,ii-1
    do j=0,jj-1

      xtmp = x(:,i,j)
      ytmp = y(:,i,j)

      ;*** compute time-average and variance ***
      xAve = dim_avg_Wrap (xtmp)              ; calculate means at each grid point
      yAve = dim_avg_Wrap (ytmp)
      xVar = dim_variance_Wrap (xtmp)         ; calculate variances
      yVar = dim_variance_Wrap (ytmp)

      ;*** get size ***
      xN = dimsizes(xtmp)
      yN = dimsizes(ytmp)

      prob(i,j) = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, tval_opt)
    end do
    end do
  end if

;---------
; return
;---------
  return(prob)

end
;### end of ignorant_ttest ###
;#############################################

;#############################################
undef("ignorant_convert_varname")
function ignorant_convert_varname(varname_in[1]:string, opt[1]:logical)
begin

  ;opt@convert_dataset = "scm2cesm"

;----------------
; program start
;----------------
  varname_dataset = (/"cesm","scm"/)	; name of dataset
  ;var_strings = (/cesm, scm/)
  varname_strings = (/  \
                  (/"CLDLOW"  , "lowcld"/)    , \
                  (/"CLDMED"  , "midcld"/)    , \
                  (/"CLDHGH"  , "hghcld"/)    , \
                  (/"CLDTOT"  , "totcld"/)    , \
                  (/"TGCLDLWP", "cldliq"/)    , \
                  (/"EVAP"    , "evppsrf"/)   , \
                  (/"FLUT"    , "TOA_LWup"/)  , \
                  (/"FLUS"    , "srflwup"/)   , \
                  (/"FLDS"    , "srflwdn"/)   , \
                  (/"FSUS"    , "srfswup"/)   , \
                  (/"FSDS"    , "srfswdn"/)   , \
                  (/"FSNTOA"  , "TOA_SWdn"/)  , \
                  (/"LHFLX"   , "lhflx" /)    , \
                  (/"OMEGA"   , "omega"/)     , \
                  (/"PRECT"   , "Prec"/)      , \
                  (/"Q"       , "q"/)         , \
                  (/"TMQ"     , "prew"/)      , \
                  (/"RELHUM"  , "REHUM"/)     , \
                  (/"SHFLX"   , "shflx" /)    , \
                  (/"SOLIN"   , "TOA_ins"/)   , \
                  (/"TREFHT"  , "Tsair"/)     , \
                  (/"TS"      , "Tskin"/)     , \
                  (/"U"       , "u"/)         , \
                  (/"V"       , "v"/)         , \
                  (/"U10"     , "windsrf"/)   , \
                  (/"zzldjss" , "skdjfksjdfsdf"/) \   ; random varnameiable, just set to the last varnameiable
                  ;(/""       , ""/)  , \
                /)

  nstring = dimsizes(varname_strings(:,0))
  ndata  = dimsizes(varname_dataset)

  ;----------------------------
  ; find input/output dataset
  ;----------------------------
  if ( isatt(opt, "convert_dataset")) then
    convert_dataset = opt@convert_dataset 
  else 
    print("ERROR: (ignorant_convert_varname) opt@convert_dataset is not given, e.g. opt@convert_dataset = scm2cesm")
    ;convert_dataset = "zzkkdd2nnssoo"
    convert_dataset = "cesm2cesm"
  end if

  string1 = str_split(convert_dataset,"2")
  ;print(string1)
  do i=0,ndata-1
    if (string1(0) .eq. varname_dataset(i)) then
      dataset_in = varname_dataset(i)
      index_in   = i
    end if
    if (string1(1) .eq. varname_dataset(i)) then
      dataset_out = varname_dataset(i)
      index_out   = i
    end if
  end do
  ;print(dataset_in+" , "+dataset_out)
 
  ;----------------------------
  ; find corresponding varname
  ;----------------------------
  varname_out = varname_in
  if (isdefined("index_in") .and. isdefined("index_out")) then
    do i=0,nstring-1
      if (varname_in .eq. varname_strings(i,index_in)) then
        varname_out = varname_strings(i,index_out)
      end if
    end do
    ;print(opt@convert_dataset+" , "+varname_in+" , "+varname_out)
  
    if (varname_out .ne. varname_in) then
      print("Done. (ignorant_convert_varname) convert ["+varname_in+"] to ["+varname_out+"] via ["+convert_dataset+"]")
  
    else
      print("WARNING: (ignorant_convert_varname) varname_in ["+varname_in+"] is not in ["+convert_dataset+"], return ["+varname_in+"]")
    end if

  else
    print("ERROR: (ignorant_convert_varname) invalid dataset ["+convert_dataset+"]")
    print("  Supported dataset")
    print("    "+varname_dataset)
    print("Return input variable name ["+varname_in+"]")
  end if

  ;print(varname_out)
  return(varname_out)

end
;### end of ignorant_convert_varname ###
;#############################################

;###################################################
undef("cnColors_PosNeg")
procedure cnColors_PosNeg (varvar:numeric, res[1]:logical, opt_cncolors[1]:logical)
begin

;-------------
; initialize
;-------------
  option_cmap = "default"

  if ( isatt(varvar,"long_name") ) then
    longname = varvar@long_name
  else
    longname = "slkdjflksjdfljasldfjlaks"
  end if

  ;if (\
  ;        longname.eq."" \
  ;   .or. longname.eq."" \
  ;   )
  ;  option_cmap = ""
  ;end if 

  ;if ( \
  ;         .not.ismissing(str_match(longname,"")) \
  ;    .or. .not.ismissing(str_match(longname,"")) \
  ;   ) then
  ;   option_cmap = ""
  ;end if

  if ( \
           .not.ismissing(str_match(longname,"cooling rate")) \
      .or. .not.ismissing(str_match(longname,"Cooling rate")) \
     ) then
     option_cmap = "warm_cold"
  end if

  if ( \
           .not.ismissing(str_match(longname,"Vertical velocity")) \
      .or. .not.ismissing(str_match(longname,"OMEGA")) \
     ) then
     option_cmap = "red_blue"
  end if

  if ( \
           .not.ismissing(str_match(longname,"Zonal wind")) \
      .or. .not.ismissing(str_match(longname,"zonal wind")) \
      .or. .not.ismissing(str_match(longname,"Meridional")) \
      .or. .not.ismissing(str_match(longname,"Meridional wind")) \
      .or. .not.ismissing(str_match(longname,"meridional wind")) \
      .or. .not.ismissing(str_match(longname,"U at")) \
      .or. .not.ismissing(str_match(longname,"V at")) \
      .or. .not.ismissing(str_match(longname,"V at")) \
      .or. .not.ismissing(str_match(longname,"meridional stream function")) \
      .or. .not.ismissing(str_match(longname,"Meridional stream function")) \
      .or. .not.ismissing(str_match(longname,"Meridional Stream Function")) \
     ) then
     option_cmap = "purple_green"
  end if
  
  if (\
           .not.ismissing(str_match(longname,"precipitation rate")) \
      .or. .not.ismissing(str_match(longname,"Precipitation rate")) \
      .or. .not.ismissing(str_match(longname,"Prec")) \
      .or. .not.ismissing(str_match(longname,"snow rate")) \
      .or. .not.ismissing(str_match(longname,"water path")) \
      .or. .not.ismissing(str_match(longname,"humidity")) \
      .or. .not.ismissing(str_match(longname,"precipitable water")) \
      .or. .not.ismissing(str_match(longname,"Precipitable water")) \
      .or. .not.ismissing(str_match(longname,"moisture flux convergence")) \
      .or. .not.ismissing(str_match(longname,"Moisture flux convergence")) \
      .or. .not.ismissing(str_match(longname,"qU")) \
      .or. .not.ismissing(str_match(longname,"sea ice")) \
      .or. .not.ismissing(str_match(longname,"sea-ice")) \
      .or. .not.ismissing(str_match(longname,"Sea ice")) \
      .or. .not.ismissing(str_match(longname,"Sea-ice")) \
     )
    option_cmap = "dry_wet"
  end if

  if (isatt(opt_cncolors,"cmap")) then
    option_cmap = opt_cncolors@cmap
  end if

;----------------
; program start
;----------------

  ;*** intialize opt_work ***
  opt_work = True
  opt_work@max_levels = 15  ; level for using function nice_mnmxintvl
  opt_work@cnlevels_predefine = False    ; True  - contour levels are given in res

  ;*** if both colormap_positive and colormap_nagative are defined ***
  if (isatt(opt_cncolors,"colormap_positive") .and. isatt(opt_cncolors,"colormap_negative")) then
    opt_work@colormap_positive_index_start = 2           ; the least positive color
    opt_work@colormap_positive_index_end = 999           ; the most  positive color

    opt_work@colormap_negative_index_start = 2           ; the least negative color
    opt_work@colormap_negative_index_end = 999           ; the most  negative color

  ;*** if only colormap_positive is defined ***
  else if (isatt(opt_cncolors,"colormap_positive") .and. .not.isatt(opt_cncolors,"colormap_negative")) then
    opt_work@colormap_negative = "MPL_GnBu"              ; colormap for plotting negative colors
    opt_work@colormap_negative_index_start = 127         ; the most  negative color
    opt_work@colormap_negative_index_end = 2             ; the least negative color

  ;*** if only colormap_negative is defined ***
  else if (.not.isatt(opt_cncolors,"colormap_positive") .and. isatt(opt_cncolors,"colormap_negative")) then
    opt_work@colormap_positive = "WhiteYellowOrangeRed"  ; colormap for plotting positive colors
    opt_work@colormap_positive_index_start = 2           ; the least positive color
    opt_work@colormap_positive_index_end = 999           ; the most  positive color

  ;*** if both colormap_positive and colormap_nagative are NOT defined ***
  else if (.not.isatt(opt_cncolors,"colormap_positive") .and. .not. isatt(opt_cncolors,"colormap_negative")) then

    ;opt_work@cnlevels_predefine = True      False - contour levels will be computed in this procedure
    opt_work@colormap_positive = "WhiteYellowOrangeRed"  ; colormap for plotting positive colors
    opt_work@colormap_positive_index_start = 2           ; the least positive color
    opt_work@colormap_positive_index_end = 999           ; the most  positive color
  
    opt_work@colormap_negative = "MPL_GnBu"              ; colormap for plotting negative colors
    opt_work@colormap_negative_index_start = 127         ; the most  negative color
    opt_work@colormap_negative_index_end = 2             ; the least negative color
  
    ;*** set individual colormap here ***
    if (option_cmap .eq. "dry_wet") then
      opt_work@colormap_positive = "MPL_YlGnBu"
      opt_work@colormap_positive_index_start = 10      ; the least positive color
      opt_work@colormap_positive_index_end = 113       ; the most  positive color
  
      opt_work@colormap_negative = "MPL_pink"
      opt_work@colormap_negative_index_start = 17      ; the most  negative color
      opt_work@colormap_negative_index_end = 97       ; the least negative color
      opt_work@white_middle = True
    end if
  
    if (option_cmap .eq. "purple_green") then
      opt_work@colormap_positive = "MPL_Purples"       ; colormap for plotting positive colors
      ;opt_work@colormap_positive_index_start = 2       ; the least positive color
      opt_work@colormap_positive_index_start = 49       ; the least positive color
      opt_work@colormap_positive_index_end = 999       ; the most  positive color
  
      opt_work@colormap_negative = "MPL_Greens"        ; colormap for plotting negative colors
      opt_work@colormap_negative_index_start = 120     ; the most  negative color
      ;opt_work@colormap_negative_index_end = 2         ; the least negative color
      opt_work@colormap_negative_index_end = 33         ; the least negative color
    end if
  
    if (option_cmap .eq. "red_blue") then
      opt_work@colormap_positive = "NCV_blu_red"         ; colormap for plotting positive colors
      opt_work@colormap_positive_index_start = 123       ; the least positive color
      opt_work@colormap_positive_index_end = 2           ; the most  positive color
  
      opt_work@colormap_negative = "NCV_blu_red"         ; colormap for plotting negative colors
      opt_work@colormap_negative_index_start = 255       ; the most  negative color
      opt_work@colormap_negative_index_end = 137         ; the least negative color
    end if

    if (option_cmap .eq. "blue_red") then
      opt_work@colormap_positive = "NCV_blu_red"         ; colormap for plotting positive colors
      opt_work@colormap_positive_index_start = 137       ; the least positive color
      opt_work@colormap_positive_index_end = 255           ; the most  positive color
  
      opt_work@colormap_negative = "NCV_blu_red"         ; colormap for plotting negative colors
      opt_work@colormap_negative_index_start = 17       ; the most  negative color
      opt_work@colormap_negative_index_end = 123         ; the least negative color
    end if

    if (option_cmap .eq. "warm_cold") then
      opt_work@colormap_positive = "MPL_GnBu"  ; colormap for plotting positive colors
      opt_work@colormap_positive_index_start = 2           ; the least positive color
      opt_work@colormap_positive_index_end = 127           ; the most  positive color

      opt_work@colormap_negative = "WhiteYellowOrangeRed"              ; colormap for plotting negative colors
      opt_work@colormap_negative_index_start = 250         ; the most  negative color
      opt_work@colormap_negative_index_end = 2             ; the least negative color
    end if

  end if 
  end if 
  end if 
  end if ; end if of colormap_positive and colormap_nagative

  copy_VarAtts(opt_cncolors, opt_work)
  if (opt_cncolors .and. isatt(res,"cnFillPalette")) then
    delete(res@cnFillPalette)
  end if

  ;*********************
  ; get contour levels
  ;*********************

  ; given res@cnLevelSelectionMode = "ExplicitLevels"
  if ( opt_work@cnlevels_predefine .and. isatt(res,"cnLevelSelectionMode") .and. res@cnLevelSelectionMode .eq. "ExplicitLevels" ) then
    if (isatt(res,"cnLevels")) then
      cnlevels = res@cnLevels
      option_cnlevels = True
    else
      option_cnlevels = False
    end if

  ; given res@cnLevelSelectionMode = "ManualLevels"
  else if ( opt_work@cnlevels_predefine .and. isatt(res,"cnLevelSelectionMode") .and. res@cnLevelSelectionMode .eq. "ManualLevels" ) then
    cnmin = res@cnMinLevelValF
    cnmax = res@cnMaxLevelValF
    cnint = res@cnLevelSpacingF

    if ( cnint .ne. 0. ) then
      n1 = toint((cnmax-cnmin)/cnint) + 1
      ;print(n1)
      cnlevels = fspan(cnmin,cnmax,n1)
      option_cnlevels = True
    else
      option_cnlevels = False
    end if

  ; compute cnlevels in this procedure
  else
    cnmin = min(varvar)    
    cnmax = max(varvar)    

    max_levels = opt_work@max_levels

    mnmxint = nice_mnmxintvl(min(varvar), max(varvar), max_levels, False)
    ;print(mnmxint)
    cnmin = mnmxint(0)
    cnmax = mnmxint(1)
    cnint = mnmxint(2)

    if ( cnint .ne. 0. ) then
      n1 = toint((cnmax-cnmin)/cnint) + 1
      ;print(n1)
      cnlevels = fspan(cnmin,cnmax,n1)
      option_cnlevels = True
    else
      option_cnlevels = False
    end if

  end if 
  end if 
  ;print(cnlevels)

if (option_cnlevels) then

  ;*****************************************************************
  ; find indexes of zero, smallest positive, and smallest negative
  ;*****************************************************************
  num_cnlevels = dimsizes(cnlevels)
  num_cncolors = num_cnlevels+1
  index_zero = closest_val(0., cnlevels)

  if ( cnlevels(index_zero) .eq. 0. ) then
    num_negative = index_zero+1
    num_positive = num_cncolors - num_negative

  else if ( cnlevels(index_zero) .gt. 0. ) then
    num_negative = index_zero+1
    num_positive = num_cncolors - num_negative

  else if ( cnlevels(index_zero) .lt. 0. ) then
    num_negative = index_zero+1
    num_positive = num_cncolors - num_negative
  end if
  end if
  end if

  ;*** if all cnlevels are > 0. ***
  if (all(cnlevels.ge.0.)) then 
    num_negative = 0
    num_positive = num_cncolors
  end if

  ;*** if all cnlevels are < 0. ***
  if (all(cnlevels.le.0.)) then 
    num_negative = num_cncolors
    num_positive = 0
  end if
 
  ;print("num_negative: "+num_negative+", num_positive: "+num_positive)

  ;************
  ; set colors
  ;************ 
  cncolors = new( (/num_cncolors,4/) , "float" )
  cncolors = 1  ; initialize colors to white 

  cmap_positive = read_colormap_file(opt_work@colormap_positive)
  cmap_negative = read_colormap_file(opt_work@colormap_negative)

  num_cmap_positive = dimsizes(cmap_positive(:,0))
  num_cmap_negative = dimsizes(cmap_negative(:,0))

  ; set positive colors
  if (isatt(opt_work, "colormap_positive_index_start") .and. opt_work@colormap_positive_index_start .lt. num_cmap_positive) then
    cmap_index_pos_start = opt_work@colormap_positive_index_start
  else
    cmap_index_pos_start = 2
  end if

  if (isatt(opt_work, "colormap_positive_index_end") .and. opt_work@colormap_positive_index_end .lt. num_cmap_positive) then
    cmap_index_pos_end = opt_work@colormap_positive_index_end
  else
    cmap_index_pos_end = num_cmap_positive-1
  end if

  ; set negative colors
  if (isatt(opt_work, "colormap_negative_index_start") .and. opt_work@colormap_negative_index_start .lt. num_cmap_negative) then
    cmap_index_neg_start = opt_work@colormap_negative_index_start
  else
    cmap_index_neg_start = 2
  end if

  if (isatt(opt_work, "colormap_negative_index_end") .and. opt_work@colormap_negative_index_end .lt. num_cmap_negative) then
    cmap_index_neg_end = opt_work@colormap_negative_index_end
  else
    cmap_index_neg_end = num_cmap_negative-1
  end if

  ; if number of negative values is too few, use 1/3 color index
  if (.not.isatt(opt_cncolors, "colormap_negative_index_start")) then
    if (num_negative .le. 3) then
      cmap_index_neg_start = (num_cmap_negative-1)/3
    end if
  end if

  ; if number of positive values is too few, use 1/3 color index
  if (.not.isatt(opt_cncolors, "colormap_positive_index_end")) then
    if (num_positive .le. 3) then
      cmap_index_pos_end = (num_cmap_positive-1)/3
    end if
  end if

  ; get index of colormap
  if (num_positive .gt. 1) then
    cmap_index_positive = toint (fspan(tofloat(cmap_index_pos_start),tofloat(cmap_index_pos_end),num_positive))
  else
    cmap_index_positive = cmap_index_pos_start
  end if

  if (num_negative .gt. 1) then
    cmap_index_negative = toint (fspan(tofloat(cmap_index_neg_start),tofloat(cmap_index_neg_end),num_negative))
  else
    cmap_index_negative = cmap_index_neg_start
  end if

  ; assign colors
  do i=0,num_negative-1
    ii1 = cmap_index_negative(i)
    cncolors(i,:) = cmap_negative(ii1,:)
  end do

  ii1 = 0
  ii2 = 0
  do i=num_negative,num_cncolors-1
    ii2 = cmap_index_positive(ii1)
    cncolors(i,:) = cmap_positive(ii2,:)
    ii1 = ii1+1
  end do

  ; set white color in both sides
  if (       isatt(opt_work,"white_middle") .and. opt_work@white_middle \
       .and. num_negative.ge.2 .and. abs(num_cncolors-num_negative).gt.1) then  ; add num_negative.ge.2 to avoid when there are no negative values
                                                                                ; add abs(num_cncolors-num_negative).gt.1 to avoid only one positive value
    ;*** shift negative/positive colors ***
    cncolors(0:num_negative-2,:) = cncolors(1:num_negative-1,:)
    cncolors(num_negative+1:num_cncolors-1,:) = cncolors(num_negative:num_cncolors-2,:)

    ;*** set least positive/negative to white ***
    i = num_negative-1
    cncolors(i,:) = 1.  ; white
    i = num_negative
    cncolors(i,:) = 1.  ; white
  end if

  ;*****************
  ; return NCL res
  ;*****************
  if ( isatt(res, "cnFillColors")) then
    delete(res@cnFillColors)
  end if

  if ( opt_work@cnlevels_predefine .and. isatt(res,"cnLevelSelectionMode") .and. res@cnLevelSelectionMode .eq. "ExplicitLevels" ) then
    res@cnFillColors = cncolors

  else if ( opt_work@cnlevels_predefine .and. isatt(res,"cnLevelSelectionMode") .and. res@cnLevelSelectionMode .eq. "ManualLevels" ) then
    res@cnFillColors = cncolors
      
  else
    if ( isatt(res, "cnMinLevelValF")) then
      delete(res@cnMinLevelValF)
    end if
    if ( isatt(res, "cnMaxLevelValF")) then
      delete(res@cnMaxLevelValF)
    end if
    if ( isatt(res, "cnLevelSpacingF")) then
      delete(res@cnLevelSpacingF)
    end if
  
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = cnmin
    res@cnMaxLevelValF = cnmax
    res@cnLevelSpacingF = cnint
    res@cnFillColors = cncolors
  end if
  end if

  ;**************
  ; return opt 
  ;**************
  if ( isatt(opt_cncolors,"colormap_positive") ) then
    delete(opt_cncolors@colormap_positive)
    opt_cncolors@colormap_positive = opt_work@colormap_positive
  else
    opt_cncolors@colormap_positive = opt_work@colormap_positive
  end if

  if ( isatt(opt_cncolors,"cmap_index_positive") ) then
    delete(opt_cncolors@cmap_index_positive)
    opt_cncolors@cmap_index_positive = cmap_index_positive
  else
    opt_cncolors@cmap_index_positive = cmap_index_positive
  end if

  if ( isatt(opt_cncolors,"colormap_negative") ) then
    delete(opt_cncolors@colormap_negative)
    opt_cncolors@colormap_negative = opt_work@colormap_negative
  else
    opt_cncolors@colormap_negative = opt_work@colormap_negative
  end if

  if ( isatt(opt_cncolors,"cmap_index_negative") ) then
    delete(opt_cncolors@cmap_index_negative)
    opt_cncolors@cmap_index_negative = cmap_index_negative
  else
    opt_cncolors@cmap_index_negative = cmap_index_negative
  end if

  delete(cncolors)
  delete(opt_work)
  delete(cmap_positive)
  delete(cmap_negative)
  delete(num_cmap_positive)
  delete(num_cmap_negative)

end if  ; end if of option_cnlevels

end
;##### end of cnColors_PosNeg ####

;###################################################
undef("binmap2")
function binmap2(xdata[*]:numeric, ydata:numeric, nbin:integer, opt:logical)
begin

;--------
; check
;--------
;

;--- check dimensions ---
  dimx = dimsizes(xdata)
  dimy1 = dimsizes(dimsizes(ydata))

  if (dimy1 .eq. 1) then
    dimy = dimsizes(ydata)
    dimy2 = 1
  else if (dimy1 .eq. 2) then
    dimy  = dimsizes(ydata(0,:))
    dimy2 = dimsizes(ydata(:,0))
  else 
    print("")
    print("ERROR: (binmap2) Y data must be 1 or 2 dimensional data")
    return(0)
    print("")
  end if
  end if

  if (dimx .ne. dimy) then
    print("")
    print("ERROR: (binmap2) dimension of X and Y must be the same! (dimx= "+dimx+", dimy= "+dimy+")")
    return(0)
    print("")
  end if

;--- set work data ---
  if (dimy1 .eq. 1) then
    ydata_work = new( (/1,dimy/), typeof(ydata))
    ydata_work(0,:) = ydata
  end if

  if (dimy1 .eq. 2) then
    ydata_work = ydata
  end if
  ;printVarSummary(ydata_work)

;--- set max and min values ---
  if ( isatt(opt,"x_max")) then
    x_max = opt@x_max
  else
    x_max = max(xdata)
  end if

  if ( isatt(opt,"x_min")) then
    x_min = opt@x_min
  else
    x_min = min(xdata)
  end if
  ;print(x_max+" , "+x_min)

;--- compute intervals ---
  mnmxint = nice_mnmxintvl(x_min,x_max,nbin,True)
  nlevels = toint(((mnmxint(1)-mnmxint(0))/mnmxint(2))) + 1
  levels  = fspan(mnmxint(0),mnmxint(1),nlevels)
  layers  = levels(0:nlevels-2)
  do i=0,nlevels-2
    layers(i) = levels(i) + mnmxint(2)/2.
  end do
  ;print(levels)

;--- set return variable ---
  ydata_binx = new( (/dimy2,5,nlevels-1/) , typeof(ydata))

  ydata_binx!0 = "variable"
  ydata_binx!1 = "index"
  ydata_binx!2 = "bin"

  ;ydata_binx@bin_edge = levels
  ;ydata_binx@bin_mid  = layers
  ;ydata_binx@bin_min = x_min
  ;ydata_binx@bin_max = x_max
  ;ydata_binx@bin_number_of_data = dimx
  ydata_binx@description = "Index 0: number in each bin; Index 1: Sum in each bin; Index 2: Average value in each bin; Index 3: standard deviation in each bin; Index 4: PDF of each bin"
  ydata_binx@bin_number_of_data = num(.not.ismissing(xdata))
  ydata_binx@bin_max = x_max
  ydata_binx@bin_min = x_min
  ydata_binx@bin_mid  = layers
  ydata_binx@bin_edge = levels
  ydata_binx@bin_interval = mnmxint(2)

  ydata_binx(:,0,:) = 0.  ; set the number in all bins to zero
  ydata_binx(:,4,:) = 0.  ; set the PDF to zero
  ;printVarSummary(ydata_binx)

;--- sort ydata based on x bin ---
  do i=0,nlevels-2
    ii = ind(xdata.ge.levels(i).and.xdata.lt.levels(i+1))
    if (.not.all(ismissing(ii))) then
      do j=0,dimy2-1
        ynum1 = dimsizes(ii)
        ysum1 = sum(ydata_work(j,ii))
        yavg1 = avg(ydata_work(j,ii))
        ystd1 = stddev(ydata_work(j,ii))
        ypdf1 = tofloat(ynum1)/tofloat(dimy)

        ydata_binx(j,0,i) = ynum1
        ydata_binx(j,1,i) = ysum1
        ydata_binx(j,2,i) = yavg1
        ydata_binx(j,3,i) = ystd1
        ydata_binx(j,4,i) = ypdf1
      end do
    end if
    delete(ii)
  end do

  return(ydata_binx)

end

;### end of binmap2 ###

;###################################################
undef("sort_for_xyplot")
function sort_for_xyplot(xdata[*]:numeric, ydata[*]:numeric, zdata[*]:numeric, nbin:integer)
begin

;-----------------
; check dimension
;-----------------
  dimx = dimsizes(xdata)
  dimy = dimsizes(ydata)
  dimz = dimsizes(zdata)
  if ( dimx .eq. dimy .and. dimy .eq. dimz ) then
    npts = dimx
  else
    print("")
    print("ERROR: (sort_for_xyplot) input dimensions are not the same! ("+dimx+" , "+dimy+" , "+dimz+")")
    print("")
    return(0)
  end if

;---------------------------------------------------
; sort zdata and create a array for XY plot
;   ref: scatter_4.ncl, http://www.ncl.ucar.edu/Applications/scatter.shtml
;---------------------------------------------------

;---Create roughly nbin equally spaced levels through the data.
  mnmxint = nice_mnmxintvl(min(zdata),max(zdata),nbin,True)
  nlevels = toint(((mnmxint(1)-mnmxint(0))/mnmxint(2))) + 1
  levels  = fspan(mnmxint(0),mnmxint(1),nlevels)
  layers  = levels(0:nlevels-2)
  do i=0,nlevels-2
    layers(i) = levels(i) + mnmxint(2)/2.
  end do

  ;print("min/max data = " + min(zdata) + "/" + max(zdata))
  ;print("# of groups  = " + (nlevels-1))
  ;print("levels start at = " + levels(0) + ", end at " + levels(nlevels-1) + \
  ;      ", with a spacing of " + mnmxint(2))
  ;print("")
  ;print("x= "+xdata+" , y= "+ydata+" , z= "+zdata)

;---Create new 2D array to hold groupings of values
  xplotvar = new((/nlevels-1,npts/),typeof(xdata))
  xplotvar = -999.99
  xplotvar@_FillValue = -999.99

  yplotvar = xplotvar

;---Group the values and put in 2D array.
  labels = new(nlevels-1,string)
  do i=0,nlevels-2
    ii = ind(zdata.ge.levels(i).and.zdata.lt.levels(i+1))
    if (.not.all(ismissing(ii))) then
      xplotvar(i,ii) = (/xdata(ii)/)
      yplotvar(i,ii) = (/ydata(ii)/)
    end if

    labels(i) = levels(i) + ":" + levels(i+1)
    delete(ii)
  end do

;---Set return variable
  plotvars = (/xplotvar,yplotvar/)
  plotvars!0 = "index"
  plotvars!1 = "nbin"
  if (.not.ismissing(getVarDimNames(xdata))) then
    plotvars!2 = getVarDimNames(xdata)
  else
    plotvars!2 = "size"
  end if
  plotvars@nbin_labels = labels
  plotvars@nbin_mid  = layers
  plotvars@nbin_edge = levels
  plotvars@nbin_mid  = layers
  plotvars@index = "0: x_plot_var, 1: y_plot_var"

  return(plotvars)

end

;### end of sort_for_xyplot ###

;###################################################
undef("unit_convert")
function unit_convert(var_from:numeric, unit_from:string, unit_to:string)
begin

;--------------
; check input
;--------------
  ;if units is an attribute in var_from, directly use it
  if ( unit_from .eq. "" .and. isatt(var_from,"units")) then
    unit_from = var_from@units
  end if

  if ( isatt(var_from,"units") .and. var_from@units .ne. unit_from) then
    print("")
    print("WARNING: (unit_convert) input unit ["+unit_from+"] differs from the unit attribute ["+var_from@units+"]")
    print("NOT DOING unit conversion")
    print("")
    return(var_from)
  end if

  unit_conversion = unit_from+"_TO_"+unit_to

;-------------------------
; set the return variable 
;-------------------------
  var_to = var_from
  var_to@units = unit_to
  
  flag_conversion = False
;---------------------------
; unit conversion variables 
;---------------------------
  m2mm = 1000.		; meter to milimeter
  m2cm = 100.		; meter to centimenter
  m2km = 0.001		; meter to kilometer
  mm2m = 1./m2mm	; milimeter  to meter
  cm2m = 1./m2cm	; centimeter to meter
  km2m = 1./m2km	; kilometer  to meter

  sec2day = 1./86400.	; seconds to days
  sec2hr  = 1./3600.	; seconds to hours
  hr2day  = 1./24.	; hours to days 
  day2sec = 1./sec2day	; days to seconds
  day2hr  = 1./hr2day	; hours to days
  hr2sec  = 1./sec2hr	; hours to seconds

  g2kg  = 0.001		; gram to kilogram
  kg2g  = 1./g2kg	; kilogram to gram
  mg2kg = 10.^-6	; milligram to kilogram
  kg2mg = 1./mg2kg	; kilogram to milligram

  fraction2percent = 100. ; fraction (0-1) to percent (%)

  hpa2pa  = 100.	; hPa to Pa
  pa2hpa  = 1./hpa2pa   ; Pa to hPa

  ms2knots = 1.94384    ; m/s to knots
  knots2ms = 1./ms2knots; knots to m/s

  rho_water = 1000.  ; water density, 1000 kg/m3
  latentheat_evaporation  =  2.5 * 10.^6	      ; latent heat of vaporization for water, J/kg  (http://glossary.ametsoc.org/wiki/Latent_heat)
  latentheat_condensation = 1./latentheat_evaporation ; latent heat of condensation for water, J/kg

;-----------------
; unit conversion
;-----------------
;
; *** ***
;  if (unit_conversion .eq. "") then
;    var_to = var_to
;  end if

; *** fractin ***
  if (unit_conversion .eq. "fraction_TO_percent" .or. unit_conversion .eq. "fraction_TO_%") then
    var_to = var_to * fraction2percent
    flag_conversion = True
  end if

  ;*** longitude unit conversion ***
  if (unit_conversion .eq. "degrees_west_TO_degrees_east") then
    var_to = 360. - var_to
    flag_conversion = True
  end if
  
  ;*** length ***
  if (unit_conversion .eq. "cm_TO_m") then
    var_to = var_to * cm2m
    flag_conversion = True
  end if

  if (unit_conversion .eq. "km_TO_m") then
    var_to = var_to * km2m
    flag_conversion = True
  end if

  ;*** mixing ratio ***
  if (unit_conversion .eq. "g/kg_TO_kg/kg") then
    var_to = var_to * g2kg
    flag_conversion = True
  end if

  if (unit_conversion .eq. "kg/kg_TO_g/kg") then
    var_to = var_to * kg2g
    flag_conversion = True
  end if

  if (unit_conversion .eq. "kg/kg_TO_mg/kg") then
    var_to = var_to * kg2mg
    flag_conversion = True
  end if
  
  ;*** Pressure ***
  if (unit_conversion .eq. "mb_TO_Pa") then
    var_to = var_to * hpa2pa
    flag_conversion = True
  end if

  ;*** Pressure ***
  if (unit_conversion .eq. "Pa_TO_hPa") then
    var_to = var_to * pa2hpa
    flag_conversion = True
  end if

  ;*** precipitation ***
  if (unit_conversion .eq. "m/s_TO_mm/day") then
    var_to = var_to * m2mm * day2sec
    flag_conversion = True
  end if

  if (unit_conversion .eq. "mm/day_TO_m/s") then
    var_to = var_to * mm2m * sec2day
    flag_conversion = True
  end if

  if (unit_conversion .eq. "mm/hour_TO_m/s") then
    var_to = var_to * mm2m * sec2hr
    flag_conversion = True
  end if

  if (unit_conversion .eq. "W/m2_TO_mm/day") then
    var_to = var_to * latentheat_condensation / rho_water * m2mm * day2sec
    flag_conversion = True
  end if
  
  ;*** tendency - heating rate ***
  if (unit_conversion .eq. "K/hour_TO_K/s") then
    var_to = var_to * sec2hr
    flag_conversion = True
  end if

  if (unit_conversion .eq. "K/s_TO_K/day") then
    var_to = var_to * day2sec
    flag_conversion = True
  end if

  ;*** tendency - mass mixing ratio ***
  if (unit_conversion .eq. "g/kg/hour_TO_kg/kg/s") then
    var_to = var_to * g2kg * sec2hr
    flag_conversion = True
  end if

  ;*** temperature ***
  if (unit_conversion .eq. "C_TO_K") then
    var_to = var_to + 273.15
    flag_conversion = True
  end if

; *** wind velocity *** 
  if (unit_conversion .eq. "m/s_TO_knots") then
    var_to = var_to * ms2knots
    flag_conversion = True
  end if

; *** vertical velocity *** 
  if (unit_conversion .eq. "mb/hour_TO_Pa/s") then
    var_to = var_to * hpa2pa * sec2hr
    flag_conversion = True
  end if

  if (unit_conversion .eq. "Pa/s_TO_hPa/day") then
    var_to = var_to * pa2hpa / sec2day
    flag_conversion = True
  end if

; *** water path *** 
  if (unit_conversion .eq. "kg/m2_TO_g/m2") then
    var_to = var_to * kg2g 
    flag_conversion = True
  end if

  if (unit_conversion .eq. "kg/m2_TO_mm") then
    var_to = var_to / rho_water * m2mm
    flag_conversion = True
  end if

;-------------------------------------------
; check whether executing unit conversion
;-------------------------------------------
  if (.not.flag_conversion) then
    var_to = var_from
    print("")
    print("WARNING: (unit_convert) unit conversion ["+unit_conversion+"] is unavilable")
    print("NOT DOING unit conversion")
    print("")
  end if

  return(var_to)
  
  ;printVarSummary(var_from)
  ;printVarSummary(var_to)
  ;print(var_from)
  ;print(var_to)

end

;### end of unit_convert ###

;###################################################
undef("ignorant_test")
function ignorant_test(a:numeric, b:numeric)

;----------------
; test function
;----------------

begin

  c = a + b
  return(c)

  c = a*b
  return(c)

end
;### end of ignorant_test ###

;###################################################
undef("taylor_diagram")
function taylor_diagram (wks:graphic ,RATIO[*][*]:numeric, CC[*][*]:numeric \
                                     ,rOpts:logical)
;--------------------------------------------------------------------
; This version of taylor_diagram supports "paneling"
; It requires NCL version 4.2.0.a034 because it uses "gsn_create_legend"
;--------------------------------------------------------------------

;
; Generate a Taylor Diagram:
; Generate Multiple Aspects of Model Performance in a Single Diagram
; Taylor, K. E., J. Geophys. Res., 106, D7, 7183-7192, 2001
;
; An example:
; http://www.grida.no/climate/ipcc_tar/wg1/fig8-4.htm
;
; This expects one or more datasets. The left dimension 
; is the number of datasets. The rightmost is the number of pts.
;
; Markers are at: 
; http://www.ncl.ucar.edu/Document/Graphics/Resources/gs.shtml
;
; By default, the function can handle up to 10 variable comparisons..
; To expand ...  modify the 'Colors' and 'Markers' attributes.
; The user can change / add some default settings.
;
; The defaults that the user can modify:
;
; rOpts                 = True 
;                                  ; 'made-up' resources
; rOpts@Colors          =  (/ "blue" , "red", "green", "cyan", "black" \
;                           , "turquoise", "brown", "yellow"/)
; rOpts@Markers         =  (/ 2, 3, 6, 14, 9, 12, 7, 4/) ; Marker Indices
; rOpts@markerTxOffset  = 0.0175   ; offset for text above marker
; rOpts@stnRad          = (/ 1. /) ;  (/ 0.50, 0.75, 1.5 /) 
; rOpts@centerDiffRMS   = False    ;  True mean draw additional radii from REF 
; rOpts@caseLabelsFontHeightF = 0.05
; rOpts@varLabelsFontHeightF  = 0.013
; rOpts@varLabelsYloc         = 0.65
; rOpts@legendWidth           = 0.015
; rOpts@legendHeight          = 0.030*nCase
; rOpts@taylorDraw            = True
; rOpts@taylorFrame           = True
;
;                                  ; standard NCL resources
; rOpts@tiMainString    = "Taylor" ; not using title makes plot bigger
; rOpts@gsMarkerSizeF   = 0.0085   ; marker size   
; rOpts@gsMarkerThicknessF = 1.0
; rOpts@txFontHeightF   = 0.0125   ; text size 
; rOpts@tiMainFontHeightF = 0.0225 ; tiMainString size
;
; It returns to the user a graphic object containing the 
; Taylor background and plotted x/y pts.
; This graphic object contains a simple Taylor background appropriate
; for standardized data and the markers for the datasets.
; ==================================================================
; This version allows paneling:
;      The 'cumbersome' "dum" variables were added by 
;      Adam Phillips to allow paneling via "gsn_add_?".
; ==================================================================
begin
  dimR                  = dimsizes(RATIO)
  nCase                 = dimR(0)    ; # of cases [models] 
  nVar                  = dimR(1)    ; # of variables

                                     ; x/y coordinates for plotting
  X    = new ( (/nCase,nVar/) , typeof(RATIO) )
  Y    = new ( (/nCase,nVar/) , typeof(RATIO) )

  do nc=0,nCase-1
     angle      = acos( CC(nc,:) )   ; array operation                                    
     X(nc,:)    = RATIO(nc,:)*cos( angle )     
     Y(nc,:)    = RATIO(nc,:)*sin( angle )    
  end do

  xyMin                 = 0.  
  xyOne                 = 1.00
  xyMax                 = 1.65
  xyMax_Panel           = xyMax+ 0.10            ; paneling purposes
 
  if (rOpts .and. isatt(rOpts,"txFontHeightF"))  then 
      FontHeightF       = rOpts@txFontHeightF    ; user wants to specify size
  else
      FontHeightF       = 0.0175
  end if
 
; ----------------------------------------------------------------
; Part 1:
; base plot: Based upon request of Mark Stevens
; basic x-y and draw the 1.0 observed and the outer curve at 1.65
; ----------------------------------------------------------------
  
  rxy                   = True       
  rxy@gsnDraw           = False
  rxy@gsnFrame          = False
  rxy@vpHeightF         = 0.65
  rxy@vpWidthF          = 0.65
  rxy@tmYLBorderOn      = False
  rxy@tmXBBorderOn      = False

  rxy@tiYAxisString     = "Standardized Deviations (Normalized)"
  rxy@tiYAxisFontHeightF= FontHeightF                        ; default=0.025 
  
  rxy@tmXBMode          = "Explicit" 
  rxy@tmXBValues        = (/0.0,0.25,0.50,0.75,1.00,1.25,1.5/)    ; major tm
                                                                  ; default  "OBS" or "REF"
 ;rxy@tmXBLabels        = (/"0.00","0.25","0.50","0.75","REF" ,"1.25","1.50"/)
  rxy@tmXBLabels        = (/"    ","0.25","0.50","0.75","REF" ,"1.25","1.50"/)
  if (rOpts .and. isatt(rOpts,"OneX") )  then                     ; eg: rOpts@OneX="1.00" 
     ;rxy@tmXBLabels        = (/"0.00","0.25","0.50","0.75",rOpts@OneX,"1.25","1.50"/)
      rxy@tmXBLabels        = (/"    ","0.25","0.50","0.75",rOpts@OneX,"1.25","1.50"/)
  end if

  rxy@tmXBMajorLengthF  = 0.015      ; default=0.02 for a vpHeightF=0.6
  rxy@tmXBLabelFontHeightF = FontHeightF
  rxy@tmXBMinorOn       = False
  rxy@trXMaxF           = xyMax_Panel

  rxy@tmYLMode          = "Manual"
  rxy@tmYLMinorOn       = False
  rxy@tmYLMajorLengthF  = rxy@tmXBMajorLengthF
  rxy@tmYLLabelFontHeightF = FontHeightF
  rxy@tmYLMode          = "Explicit" 
  rxy@tmYLValues        = (/0.0, .25,0.50, 0.75, 1.00, 1.25, 1.5/) ; major tm
  rxy@tmYLLabels        = (/"0.00","0.25","0.50","0.75","1.00","1.25","1.50"/)
 ;rxy@tmYLLabels        = (/"    ","0.25","0.50","0.75","1.00","1.25","1.50"/)
  rxy@trYMaxF           = xyMax_Panel

  rxy@tmYRBorderOn      = False
  rxy@tmYROn            = False      ; Turn off right tick marks.

  rxy@tmXTBorderOn      = False
  rxy@tmXTOn            = False      ; Turn off right tick marks.

  rxy@xyDashPatterns    = (/ 0 /)    ; line characteristics (dash,solid)
  rxy@xyLineThicknesses = (/ 2./)    ; choose line thickness

  rxy@gsnFrame          = False      ; Don't advance the frame.

                                            ; create outer 'correlation axis'
  npts    = 200                        ; arbitrary
  xx      = fspan(xyMin,xyMax,npts) 
  yy      = sqrt(xyMax^2 - xx^2    )   ; outer correlation line (xyMax)

  sLabels = (/"0.0","0.1","0.2","0.3","0.4","0.5","0.6" \ ; correlation labels
             ,"0.7","0.8","0.9","0.95","0.99","1.0"     /); also, major tm
  cLabels = stringtofloat(sLabels)
  rad     = 4.*atan(1.0)/180.
  angC    = acos(cLabels)/rad                     ; angles: correlation labels
                                                                       
  if (rOpts .and. isatt(rOpts,"tiMainString")) then
      rxy@tiMainString      = rOpts@tiMainString
     ;rxy@tiMainOffsetYF    = 0.015               ; default  0.0
      if (isatt(rOpts,"tiMainFontHeightF")) then
           rxy@tiMainFontHeightF = rOpts@tiMainFontHeightF
      else
           rxy@tiMainFontHeightF = 0.0225         ; default  0.025              
      end if
  end if
;;if (rOpts .and. isatt(rOpts,"gsnCenterString")) then
;;    rxy@gsnCenterString  = rOpts@gsnCenterString      ; only gsn_csm_xy
;;end if

  taylor  = gsn_xy(wks,xx,yy,rxy)                 ; Create and draw XY plot.

  rsrRes  = True
  rsrRes@gsLineThicknessF  = rxy@xyLineThicknesses(0)  ; line thickness
  rsrRes@gsLineDashPattern = 0                    ; solid line pattern
                                                  ; draw x and y to xyMax
  dum0 = gsn_add_polyline(wks,taylor,(/0.,  0. /),(/0.,xyMax/), rsrRes)
  dum1 = gsn_add_polyline(wks,taylor,(/0.,xyMax/),(/0.,  0. /), rsrRes)

  xx   = fspan(xyMin, xyOne ,npts)                ; draw 1.0 standard radius
  yy   = sqrt(xyOne - xx^2)   
  rsrRes@gsLineDashPattern = 1                    ; dashed line pattern
  rsrRes@gsLineThicknessF  = rxy@xyLineThicknesses(0)  ; line thickness
  dum2 = gsn_add_polyline(wks,taylor,xx,yy, rsrRes)
  delete(xx)
  delete(yy)
                                                  
  if (rOpts .and. isatt(rOpts,"stnRad") ) then
      rsrRes@gsLineThicknessF  = 1   ; rxy@xyLineThicknesses(0)  
      nStnRad = dimsizes(rOpts@stnRad)

      dum3  = new(nStnRad,graphic)
      do n=0,nStnRad-1
         rr = rOpts@stnRad(n)
         xx = fspan(xyMin, rr ,npts) 
         yy = sqrt(rr^2   - xx^2)   
         dum3(n) = gsn_add_polyline(wks,taylor,xx,yy, rsrRes)
      end do
      taylor@$unique_string("dum")$ = dum3

      delete(xx)
      delete(yy)
  end if

  getvalues taylor                                ; get style info from taylor
    "tmYLLabelFont"        : tmYLLabelFont        ; use for correlation axis
    "tmYLLabelFontHeightF" : tmYLLabelFontHeightF
  end getvalues

; ----------------------------------------------------------------
; Part 2:
; Correlation labels
; ----------------------------------------------------------------
  radC    = xyMax                                  ; for correlation labels
  xC      = radC*cos(angC*rad)
  yC      = radC*sin(angC*rad)
; added to get some separation
  xC      = xC + 0.020*cos(rad*angC)
  yC      = yC + 0.060*sin(rad*angC)

  txRes               = True                      ; text mods desired
  txRes@txFontHeightF = FontHeightF               ; match YL 
  txRes@tmYLLabelFont = tmYLLabelFont             ; match YL
  txRes@txAngleF      = -45.
  if (.not.isatt(rOpts,"drawCorLabel") .or. rOpts@drawCorLabel) then 
      dum4 = gsn_add_text(wks,taylor,"Correlation",1.30,1.30,txRes)
	 taylor@$unique_string("dum")$ = dum4
  end if
  txRes@txAngleF      = 0.0 
  txRes@txFontHeightF = FontHeightF*0.50          ; bit smaller

;;dum0 = gsn_add_text(wks,taylor,"OBSERVED",1.00,0.075,txRes)

  plRes               = True
  plRes@gsLineThicknessF = 2.
  
  txRes@txJust        = "CenterLeft"              ; Default="CenterCenter".
  txRes@txFontHeightF = FontHeightF               ; match YL 
 ;txRes@txBackgroundFillColor = "white"

  tmEnd = 0.975
  radTM = xyMax*tmEnd                             ; radius end: major TM 
  xTM   = new( 2 , "float")
  yTM   = new( 2 , "float")

  dum5 = new(dimsizes(sLabels),graphic)
  dum6 = dum5

  do i=0,dimsizes(sLabels)-1                      ; Loop to draw strings
    txRes@txAngleF = angC(i)
    dum5(i) = gsn_add_text(wks, taylor, sLabels(i),xC(i),yC(i),txRes) ; cor label
    xTM(0)   = xyMax*cos(angC(i)*rad)             ; major tickmarks at
    yTM(0)   = xyMax*sin(angC(i)*rad)             ; correlation labels
    xTM(1)   = radTM*cos(angC(i)*rad)             
    yTM(1)   = radTM*sin(angC(i)*rad)
    dum6(i) = gsn_add_polyline(wks,taylor,xTM,yTM,plRes)
  end do
                                                  ; minor tm locations
  mTM     = (/0.05,0.15,0.25,0.35,0.45,0.55,0.65 \ 
             ,0.75,0.85,0.91,0.92,0.93,0.94,0.96,0.97,0.98  /)
  angmTM  = acos(mTM)/rad                         ; angles: correlation labels
  radmTM  = xyMax*(1.-(1.-tmEnd)*0.5)             ; radius end: minor TM 

  dum7 = new(dimsizes(mTM),graphic)

  do i=0,dimsizes(mTM)-1                          ; manually add tm
    xTM(0)   = xyMax*cos(angmTM(i)*rad)           ; minor tickmarks
    yTM(0)   = xyMax*sin(angmTM(i)*rad)
    xTM(1)   = radmTM*cos(angmTM(i)*rad)          
    yTM(1)   = radmTM*sin(angmTM(i)*rad)
    dum7(i)  = gsn_add_polyline(wks,taylor,xTM,yTM,plRes)
  end do
                                                  ; added for Wanli
  if (rOpts .and. isatt(rOpts,"ccRays") ) then
      angRL = acos(rOpts@ccRays)/rad             ; angles: radial lines

      rlRes = True
      rlRes@gsLineDashPattern= 2  ; line pattern
      rlRes@gsLineThicknessF = 1  ; choose line thickness
      if (isatt(rOpts,"ccRays_color")) then
          rlRes@gsLineColor    =  "LightGray"
      end if

      dum8 = new(dimsizes(angRL),graphic)
      do i=0,dimsizes(angRL)-1
         xRL     = xyMax*cos(angRL(i)*rad)
         yRL     = xyMax*sin(angRL(i)*rad)
         dum8(i) = gsn_add_polyline(wks,taylor,(/0, xRL /),(/0,  yRL  /),rlRes)
      end do
      taylor@$unique_string("dum")$ = dum8
  end if
  
; ----------------------------------------------------------------
; Part 3:
; Concentric about 1.0 on XB axis
; ----------------------------------------------------------------
  if (rOpts .and. isatt(rOpts,"centerDiffRMS") \
            .and. rOpts@centerDiffRMS) then
      respl                    = True                ; polyline mods desired
      respl@xyLineThicknessF   = 1.0                 ; line thickness
      respl@xyLineDashPattern  = 2                   ; short dash lines
      respl@gsLineColor        = "Black"             ; line color     
      if (isatt(rOpts,"centerDiffRMS_color")) then
          respl@gsLineColor    =  "LightGray"
      end if
      
      dx   = 0.25
      ncon = 4                                       ; 0.75, 0.50, 0.25, 0.0
      npts = 100                                     ; arbitrary
      ang  = fspan(180,360,npts)*rad

      dum9 = new(ncon,graphic)

      do n=1,ncon 
         rr  = n*dx            ; radius from 1.0 [OBS] abscissa
         xx  = 1. + rr*cos(ang)
         yy  = fabs( rr*sin(ang) )
         if (n.le.2) then
             dum9(n-1) = gsn_add_polyline(wks,taylor,xx,yy,respl)
         end if
         if (n.eq.3) then
             n3 = floattointeger( 0.77*npts ) 
             dum9(n-1) = gsn_add_polyline(wks,taylor,xx(0:n3),yy(0:n3),respl)
         end if
         if (n.eq.4) then
             n4 = floattointeger( 0.61*npts ) 
             dum9(n-1) = gsn_add_polyline(wks,taylor,xx(0:n4),yy(0:n4),respl)
         end if
      end do
      delete(ang)
      delete(xx)
      delete(yy)
      taylor@$unique_string("dum")$ = dum9

  end if
; ---------------------------------------------------------------
; Part 4:
; generic resources that will be applied to all users data points
; of course, these can be changed 
; http://www.ncl.ucar.edu/Document/Graphics/Resources/gs.shtml
; ---------------------------------------------------------------
  if (rOpts .and. isatt(rOpts,"Markers")) then
      Markers = rOpts@Markers
  else
      Markers = (/ 4, 6, 8,  0, 9, 12, 7, 2, 11, 16/) ; Marker Indices
  end if

  if (rOpts .and. isatt(rOpts,"Colors")) then
      Colors  = rOpts@Colors
  else
      Colors  = (/ "red", "blue", "green", "cyan", "orange" \
                 , "turquoise", "brown", "yellow", "purple", "black"/)
  end if

  if (rOpts .and. isatt(rOpts,"gsMarkerThicknessF")) then
      gsMarkerThicknessF = rOpts@gsMarkerThicknessF
  else
      gsMarkerThicknessF = 1.0
  end if

  if (rOpts .and. isatt(rOpts,"gsMarkerSizeF")) then
      gsMarkerSizeF      = rOpts@gsMarkerSizeF
  else
      gsMarkerSizeF      = 0.0085                  ; Default: 0.007
  end if

  gsRes = True
  gsRes@gsMarkerThicknessF = gsMarkerThicknessF      ; default=1.0
  gsRes@gsMarkerSizeF      = gsMarkerSizeF           ; Default: 0.007 

  ptRes = True                        ; text options for points
  ptRes@txJust             = "BottomCenter"; Default="CenterCenter".
  ptRes@txFontThicknessF   = 1.2      ; default=1.00
  ptRes@txFontHeightF      = 0.0125   ; default=0.05
  if (rOpts .and. isatt(rOpts,"txFontHeightF")) then
      ptRes@txFontHeightF  = rOpts@txFontHeightF  
  end if

  markerTxYOffset          = 0.0175   ; default
  if (rOpts .and. isatt(rOpts,"markerTxYOffset")) then
      markerTxYOffset = rOpts@markerTxYOffset             ; user defined offset
  end if

  dum10 = new((nCase*nVar),graphic)
  dum11 = dum10

  do n=0,nCase-1
     gsRes@gsMarkerIndex   = Markers(n)             ; marker style (+)
     gsRes@gsMarkerColor   = Colors(n)              ; marker color
     ptRes@txFontColor     = gsRes@gsMarkerColor
    do i=0,nVar-1
       dum10(n*nVar+i) = gsn_add_polymarker(wks,taylor,X(n,i),Y(n,i),gsRes) 
       dum11(n*nVar+i) = gsn_add_text(wks,taylor,(i+1),X(n,i),Y(n,i)+markerTxYOffset,ptRes)
    end do
  end do

; ---------------------------------------------------------------
; Part 5: ; add case legend and variable labels
; ---------------------------------------------------------------

  if (rOpts .and. isatt(rOpts,"caseLabels")) then 

      if (isatt(rOpts,"caseLabelsFontHeightF")) then
          caseLabelsFontHeightF = rOpts@caseLabelsFontHeightF
      else
          caseLabelsFontHeightF = 0.05  
      end if

      lgres                    = True
      lgres@lgMarkerColors     = Colors        ; colors of markers
      lgres@lgMarkerIndexes    = Markers       ; Markers 
      lgres@lgMarkerSizeF      = gsMarkerSizeF ; Marker size
      lgres@lgItemType         = "Markers"     ; draw markers only
      lgres@lgLabelFontHeightF = caseLabelsFontHeightF  ; font height of legend case labels

      if (isatt(rOpts,"legendWidth")) then
          lgres@vpWidthF       = rOpts@legendWidth
      else
          lgres@vpWidthF       = 0.15           ; width of legend (NDC)
      end if

      if (isatt(rOpts,"legendHeight")) then
          lgres@vpHeightF      = rOpts@legendHeight
      else   
          lgres@vpHeightF      = 0.030*nCase   ; height of legend (NDC)
      end if

      lgres@lgPerimOn          = False         ; turn off perimeter
      nModel                   = dimsizes( rOpts@caseLabels )
      lbid = gsn_create_legend(wks,nModel,rOpts@caseLabels,lgres)
	 
      amres = True
      amres@amParallelPosF     =  0.35           
      amres@amOrthogonalPosF   = -0.35             
      annoid1 = gsn_add_annotation(taylor,lbid,amres)	; add legend to plot
  end if

  if (rOpts .and. isatt(rOpts,"varLabels")) then 
      nVar    = dimsizes(rOpts@varLabels)

      if (isatt(rOpts,"varLabelsFontHeightF")) then
          varLabelsFontHeightF = rOpts@varLabelsFontHeightF
      else
          varLabelsFontHeightF = 0.013
      end if

      txres = True
      txres@txFontHeightF = varLabelsFontHeightF
      txres@txJust = "CenterLeft"              ; justify to the center left

     ;delta_y = 0.02       
     delta_y = 0.06   
      if (rOpts .and. isatt(rOpts,"varLabelsYloc")) then
          ys  = rOpts@varLabelsYloc            ; user specified
      else
          ys  = max( (/nVar*delta_y , 0.30/) )
      end if

      
      do i = 1,nVar     
         if (i.eq.1) then
             dum12 = new(nVar,graphic)
	 end if

         dum12(i-1) = gsn_add_text(wks,taylor,i+" - "+rOpts@varLabels(i-1), .125,ys,txres)
         ys = ys- delta_y
      end do

      taylor@$unique_string("dum")$ = dum12
  end if

  taylor@$unique_string("dum")$ = dum0   ; x-axis
  taylor@$unique_string("dum")$ = dum1   ; y-axis
  taylor@$unique_string("dum")$ = dum2   ; 1.0 std curve
  taylor@$unique_string("dum")$ = dum5   ; labels [COR]
  taylor@$unique_string("dum")$ = dum6   ; major tm [COR]
  taylor@$unique_string("dum")$ = dum7   ; minor tm
  taylor@$unique_string("dum")$ = dum10  ; markers
  taylor@$unique_string("dum")$ = dum11  ; text
  
  if (.not.isatt(rOpts,"taylorDraw") .or. \
     (isatt(rOpts,"taylorDraw") .and. rOpts@taylorDraw)) then 
	draw(taylor)
  end if
  if (.not.isatt(rOpts,"taylorFrame") .or. \
     (isatt(rOpts,"taylorFrame") .and. rOpts@taylorFrame)) then 
	frame(wks)
  end if
  return(taylor)
end

;### end of taylor_diagram ###

;###################################################
undef("binmap")
function binmap(xx:numeric,yy:numeric,xx_min[1]:numeric,xx_max[1]:numeric,xx_int[1]:numeric)

begin 

;--------------
; check part
;--------------

; check dimension sizes
  xx_dim = dimsizes(xx)
  yy_dim = dimsizes(yy)

  dd = dimsizes(xx_dim)

  do i=0,dd-1
    if (xx_dim(i).ne.yy_dim(i)) then
      print("")
      print("ERROR: (func binmap) variable dimension index [" + i + "] does not match [" + xx_dim(i) + "," + yy_dim(i) + "]")
      print("function binmap stop")
      return(0)
    end if
  end do

; check types
  xx_type = typeof(xx)
  yy_type = typeof(yy)

  if (xx_type.ne.yy_type) then
    print("")
    print("ERROR: (func binmap) variable type [" + xx_type + "] does not match [" + yy_type + "]")
    print("function binmap stop")
    return(0)
  end if

; check min,max,interval
  if (xx_min.ge.xx_max) then
    print("")
    print("ERROR: (func binmap) MIN value [" + xx_min + "] is greater than MAX value [" + xx_max + "]")
    print("function binmap stop")
    return(0)
  else if (xx_int.gt.xx_max-xx_min) then
    print("")
    print("ERROR: (func binmap) INTERVAL value [" + xx_int + "] is greater than the value range [" + xx_min + "," + xx_max + "]")
    print("function binmap stop")
    return(0)
  end if
  end if

;----------------------
; calculate bin for xx
;----------------------

;*** merge multidiemsional input array into 1-D temp array ***
  length = 1
  do i=0,dd-1
    length = length * xx_dim(i)  ; new 1-D array length
  end do

  xtemp  = new(length, xx_type)
  ytemp  = new(length, yy_type)
  ytemp0 = new(length, yy_type)

  if (isatt(xx,"_FillValue")) then
    xtemp@_FillValue  = xx@_FillValue
    ytemp@_FillValue  = xx@_FillValue
    ytemp0@_FillValue = xx@_FillValue
  else
    xtemp@_FillValue  = -999.99
    ytemp@_FillValue  = -999.99
    ytemp0@_FillValue = -999.99
  end if

  n=0

  if (dd.eq.1) then
    do i=0,xx_dim(0)-1
      if (.not.ismissing(xx(i))) then
        xtemp(n) = xx(i)
        ytemp0(n) = yy(i)
        n=n+1
      end if
    end do
  else if (dd.eq.2) then
    do i=0,xx_dim(0)-1
    do j=0,xx_dim(1)-1
      if (.not.ismissing(xx(i,j)) .and. .not.ismissing(yy(i,j))) then
        xtemp(n) = xx(i,j)
        ytemp0(n) = yy(i,j)
        n=n+1
      end if
    end do
    end do

  else if (dd.eq.3) then
    do i=0,xx_dim(0)-1
    do j=0,xx_dim(1)-1
    do k=0,xx_dim(2)-1
      if (.not.ismissing(xx(i,j,k)) .and. .not.ismissing(yy(i,j,k))) then
        xtemp(n) = xx(i,j,k)
        ytemp0(n) = yy(i,j,k)
        n=n+1
      end if
    end do
    end do
    end do

  else if (dd.eq.4) then
    do i=0,xx_dim(0)-1
    do j=0,xx_dim(1)-1
    do k=0,xx_dim(2)-1
    do t=0,xx_dim(3)-1
      if (.not.ismissing(xx(i,j,k,t)) .and. .not.ismissing(yy(i,j,k,t))) then
        xtemp(n) = xx(i,j,k,t)
        ytemp0(n) = yy(i,j,k,t)
        n=n+1
      end if
    end do
    end do
    end do
    end do

  end if
  end if
  end if
  end if

;--------------
; bin input xx
;--------------

; output arrays
  num_z = toint((xx_max-xx_min)/xx_int)

  zz_bin_int  = new(num_z+1,xx_type)  ; bin value at interface

  zz_bin_mid   = new(num_z,xx_type)   ; bin value at midpoint
  if (isatt(xx,"_FillValue")) then
    zz_bin_mid@_FillValue = xx@_FillValue
  else
    zz_bin_mid@_FillValue = -999.99
  end if

  zz_mean = zz_bin_mid  ; bin average value
  zz_std  = zz_bin_mid  ; bin standard deviation

  zz_count = new(num_z,"integer")  ; bin number count
  zz_count = 0

; set bin values at interface
  zz_bin_int(0)       = xx_min
  zz_bin_int(num_z) = xx_max

  do i=1,num_z-1
    zz_bin_int(i) = zz_bin_int(i-1) + xx_int
  end do

; set bin midpoint value
  do i=0,num_z-1
    zz_bin_mid(i) = (zz_bin_int(i)+zz_bin_int(i+1))/2.
  end do

; sort x and get indexes
  kflg = 2
  idx = dim_pqsort(xtemp(0:n-1),kflg)  ; sort xx in increasing order and get indexes

; count number in each bin 
  do i=0,n-1
    if (.not.ismissing(xtemp(i)) .and. xtemp(i) .le. zz_bin_int(num_z)) then
      j = closest_val(xtemp(i),zz_bin_mid)
      if (j .ge. num_z) then
        j=num_z-1
      end if
      zz_count(j) = zz_count(j)+1
    end if
    print(i + "/" + n + " bin processing...")
  end do

; assign yy into each bin
  do i=0,n-1
    ytemp(i) = ytemp0(idx(i))
  end do

; compute mean and standrad deviation of the first bin
  i=0
  if (zz_count(i).eq.0) then
    zz_mean(i) = zz_bin_mid@_FillValue
    zz_std (i) = zz_bin_mid@_FillValue
  else
    i1=0
    i2=0+zz_count(i)-1
    zz_mean(i) = sum(ytemp(i1:i2))/zz_count(i)
    zz_std (i) = stddev(ytemp(i1:i2))
  end if

; compute mean and standard deviation at each bin
  do i=1,num_z-1

    ; if no value in the bin
    if (zz_count(i).eq.0) then
      zz_mean(i) = zz_bin_mid@_FillValue
      zz_std (i) = zz_bin_mid@_FillValue

    else
      i1 = sum(zz_count(0:i-1))
      i2 = i1 + zz_count(i) - 1
      zz_mean(i) = sum(ytemp(i1:i2))/zz_count(i)
      zz_std (i) = stddev(ytemp(i1:i2))
    end if

  end do

;----------------
; return arrays
;----------------

  zz_all = new((/3,num_z/),xx_type) 

;  zz_all(0,:) = zz_bin_mid(:)

  zz_all(0,:) = zz_count  (:)
  zz_all(1,:) = zz_mean   (:)
  zz_all(2,:) = zz_std    (:)

  ; coordinate variable
  zz_all!0 = "index"
  zz_all!1 = "bin"
  zz_all&bin = zz_bin_mid

  ; attributes 
  zz_all@lingname = "Sort input data to bin"
  zz_all@description = "Index 0: number in each bin; Index 1: Average value in each bin; Index 2: standard deviation in each bin"
  zz_all@bin_max      = xx_max
  zz_all@bin_min      = xx_min
  zz_all@bin_interval = xx_int
  zz_all@bin_total_number = sum(zz_count(:))
  zz_all@bin_interface = zz_bin_int
  zz_all@pdf = tofloat(zz_count)/tofloat(zz_all@bin_total_number)

  return(zz_all)
 
end

;### end of bimmap ###

;###################################################
undef("lat_wgt_avg")
function lat_wgt_avg(data_in:numeric)

; weighting: 
; https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2004/msg00023.html
; http://mathforum.org/dr.math/faq/formulas/faq.sphere.html

begin

  pi = acos(-1.)
  lat = data_in&lat  ; latitude must in -90 to 90

  ;coslat = abs(cos(lat/180.*pi)) ; wrong weighting
  gwlat = latRegWgt(lat,typeof(lat),0)  ; no gw in file

  ; function wgt_areaave (q[(...),y,x]:numeric, wgty [*]:numeric, wgtx[*]:numeric, 
  ;                       opt:integer, 0: calculate only non-missing data; 1:return missing value if there is any missing data)
  ; ** yihsuan note 2019-01-24 *** : opt=0,1 works for wgt_areaave, but NOT work for wgt_areaave_Wrap
  data_out = wgt_areaave_Wrap(data_in,gwlat,1.0,1)

  ;data_out00 = wgt_areaave(data_in,coslat,1.0,1)  ; if return missing value if there is any missing data
  ;data_out = (/data_out00/)

  if (isatt(data_in,"long_name")) then
    data_out@long_name = data_in@long_name ;+ " (latitude-weighted area mean)"
  ;else
  ;  data_out@long_name = "latitude-weighted area mean variable"
  end if

  return(data_out)

;*** brute force way ***
;  pi = acos(-1.)
;  lat = data_in&lat  ; latitude must in -90 to 90
;  coslat = tofloat(abs(cos(lat/180.*pi)))
;
;  dim = dimsizes(dimsizes(data_in))
;
;  print("num of missing data: "+num(ismissing(data_in)))
;
;  if (dim.eq.3) then
;    sum_datain = data_in(:,0,0)
;    sum_datain = 0.
;
;    avg_datain = sum_datain
;
;    nyy = dimsizes(data_in(0,:,0))
;    nxx = dimsizes(data_in(0,0,:))
;    nn1 = dimsizes(data_in(:,0,0))
;
;    do n1=0,nn1-1
;      count = 0
;      do xx=0,nxx-1
;      do yy=0,nyy-1
;        if (.not.ismissing(data_in(n1,yy,xx))) then
;          sum_datain(n1) = sum_datain(n1) + data_in(n1,yy,xx)*coslat(yy)
;          count = count+1
;        end if
;      end do ; end loop of yy
;      end do ; end loop of xx
;
;      avg_datain(n1) = sum_datain(n1)/tofloat(count)
;    end do ; end loop of nn1
;
;    print("count= "+count)
;  end if  ; end if of dim=3
;
;  ;print(avg_datain)
;  data_out = avg_datain
;*** brute force way ***


end
;### end of lat_wgt_avg ###

;###################################################
undef("plus")
function plus(a:numeric, b:numeric)

begin

  a_plus_b = a
  a_plus_b = a + b
  return(a_plus_b)

end
;### end of plus ###

;###################################################
undef("divide")
function divide(a:numeric, b:numeric)

begin

  a_divide_b = a
  a_divide_b = a / b
  return(a_divide_b)

end
;### end of divide ###

;###################################################
undef("multiply")
function multiply(a:numeric, b[1]:numeric)

begin

  a_multi_b = a
  a_multi_b = a * b
  return(a_multi_b)

end
;### end of multiply ###

;###################################################
undef("diff")
function diff(a:numeric, b:numeric)

begin

  a_minus_b = a
  a_minus_b = a - b
  return(a_minus_b)

end
;### end of diff ###

;###################################################
undef("diff_group")
function diff_group(var_in:numeric, option_group[1]:logical)

begin

  ;*** get var_in dimension ***
  nn1 = dimsizes(var_in)
  nfile = nn1(0)
  ndim = dimsizes(dimsizes(var_in))
  ;print(ndim)

  ;*** set base index array ***
  base_idx = new(nfile, "integer")
  ;base_idx = ispan(0, nfile-1, 1)
  ;print(base_idx)

  ;*** set base var ***
  var_base = var_in
  var_base = 0.
  
  ;if (option_group) then

    if (isatt(option_group,"nfile_per_group") .and. isint(option_group@nfile_per_group)) then
      nfile_per_group = option_group@nfile_per_group
    else
      print("ERROR (diff_group): option_group@nfile_per_group is either not set or is not an integer")
      print("Program stop")
      exit()
    end if

   ; make sure the number of group is an integer
    if ( nfile_per_group.eq.nfile .or. mod(nfile,nfile_per_group).ne.0 ) then
      print("ERROR (diff_group): nfile is ["+nfile+"] and nfile_per_group is ["+nfile_per_group+"], # of group MUST be an interger and not equal to nfile")
      print("Program stop")
      exit()
    end if

    ;*** get base index ***
    do i=0,nfile-1,nfile_per_group
      ii = i+nfile_per_group-1
      base_idx(i:ii) = ii
    end do
    ;print(base_idx)

    ;*** set base var ***
    do i=0,nfile-1
      bb1 = base_idx(i)
      if (ndim.eq.1) then
        var_base(i) = var_in(bb1)
      end if
      if (ndim.eq.2) then
        var_base(i,:) = var_in(bb1,:)
      end if
      if (ndim.eq.3) then
        var_base(i,:,:) = var_in(bb1,:,:)
      end if
      if (ndim.eq.4) then
        var_base(i,:,:,:) = var_in(bb1,:,:,:)
      end if
      if (ndim.eq.5) then
        var_base(i,:,:,:,:) = var_in(bb1,:,:,:,:)
      end if
      if (ndim.eq.6) then
        var_base(i,:,:,:,:,:) = var_in(bb1,:,:,:,:,:)
      end if
    end do
    ;print(var_in+","+var_base)

    ;*** minus base var ***
    var_inMbase = diff(var_in, var_base)
    var_inMbase@nfile_per_group = nfile_per_group
    var_inMbase@base_index = base_idx
    ;print(var_in+" MINUS "+var_base+" = "+var_inMbase)

    return(var_inMbase)

  ;end if ; end if of option_group

end
;### end of diff_group ###

;###################################################
undef("gsn_csm_contour_map_avg")
function gsn_csm_contour_map_avg(wks[1]:graphic, data_input0[*][*]:numeric, \
                                 nclres[1]:logical, opt[1]:logical)

begin

;---------------------------------------
; read data_input in certain lat/lon
;---------------------------------------
data_input = get_VarLatLon(data_input0,nclres)

;-------------
; initialize
;-------------
mean_format = "%7.3f"
mean_prefix = "mean: "
LatWeightCheck = True
pi = acos(-1.)

;-------------------
; set mean format
;-------------------
if ( isatt(opt, "MeanFormat")) then
  mean_format = opt@MeanFormat
end if

if ( isatt(opt, "MeanPrefix")) then
  mean_prefix = opt@MeanPrefix
end if

if ( isatt(opt, "LatWeighted")) then 
  LatWeightCheck = opt@LatWeighted
;else
;  LatWeightCheck = LatWeightCheck_default
end if

;--------------------------------------
; check input data has right lat coordinate
;--------------------------------------

if ( LatWeightCheck ) then
  if ( iscoord(data_input, "lat") .and. data_input&lat@units .eq. "degrees_north" ) then
      LatWeightCheck = True
      num_lat = dimsizes(data_input&lat)
      lat_temp = data_input&lat

  else if ( iscoord(data_input, "latitude") .and. data_input&latitude@units .eq. "degrees_north" ) then
      LatWeightCheck = True
      num_lat = dimsizes(data_input&latitude)
      lat_temp = data_input&latitude

  else
    LatWeightCheck = False

    print("")
    print("WARNING: (gsn_csm_contour_map_avg) Invalid coordinate variable's setup. Its name MUST be [lat OR latitude] and units be [degrees_north]")
    print("")

  end if    
  end if  
end if

;--------------------------
; compute global average
;   1. DO     latitude correction
;   2. NOT DO latitude correction
;--------------------------

;
;*** DO latitude correction ***
;
if (LatWeightCheck) then

  coslat = abs(cos(lat_temp/180.*pi))
  ; function wgt_areaaveq (q[(...),y,x]:numeric, wgty [*]:numeric, wgtx[*]:numeric, 
  ;                        opt:integer, 0: calculate only non-missing data; 1:return missing value if there is any missing data)
  data_input_avg = wgt_areaave(data_input,coslat,1.0,1)
;
;*** NOT DO latitude correction ***
;
else

  data_input_avg = wgt_areaave(data_input,1.0,1.0,1)

end if

;------------------------
; set ncl res and plot
;------------------------
if ( isatt(nclres,"gsnCenterString") ) then
  wres = nclres

else
  wres = nclres
  wres@gsnCenterString = (mean_prefix + sprintf(mean_format,data_input_avg))

end if

;-------
; plot
;-------

if (isatt(opt,"mapview") .and. opt@mapview .eq. "polar_NH") then
  wres@gsnPolar   = "NH"
  plot = gsn_csm_contour_map_polar(wks,data_input,wres)

else if (isatt(opt,"mapview") .and. opt@mapview .eq. "polar_SH") then
  wres@gsnPolar   = "SH"
  plot = gsn_csm_contour_map_polar(wks,data_input,wres)

else
  plot = gsn_csm_contour_map(wks,data_input,wres)

end if
end if

;--------
; return 
;--------
return(plot)

end

;### end of gsn_csm_contour_map_avg ###

;######################################################
undef("ignorant_combine_colormaps")
function ignorant_combine_colormaps( cmap1[*][3]:numeric, cmap2[*][3]:numeric, \
                                     choice:integer  )

begin

;---------------------------
; save input cmap dimension
;----------------------------
  cmap1_dim = dimsizes(cmap1(:,0))
  cmap2_dim = dimsizes(cmap2(:,0))

;-----------------------------
; determine output color map
;-----------------------------

;*** combine total cmap1 & cmap2 *** 

  if (choice .eq. 1) then

    ;---------------------
    ; create output cmap
    ;---------------------
      ocamp_dim = cmap1_dim + cmap2_dim
      ocmap = new((/ocamp_dim,3/),"float")

    ;---------------------------
    ; determine out cmap value
    ;---------------------------
      c1 = cmap1_dim-1  ; cmap1 stop index in ocmap
      c2 = ocamp_dim-1  ; cmap2 stop index in ocmap

      ocmap (0:c1,:)    = cmap1
      ocmap (c1+1:c2,:) = cmap2

;*** don't combine cmap2 foreground/background color ***

  else

    ;---------------------
    ; create output cmap
    ;---------------------
      ocamp_dim = cmap1_dim + cmap2_dim - 2
      ocmap = new((/ocamp_dim,3/),"float")

    ;---------------------------
    ; determine out cmap value
    ;---------------------------
      c1 = cmap1_dim-1  ; cmap1 stop index in ocmap
      c2 = ocamp_dim-1  ; cmap2 stop index in ocmap

      ocmap (0:c1,:)    = cmap1
      ocmap (c1+1:c2,:) = cmap2(2:,:)

  end if

;----------
; return
;----------
  return(ocmap)

end

;### end of ignorant_combine_colormaps ###

;###############################################
undef("ignorant_xy")
function ignorant_xy(wks:graphic, xx:numeric, yy:numeric , \
                     labels:string, nclres:logical, opt:logical)
begin

  xxx = xx  ; plot X variable
  yyy = yy  ; plot Y variable

;--------------------
; set labels & plot
;--------------------
  if ( isatt(opt, "MeanFormat")) then
    mean_format = opt@MeanFormat
  else 
    mean_format = "%7.3f"
  end if

  if (isatt(opt,"legend_mean") .and. opt@legend_mean.eq.False) then
    labels_work = "  "+labels

  else
    yyy_avg = dim_avg_Wrap(yyy)
    nlabels = dimsizes(labels)    ; avoid nlabels .ne. nyyy, e.g. ensemble plots
    nyyy    = dimsizes(yyy_avg)

    if (nlabels.eq.nyyy) then
      labels_work = "  "+labels + ", mean: "+sprintf(mean_format,yyy_avg)
    else
      labels_work = "  "+labels
    end if

  end if

;------------------------
; set default resources
;------------------------

  res_default = True
  res_default@gsnMaximize = True

;*******************
; set XY resources
;*******************

;*** set XY line mode ***
  res_default@xyMonoMarkLineMode = True
  res_default@xyMarkLineMode     = "Lines"

;*** set XY line color ***
  if (isatt(opt,"xy_colors1")) then
    xy_colors1 = opt@xy_colors1
  else
    xy_colors1 = (/ "red", \               
                    "blue",  \
                    "orange", \
  ;                  "yellow",  \
                    "green", \
                    "purple",  \
                    "pink",  \
                    "brown",  \
                    "gray",  \
                    "purple",  \
                    "black"  \
                  /)
  end if

  res_default@xyMonoLineColor    = False
  res_default@xyLineColors       = xy_colors1

;*** set XY line dash ***
  res_default@xyMonoDashPattern  = True
  res_default@xyDashPattern      = 0

;*** set XY line thickness ***
  res_default@xyMonoLineThickness = True
  res_default@xyLineThicknessF    = 2.5

;***********************
; set Legend resources
;***********************

;*** set Legend size
  res_default@pmLegendDisplayMode    = "Always"
  res_default@pmLegendWidthF      = 0.2
  res_default@pmLegendHeightF     = 0.1

;*** set Legend font ***
  res_default@lgAutoManage = False
  res_default@lgLabelFontHeightF = 0.02
  res_default@xyExplicitLegendLabels = labels_work

;*** set Legend perim
  res_default@lgPerimOn              = True

;**********************
; set Title resources
;**********************

;*** XY title default value ***
  x_title = " "
  y_title = " "

;*** set X-Axis title ***
  if ( isatt(xxx,"long_name") ) then
    x_title = xxx@long_name
  end if

  if ( isatt(xxx,"units") ) then
    x_title = x_title + " (" + xxx@units + ")"
  end if

;*** set Y-Axis title ***
  if ( isatt(yyy,"long_name") ) then
    y_title = yyy@long_name
  end if

  if ( isatt(yyy,"units") ) then
    y_title = y_title + " (" + yyy@units + ")"
  end if

;*** call title resources ***
  res_default@tiXAxisString = x_title
  res_default@tiYAxisString = y_title

;-------------------------------------------
; set work resource & plot
; replace default by input NCL resources
;-------------------------------------------

  res_work = res_default
  copyatt(res_work,nclres)

  if (isatt(opt,"marker")) then
    res_work@gsnFrame = False
  end if 
  plot = gsn_csm_xy(wks,xxx,yyy,res_work)

;-----------------------------
; draw opt@marker = p_values
;-----------------------------

  if (isatt(opt,"marker") .and. opt@marker .eq. "p_values" .and. isatt(opt,"prob")) then 

    ; read p_values
    prob = opt@prob
    nn = dimsizes(dimsizes(prob))
    if (isatt(opt,"siglvl")) then
      siglvl = opt@siglvl
    else
      siglvl = 0.05
    end if

    ; set x coordiante of marker
    if (isatt(opt,"marker_x")) then
      marker_x = xxx
      marker_x = opt@marker_x
    else
      marker_x = xxx
    end if

    ; set y coordiante of marker
    if (isatt(opt,"marker_y")) then
      marker_y = yyy
      marker_y = opt@marker_y
    else
      marker_y = yyy
    end if

    ; set marker res
    polyres = True
    polyres@gsMarkerIndex = 3
    polyres@gsMarkerSizeF = 15.
    copyatt(polyres,opt)

    ;*** if prob is a 1D array ***
    if (nn.eq.1) then
      ii = dimsizes(prob)
      do i=0,ii-1
        if (prob(i) .ge. 0. .and. prob(i) .le. siglvl) then
          x1 = marker_x(i)
          y1 = marker_y(i)
          gsn_polymarker(wks,plot,x1,y1,polyres)
        end if
      end do
    end if

    ;*** if prob is a 2D array ***
    if (nn.eq.2) then
      qq = dimsizes(prob)
      dum = new(qq, "graphic")
      ii = qq(0)
      jj = qq(1)
      do i=0,ii-1
      do j=0,jj-1
        if (prob(i,j) .ge. 0. .and. prob(i,j) .le. siglvl) then

          if (dimsizes(dimsizes(xxx)) .eq. 1) then
            x1 = marker_x(j)
          else if (dimsizes(dimsizes(xxx)) .eq. 2) then
            x1 = marker_x(i,j)
          end if
          end if

          if (dimsizes(dimsizes(yyy)) .eq. 1) then
            y1 = marker_y(j)
          else if (dimsizes(dimsizes(yyy)) .eq. 2) then
            y1 = marker_y(i,j)
          end if
          end if

          polyres@gsMarkerColor = xy_colors1(i)
          gsn_polymarker(wks,plot,x1,y1,polyres)
        end if
      end do
      end do
    end if

    ;*** frame plot *** 
    if (isatt(nclres,"gsnFrame") .and. .not.(nclres@gsnFrame)) then
      aa=0
    else
      bb=0
      frame(wks)
    end if

  end if  ; end if of opr@marker = "p_values"

;---------
; return
;---------
  return(plot)
  
end

;### end of ignorant_xy ###

;###############################################
undef("ignorant_xy_ensemble")
function ignorant_xy_ensemble (wks:graphic, xx:numeric, yy:numeric , \
                               labels:string, nclres:logical, opt:logical)
  ;opt = True
  ;opt@return_type = "mxmnag"  ; (return group,max/min/avg,npoints)
  ;opt@option_group = True     ; if grouping the first dimension
  ;opt@nmember_per_group = 3   ; number of member per group
  ;opt@XYStyle = "fill", "solid", "dash" (default is solid)
  ;nclres@gsnFrame = True      ; set gsnFrame to True if no panel plots

begin

;--------------------------------
; set line colors & fill colors
;--------------------------------
  colors_main = (/"red", \
                  "blue", \
                  "goldenrod", \
                  "green", \
                  "purple", \
                  "black" /)
                 
  colors_light = (/"lightpink", \
                  "lightblue", \
                  "lightgoldenrod",  \
                  "darkolivegreen1", \
                  "plum", \
                  "grey" /) 

;------------------
; create ensemble
;------------------
  opt_ens = opt

  xdim  = dimsizes(dimsizes(xx))
  ydim  = dimsizes(dimsizes(yy))
  if (.not.isatt(opt_ens,"ensemble") .and. xdim.eq.2 .and. ydim.eq.1) then
    opt_ens@ensemble = "x"

  else if (.not.isatt(opt_ens,"ensemble") .and. xdim.eq.1 .and. ydim.eq.2) then
    opt_ens@ensemble = "y"

  else
    opt_ens@ensemble = "error"
    print("ERROR (): either xx or yy MUST be a two-dimension array ")
    exit
  end if
  end if

  ;*** don't use mxmnag because it will produce wrong area, e.g. X -> O ***
  ;if (isatt(opt,"XYFill") .and. opt@XYFill) then
  ;  opt_ens@return_type = "mxmnag"
  ;end if
  ;print(opt_ens)

  if (isatt(opt_ens,"ensemble") .and. opt_ens@ensemble .eq. "x") then
    xx_ens = ensemble_compute(xx, opt_ens)
    yy_ens = yy
    ngroup = dimsizes(xx_ens(:,0,0))
    nvalue = dimsizes(xx_ens(0,:,0))
  else if (isatt(opt_ens,"ensemble") .and. opt_ens@ensemble .eq. "y") then
    xx_ens = xx
    yy_ens = ensemble_compute(yy, opt_ens)
    ngroup = dimsizes(yy_ens(:,0,0))
    nvalue = dimsizes(yy_ens(0,:,0))

  else
    xx_ens = 0.
    yy_ens = 0.
    ngroup = 0.
    nvalue = 0.
  end if
  end if

  ;printVarSummary(xx_ens)
  ;printVarSummary(yy_ens)
  ;print(xx_ens)

;---------
; set res
;---------
  res_work = True
  ;res_work@gsnFrame = False
  ;res_work@gsnDraw = False
  
  res_work@pmLegendDisplayMode = "noCreate"
  res_work@xyMonoLineThickness = True
  res_work@xyLineThicknessF = 3.

  res_work@xyLineColors = colors_main
  res_work@gsnXYFillColors = colors_light
  res_work@gsnXYFillOpacities = 0.15
  
  res_work@trXMaxF = max(xx_ens)
  res_work@trXMinF = min(xx_ens)
  res_work@trYMaxF = max(yy_ens)
  res_work@trYMinF = min(yy_ens)
  
  copy_VarAtts(nclres, res_work)

  ;print(res_work)

  res_work@gsnFrame = False
  res_work@gsnDraw = False
  colors_line = res_work@xyLineColors
  colors_fill = res_work@gsnXYFillColors

  colors1_lines  = new(nvalue, "string")
  colors1_dashes = new(nvalue, "integer")

;---------
; plot
;---------
  plot_all = new(ngroup, "graphic")

  labels_work = new(ngroup, "string")
  labels_work = labels

  do i=0,ngroup-1
    color1_line = colors_line(i) 
    color1_fill = colors_fill(i) 

    delete(res_work@xyLineColors)
    delete(res_work@gsnXYFillColors)

    if (isatt(opt,"XYStyle") .and. opt@XYStyle.eq."fill") then
      colors1_lines = "transparent"
      colors1_lines(nvalue-1) = color1_line  ; ensemble mean

      res_work@gsnXYFillColors = color1_fill
      res_work@xyLineColors = colors1_lines 
      ;res_work@xyLineColors = (/"transparent","transparent",color1_line/)

    else if (isatt(opt,"XYStyle") .and. opt@XYStyle.eq."dash") then
      colors1_lines   = color1_line
      colors1_dashes = 16 ; dash
      colors1_dashes(nvalue-1) = 0 ; solid line for ensemble mean

      res_work@gsnXYFillColors = "transparent"
      res_work@xyLineColors = colors1_lines
      res_work@xyMonoDashPattern = False
      res_work@xyDashPatterns = colors1_dashes

    else
      colors1_lines = color1_fill
      colors1_lines(nvalue-1) = color1_line
      res_work@gsnXYFillColors = "transparent"
      res_work@xyLineColors = colors1_lines

      colors1_opac = new( dimsizes(colors1_lines), "float")
      colors1_opac = 0.7     
      colors1_opac(nvalue-1) = 1.
      res_work@xyLineOpacities = colors1_opac
      delete(colors1_opac)
    end if
    end if

    if (isatt(opt_ens,"ensemble") .and. opt_ens@ensemble .eq. "x") then
      plot_all(i) =  ignorant_xy(wks, xx_ens(i,:,:), yy_ens, labels_work, res_work, opt)
    else if (isatt(opt_ens,"ensemble") .and. opt_ens@ensemble .eq. "y") then
      plot_all(i) =  ignorant_xy(wks, xx_ens, yy_ens(i,:,:), labels_work, res_work, opt)
    end if
    end if

  end do

  do i=1,ngroup-1
    overlay(plot_all(0),plot_all(i))
  end do

  plot = plot_all(0)

  if (isatt(nclres, "gsnFrame") .and. nclres@gsnFrame ) then
    draw(plot)
    frame(wks)
  end if  

  return(plot)

end
;### end of ignorant_xy_ensemble ###

;###################################################
undef("ignorant_merge_colormaps")
function ignorant_merge_colormaps(icmap)

begin

;------------------------------
; define work color map array
;------------------------------
  wcmap_size = 10000
  wcmap_val = new( (/wcmap_size,3/), "float" )
  wcmap_val = -999.99

;*********************************
; merge NCL predefined color map
;*********************************

    if ( isstring(icmap) ) then

      icamp_number = dimsizes(icmap)

      do j = 0, icamp_number-1

;--------------------------------
; open wks & get colormap value
;--------------------------------
        wks = gsn_open_wks("ncgm","owrsfnmks")
        gsn_define_colormap( wks,icmap(j) )

        getvalues wks
          "wkColorMap" : icmap_j_val
        end getvalues

        icmap_j_size = dimsizes(icmap_j_val(:,0))

;---------------------------------
; decide work colormap position
; to put read-in color map value
;---------------------------------
        wcmap_j_index_start = 0
        do k = 0, wcmap_size-1
          if ( wcmap_val(k,0) .ne. -999.99 ) then
            wcmap_j_index_start = k+1
          end if
        end do

        wcmap_j_index_end = wcmap_j_index_start + icmap_j_size - 1

        wstart = wcmap_j_index_start
        wend   = wcmap_j_index_end

        wcmap_val(wstart:wend ,:) = icmap_j_val

;----------------------
; delete used varbles
;----------------------
        delete(icmap_j_val)
        delete(icmap_j_size)
        delete(wcmap_j_index_start)
        delete(wcmap_j_index_end)
        delete(wstart)
        delete(wend)       

      end do

;----------------------------
; return new color map value
;----------------------------

; remove temporary wks
      system("rm -f owrsfnmks.ncgm*")

; find output colormap size
      do k = 0, wcmap_size-1
        if ( wcmap_val(k,0) .ne. -999.99 ) then
          ocmap_size = k+1
        end if
      end do

; determine output colormap value
      ocmap_val = new( (/ocmap_size,3/),"float" )

      do k = 0, ocmap_size-1
        ocmap_val(k,:) = wcmap_val(k,:)
      end do

      return(ocmap_val)

  end if

end

;### end of ignorant_merge_colormaps ###

;#####################################################3
undef("ignorant_select_colormap")
function ignorant_select_colormap(icmap_name:string, ocmap_index[*]:numeric, \
                                  option_foreback:logical)

begin

;-------------------------------
; get input colormap RGB array
;-------------------------------

  icmap_val0 = read_colormap_file(icmap_name) ; not read fore/back ground colors

  ndim = dimsizes(icmap_val0(:,0))   ; add fore/background colors
  ndim2 = ndim+2
  ndim1 = ndim2-1
  icmap_val = new( (/ndim+2,3/),typeof(icmap_val0)) ; only keep RGB values

  icmap_val(0,:) = 1.
  icmap_val(1,:) = 0.
  icmap_val(2:ndim1,:) = icmap_val0(:,0:2)

; if user select all colors in given colormap
  if ( dimsizes(ocmap_index) .eq. 1 .and. ocmap_index .eq. 999 ) then

    if (option_foreback) then
      ocmap_val = icmap_val
    else
      ocmap_val = icmap_val(2:,:)
    end if

    return(ocmap_val)
  end if

  icmap_size = dimsizes(icmap_val)

;---------------------------
; write out output colormap
;---------------------------

; add fore and background colors (white & black)
  if (option_foreback) then

    ocmap_size = dimsizes(ocmap_index) + 2
    ocmap_val = new((/ocmap_size,3/),"float")

    ocmap_val(0,:) = (/1., 1., 1./) ; set foreground color to black
    ocmap_val(1,:) = (/0., 0., 0./) ; set background color to white

    do i=2,ocmap_size-1
      j = ocmap_index(i-2)
      ocmap_val(i,:) = icmap_val(j,0:2)
    end do

; not add fore and background colors
  else

    ocmap_size = dimsizes(ocmap_index)
    ocmap_val = new((/ocmap_size,3/),"float")

    do i=0,ocmap_size-1
      j = ocmap_index(i)
      ocmap_val(i,:) = icmap_val(j,0:2)
    end do

  end if

;-----------------------------------
; return output colormap RGB array
;----------------------------------
  return(ocmap_val)

end
;
;
;------------------------------
undef("ignorant_write_nc")
procedure ignorant_write_nc(filename:string,outvarname:string,invar,opt:integer)

begin

  outfile = addfile(filename,"w")

;*******************************************************
;  outfile have unlimited time, overwrite from max time
;*******************************************************
  if (opt.eq.-1) then

    outvar = outfile->$outvarname$
    dim_outvar = dimsizes(outvar)
    dimT_out = dim_outvar(0)
    nout = dimT_out

    outfile->$outvarname$(nout,:,:) = (/invar/)

  end if

;*********************************************************
;  outfile have unlimited time, overwrite from first time
;*********************************************************
  if (opt.eq.-2) then

    nout = 0
    outfile->$outvarname$(nout,:,:) = (/invar/)

  end if

;**************************************************************
;  outfile don't have [outvarname] varible, 
;  automatically create [outvarname] varible and dimensions
;**************************************************************
  if (opt.eq.-3) then

    outfile->$outvarname$ = (/invar/)

  end if

;*********************************************************
;  Don't write in unlimited time varible
;*********************************************************
  if (opt.ge.0) then

    outfile->$outvarname$(opt,:,:) = (/invar/)

  end if

 
end

;---------------------------
undef("ignorant_LWP")
function ignorant_LWP(qc[*][*][*],p[*][*][*])

begin

;*** know input data dimesion ***
;
  dimqc = dimsizes(qc)
  dimz = dimqc(0)
  dimy = dimqc(1)
  dimx = dimqc(2)
;
;*** calculate liquid water path ***
;
  qcmean = new((/dimz-1,dimy,dimx/),float)
  dp = new((/dimz-1,dimy,dimx/),float)

  LWP = new((/dimy,dimx/),float)
  LWP!0 = qc!1
  LWP!1 = qc!2
  LWP@description = "liquid water path"
  LWP@units = "kg m-2"

  do k=0,dimz-1-1
    qcmean(k,:,:) = (qc(k+1,:,:) + qc(k,:,:))/2.
     dp(k,:,:) = p(k,:,:) - p(k+1,:,:)
  end do


  do i=0,dimx-1
    do j=0,dimy-1

      LWP(j,i) =  prcwater_dp(qcmean(:,j,i),dp(:,j,i))

    end do
  end do

  return(LWP)

end

;--------------------------------
undef("ignorant_HMLcloud")
function ignorant_HMLcloud(P[*][*][*]:numeric,Tk[*][*][*]:numeric,qv[*][*][*]:numeric,qall[*][*][*]:numeric)

begin

  dimP = dimsizes(P)
  nz = dimP(0)
  ny = dimP(1)
  nx = dimP(2)

  HML = new((/3,ny,nx/),float)
  HML = -100.
  HML@description = "High,Middle,Low cloud amount"
  HML@HighCloud = "index 2, 300-550 mb"
  HML@MiddleCloud = "index 1, 550-800 mb"
  HML@LowCloud = "index 0, 800-1000 mb"

  Pqs = 6.11*exp(19.83-5417./Tk)  ; Compute saturation water vapor pressure in millibar
  qsat = (Pqs/P)*18./28.56        ; Change Pqs to saturation water vapor mixing ratio

  do j=0,ny-1
    do i=0,nx-1
      do k=0,nz-1

        if (P(k,j,i).gt.800..and.P(k,j,i).le.1000.) then       ; Calaulate High cloud amount
          ratio0 = qall(k,j,i)/(0.01*qsat(k,j,i))
          ratio0 = ratio0 < 1.

          if (ratio0.gt.HML(0,j,i)) then
            HML(0,j,i) = ratio0
          end if
        end if

        if (P(k,j,i).gt.550..and.P(k,j,i).le.800.) then        ; Calculate Middle cloud amount
          ratio1 = qall(k,j,i)/(0.01*qsat(k,j,i))
          ratio1 = ratio1 < 1.

          if (ratio1.gt.HML(1,j,i)) then
            HML(1,j,i) = ratio1
          end if
        end if

        if (P(k,j,i).gt.300..and.P(k,j,i).le.550.) then        ; Calculate Low cloud amount
          ratio2 = qall(k,j,i)/(0.01*qsat(k,j,i))
          ratio2 = ratio2 < 1.

          if (ratio2.gt.HML(2,j,i)) then
            HML(2,j,i) = ratio2
          end if
        end if

      end do
    end do
  end do

  return(HML)

end

;-------------------------------------------------
undef("ignorant_contour")
function ignorant_contour(wks[1]:graphic,data[*][*]:numeric,opt[1]:logical)

begin

  res = opt

  res@gsnMaximize = True
  
;***********************
;  Contour Attribute
;***********************

  cn_logical = (/"cnFillOn", "cnMonoFillPattern", "cnMonoFillColor", \
                 "cnInfoLabelOn","cnLineLabelsOn"/)

  cnatt_logical =(/True,True,False,False,True/)

  cndim = dimsizes(cn_logical)

  do i=0,cndim-1
    cc =isatt(opt,cn_logical(i))
    if (cc.eq.False) then
      res@$cn_logical(i)$ =  cnatt_logical(i)
    end if
  end do

;************************************************************
;  Label Bar Attribute  (only valid when cnFillOn = True)
;************************************************************

if (isatt(opt,"cnFillOn").eq.False.or.opt@cnFillOn.eq.True) then

   res@cnLineLabelsOn = False
   res@cnLinesOn = False

   lb_string = (/"pmLabelBarDisplayMode","pmLabelBarSide","lbOrientation"/)
   lbatt_string = (/"Always","Bottom","Horizontal"/)

   lb_num = (/"pmLabelBarWidthF","pmLabelBarHeightF"/)
   lbatt_num = (/0.8,0.1/)

   lb_logical = (/"lbPerimOn"/)
   lbatt_logical = (/False/)

   dimlbstr = dimsizes(lb_string)
   dimlbnum = dimsizes(lb_num)
   dimlblog = dimsizes(lb_logical)

  do i=0,dimlbstr-1
    cc =isatt(opt,lb_string(i))
    if (cc.eq.False) then
      res@$lb_string(i)$ =  lbatt_string(i)
    end if
  end do

  do i=0,dimlbnum-1
    cc =isatt(opt,lb_num(i))
    if (cc.eq.False) then
      res@$lb_num(i)$ =  lbatt_num(i)
    end if
  end do

  do i=0,dimlblog-1
    cc =isatt(opt,lb_logical(i))
    if (cc.eq.False) then
      res@$lb_logical(i)$ =  lbatt_logical(i)
    end if
  end do

end if


  plot = gsn_contour(wks,data,res)

  return(plot)

end
;-------------------------------------------------------

undef("ignorant_nc_contour")
function ignorant_nc_contour( wks[1]:graphic, nc_file:file, data[*][*]:numeric, opt:logical, nclres:logical )

begin

  res = nclres
  res@tfDoNDCOverlay = True

  if (isatt(res,"gsnMaximize").eq.False) then
    res@gsnMaximize = True
  end if

;***********************************************
;  Determine X,Y range , default is full range
;***********************************************

  pos = (/"XLeft","XRight","YBottom","YTop"/)
  dim = dimsizes(data)
  pos_val = (/0,dim(1)-1,0,dim(0)-1/)
  do i=0,3
    if (isatt(opt,pos(i)).eq.True) then
      pos_val(i) = opt@$pos(i)$
    end if
  end do

;********************************************************
;  Determine Map Corners Latitude and Longitude value
;********************************************************

    proj_att = (/"MAP_PROJ", "TRUELAT1", "TRUELAT2", "STAND_LON",  "REF_LON", \
                 "REF_LAT" ,   "KNOWNI",   "KNOWNJ",  "POLE_LON", "POLE_LAT", \
                       "DX",       "DY",   "LATINC",    "LONINC" /)

    dim_projatt = dimsizes(proj_att)
 
    projres = True

    do i=0,dim_projatt-1

      if ( isatt(opt,proj_att(i)).eq.True) then
        projres@$proj_att(i)$ = opt@$proj_att(i)$

      else

        if (isatt(nc_file,proj_att(i)).eq.True) then
          projres@$proj_att(i)$ = nc_file@$proj_att(i)$
        end if

      end if

    end do


    if (isatt(opt,"WRFfile").eq.True.and.opt@WRFfile.eq.True) then

       LeftLon  = nc_file->XLONG(0,pos_val(2),pos_val(0))
       LeftLat  = nc_file->XLAT(0,pos_val(2),pos_val(0))
       RightLon = nc_file->XLONG(0,pos_val(3),pos_val(1))
       RightLat = nc_file->XLAT(0,pos_val(3),pos_val(1))

    else

       llLeft  = wrf_ij_to_ll(pos_val(0),pos_val(2),projres)
       llRight = wrf_ij_to_ll(pos_val(1),pos_val(3),projres)

       LeftLon  = llLeft(0)
       LeftLat  = llLeft(1)
       RightLon = llRight(0)
       RightLat = llRight(1) 
  
    end if

;************************************
;  Map Attribute
;************************************

   if (projres@MAP_PROJ.eq.1) then

     res@mpProjection = "LambertConformal"
     res@mpLambertMeridianF = projres@STAND_LON
     res@mpLambertParallel1F = projres@TRUELAT1
     res@mpLambertParallel2F = projres@TRUELAT2

   end if

   res@mpLimitMode = "Corners"
   res@mpLeftCornerLonF =  LeftLon
   res@mpLeftCornerLatF =  LeftLat
   res@mpRightCornerLonF = RightLon
   res@mpRightCornerLatF = RightLat


   if (isatt(res,"mpPerimOn").eq.False) then
     res@mpPerimOn = True
   end if

   if (isatt(res,"mpGridAndLimbOn").eq.False) then
     res@mpGridAndLimbOn = False
   end if

   if (isatt(res,"mpDataBaseVersion").eq.False) then
     res@mpDataBaseVersion = "Ncarg4_1"
   end if

   if (isatt(res,"mpDataSetName").eq.False) then
     res@mpDataSetName = "Earth..4" 
   end if

;***********************
;  Contour Attribute
;***********************

  cn_logical = (/"cnFillOn", "cnMonoFillPattern", "cnMonoFillColor", \
                 "cnInfoLabelOn","cnLineLabelsOn"/)

  cnatt_logical =(/True,True,False,False,True/)

  cndim = dimsizes(cn_logical)

  do i=0,cndim-1
    cc =isatt(nclres,cn_logical(i))
    if (cc.eq.False) then
      res@$cn_logical(i)$ =  cnatt_logical(i)
    end if
  end do

;************************************************************
;  Label Bar Attribute  (only valid when cnFillOn = True)
;************************************************************

if (isatt(nclres,"cnFillOn").eq.False.or.nclres@cnFillOn.eq.True) then

   if (isatt(nclres,"cnLineLabelsOn").eq.False) then
     res@cnLineLabelsOn = False
   end if

   if (isatt(nclres,"cnLinesOn").eq.False) then
     res@cnLinesOn = False
   end if

   lb_string = (/"pmLabelBarDisplayMode","pmLabelBarSide","lbOrientation"/)
   lbatt_string = (/"Always","Bottom","Horizontal"/)

   lb_num = (/"pmLabelBarWidthF","pmLabelBarHeightF"/)
   lbatt_num = (/0.8,0.1/)

   lb_logical = (/"lbPerimOn"/)
   lbatt_logical = (/False/)

   dimlbstr = dimsizes(lb_string)
   dimlbnum = dimsizes(lb_num)
   dimlblog = dimsizes(lb_logical)

  do i=0,dimlbstr-1
    cc =isatt(nclres,lb_string(i))
    if (cc.eq.False) then
      res@$lb_string(i)$ =  lbatt_string(i)
    end if
  end do

  do i=0,dimlbnum-1
    cc =isatt(nclres,lb_num(i))
    if (cc.eq.False) then
      res@$lb_num(i)$ =  lbatt_num(i)
    end if
  end do

  do i=0,dimlblog-1
    cc =isatt(nclres,lb_logical(i))
    if (cc.eq.False) then
      res@$lb_logical(i)$ =  lbatt_logical(i)
    end if
  end do

end if

; Label Bar Title

  name = " "
  unit = " "

  if (isatt(res,"lbTitleString").eq.False) then
    
    if (isatt(data,"description").eq.True) then
      name = data@description
    end if

    if (isatt(data,"units").eq.True) then
      unit = data@units
    end if

    res@lbTitleString = (name + "    (" + unit + ")")

  end if

;*********************************
;  Tick Mark Attribute
;*********************************

   if (isatt(res,"pmTickMarkDisplayMode").eq.False) then  
     res@pmTickMarkDisplayMode = "Conditional"
   end if

   if (isatt(res,"tmXTLabelsOn").eq.False) then
     res@tmXTLabelsOn = False
   end if

   if (isatt(res,"tmXTOn").eq.False) then
     res@tmXTOn = False
   end if

   if (isatt(res,"tmYROn").eq.False) then
     res@tmYROn = False
   end if

   if (isatt(res,"tmXBLabelFontHeightF").eq.False) then
     res@tmXBLabelFontHeightF = 0.005
   end if

;********************************
;  Plot and return
;********************************

   plot = gsn_contour_map(wks,data,res)

   return(plot)

end

;### end of ignorant_select_colormap ###

;#############################################
undef("change_var_atts")
procedure change_var_atts (varname[1]:string, varin, res[1]:logical, opt_chgatt[1]:logical)
begin

if (opt_chgatt) then

;*************
; initialize
;*************

  if ( isatt(varin,"long_name") ) then
    longname = varin@long_name
  else
    longname = "slkdjflksjdfljasldfjlaks"
  end if

;***************
; change units
;***************
  if ( \
           .not.ismissing(str_match(longname,"precipitation rate")) \
      .or. .not.ismissing(str_match(longname,"Precipitation rate")) \
     ) then
      if ( isatt(varin, "units") .and. varin@units.eq."m/s") then
        varin  = unit_convert(varin,"m/s","mm/day")
      end if
  end if

  if ( \
           .not.ismissing(str_match(longname,"Specific humidity")) \
      .or. .not.ismissing(str_match(longname,"specific humidity")) \
     ) then
    if ( isatt(varin, "units") .and. varin@units.eq."kg/kg") then
      varin  = unit_convert(varin,"kg/kg","g/kg")
    end if
  end if

  if ( \
           .not.ismissing(str_match(longname,"Precipitable water")) \
      .or. .not.ismissing(str_match(longname,"precipitable water")) \
     ) then
    if ( isatt(varin, "units") .and. varin@units.eq."kg/m2") then
      varin  = unit_convert(varin,"kg/m2","mm")
    end if
  end if

  if (varname.eq."TGCLDLWP" .or. varname.eq."TGCLDIWP") then
    if ( isatt(varin, "units") .and. varin@units.eq."kg/m2") then
      varin  = unit_convert(varin,"kg/m2","g/m2")
    end if
  end if

  if (varname.eq."CLDLIQ" .or. varname.eq."CLDICE") then
    if ( isatt(varin, "units") .and. varin@units.eq."kg/kg") then
      varin  = unit_convert(varin,"kg/kg","mg/kg")
    end if
  end if

  if (varname.eq."OMEGA") then
    if ( isatt(varin, "units") .and. varin@units.eq."Pa/s") then
      varin  = unit_convert(varin,"Pa/s","hPa/day")
    end if
  end if

  if (varname.eq."FSNS") then
    varin0 = varin
    ;varin = -1. * varin
    ;print(varin0+" , "+varin)
  end if

  if ( isatt(varin, "units") .and. varin@units.eq."K/s") then
    varin  = unit_convert(varin,"K/s","K/day")
  end if

  if ( isatt(varin, "units") .and. varin@units.eq."fraction") then
    varin  = unit_convert(varin,"fraction","%")
  end if

;*****************************************
; change fill colors, for absolute value
;*****************************************

  if (isatt(opt_chgatt,"chg_cncolors") .and. opt_chgatt@chg_cncolors) then

    opt_chgcolors = False
    option_cmap = "default"

    if ( \
             .not.ismissing(str_match(longname,"Vertical velocity")) \
        .or. .not.ismissing(str_match(longname,"OMEGA")) \
       ) then
       option_cmap = "red_blue"
       opt_chgcolors = True
    end if

    if ( \
             .not.ismissing(str_match(longname,"heating rate")) \
        .or. .not.ismissing(str_match(longname,"Heating rate")) \
       ) then
       option_cmap = "cold_warm"
       opt_chgcolors = True
    end if

    if ( \
             .not.ismissing(str_match(longname,"warming rate")) \
        .or. .not.ismissing(str_match(longname,"Warming rate")) \
       ) then
       option_cmap = "cold_warm"
       opt_chgcolors = True
    end if

    if ( \
             .not.ismissing(str_match(longname,"cooling rate")) \
        .or. .not.ismissing(str_match(longname,"Cooling rate")) \
       ) then
       option_cmap = "warm_cold"
       opt_chgcolors = True
    end if
  
    if ( \
             .not.ismissing(str_match(longname,"Zonal wind")) \
        .or. .not.ismissing(str_match(longname,"zonal wind")) \
        .or. .not.ismissing(str_match(longname,"Meridional")) \
        .or. .not.ismissing(str_match(longname,"Meridional wind")) \
        .or. .not.ismissing(str_match(longname,"meridional wind")) \
        .or. .not.ismissing(str_match(longname,"U at")) \
        .or. .not.ismissing(str_match(longname,"V at")) \
        .or. .not.ismissing(str_match(longname,"meridional stream function")) \
        .or. .not.ismissing(str_match(longname,"Meridional stream function")) \
        .or. .not.ismissing(str_match(longname,"Meridional Stream Function")) \
       ) then
       option_cmap = "purple_green"
       opt_chgcolors = True
    end if
  
    if (opt_chgcolors) then

      if (isatt(res,"cnFillPalette")) then
        delete(res@cnFillPalette)
      end if 
  
      opt_work = True
      if (isatt(res,"max_levels")) then
        opt_work@max_levels = res@max_levels  ; level for using function nice_mnmxintvl
      else
        opt_work@max_levels = 15  ; level for using function nice_mnmxintvl
      end if  

      if (isatt(res,"cnlevels_predefine") .and. res@cnlevels_predefine) then
        opt_work@cnlevels_predefine = True    ; True  - contour levels are given in res
      else
        opt_work@cnlevels_predefine = False    ; False - contour levels will be computed in this procedure
      end if
 
  if (option_cmap .eq. "dry_wet") then
    opt_work@colormap_positive = "MPL_YlGnBu"
    opt_work@colormap_positive_index_start = 10      ; the least positive color
    opt_work@colormap_positive_index_end = 113       ; the most  positive color

    opt_work@colormap_negative = "MPL_pink"
    opt_work@colormap_negative_index_start = 17      ; the most  negative color
    opt_work@colormap_negative_index_end = 129       ; the least negative color
  end if

 
;      if (option_cmap .eq. "dry_wet") then
;        opt_work@colormap_positive = "MPL_YlGnBu"
;        ;opt_work@colormap_positive_index_start = 2       ; the least positive color
;        ;opt_work@colormap_positive_index_end = 999       ; the most  positive color
;  
;        ;opt_work@colormap_negative = "MPL_YlOrBr"
;        ;opt_work@colormap_negative_index_start = 120     ; the most  negative color
;        ;opt_work@colormap_negative_index_end = 2         ; the least negative color
;        opt_work@colormap_negative = "MPL_YlGnBu"
;        opt_work@colormap_negative_index_start = 2     ; the most  negative color
;        opt_work@colormap_negative_index_end = 129         ; the least negative color
;      end if
  
      if (option_cmap .eq. "purple_green") then
        opt_work@colormap_positive = "MPL_Purples"       ; colormap for plotting positive colors
        opt_work@colormap_positive_index_start = 2       ; the least positive color
        opt_work@colormap_positive_index_end = 999       ; the most  positive color
    
        opt_work@colormap_negative = "MPL_Greens"        ; colormap for plotting negative colors
        opt_work@colormap_negative_index_start = 120     ; the most  negative color
        opt_work@colormap_negative_index_end = 2         ; the least negative color
      end if
    
      if (option_cmap .eq. "red_blue") then
        opt_work@colormap_positive = "NCV_blu_red"         ; colormap for plotting positive colors
        opt_work@colormap_positive_index_start = 123       ; the least positive color
        opt_work@colormap_positive_index_end = 2           ; the most  positive color
    
        opt_work@colormap_negative = "NCV_blu_red"         ; colormap for plotting negative colors
        opt_work@colormap_negative_index_start = 255       ; the most  negative color
        opt_work@colormap_negative_index_end = 137         ; the least negative color
      end if

      if (option_cmap .eq. "cold_warm") then
        opt_work@colormap_positive = "WhiteYellowOrangeRed"
        opt_work@colormap_positive_index_start = 17          ; the least positive color
        opt_work@colormap_positive_index_end = 230           ; the most  positive color
        opt_work@colormap_negative = "WhiteBlue"
        opt_work@colormap_negative_index_start = 230         ; the most  negative color
        opt_work@colormap_negative_index_end = 33            ; the least negative color
      end if

      if (option_cmap .eq. "warm_cold") then
        opt_work@colormap_positive = "WhiteBlue"
        opt_work@colormap_positive_index_start = 33           ; the least positive color
        opt_work@colormap_positive_index_end = 230            ; the most  positive color
        opt_work@colormap_negative = "WhiteYellowOrangeRed"
        opt_work@colormap_negative_index_start = 230          ; the least negative color
        opt_work@colormap_negative_index_end = 17             ; the most  negative color
      end if
    
      if (option_cmap .eq. "default") then
        opt_work@colormap_positive = "WhiteYellowOrangeRed"  ; colormap for plotting positive colors
        opt_work@colormap_positive_index_start = 2           ; the least positive color
        opt_work@colormap_positive_index_end = 999           ; the most  positive color
    
        opt_work@colormap_negative = "MPL_GnBu"              ; colormap for plotting negative colors
        opt_work@colormap_negative_index_start = 127         ; the most  negative color
        opt_work@colormap_negative_index_end = 2             ; the least negative color
      end if
    
      copy_VarAtts(opt_chgatt, opt_work)
  
      ndim = dimsizes(dimsizes(varin))
      if (ndim.eq.2) then
        varin_temp1 = varin
      else if (ndim.eq.3) then
        varin_temp1 = dim_avg_n_Wrap(varin,0)
      else if (ndim.eq.4) then
        varin_temp0 = dim_avg_n_Wrap(varin,0)
        varin_temp1 = dim_avg_n_Wrap(varin_temp0,0)
      end if
      end if
      end if
  
      cnColors_PosNeg(varin_temp1, res, opt_work)
      ;copy_VarAtts(opt_work, opt_chgatt)
    end if ; end if of opt_chgcolors
  
  end if ; end if of opt_chgatt@chg_cncolors


;***********************
; change units string
;***********************
;
  if (isatt(varin,"units") .and. varin@units .eq. "W/m2") then
    varin@units = "W m~S~-2~N~"
    res@gsnRightString = varin@units
  end if

  if (isatt(varin,"units") .and. varin@units .eq. "gram/m2") then
    varin@units = "g m~S~-2~N~"
    res@gsnRightString = varin@units
  end if

  if (isatt(varin,"units") .and. varin@units .eq. "g/m2") then
    varin@units = "g m~S~-2~N~"
    res@gsnRightString = varin@units
  end if

  if (isatt(varin,"units") .and. varin@units .eq. "kg/m2") then
    varin@units = "kg m~S~-2~N~"
    res@gsnRightString = varin@units
  end if

;*********************
; change fill colors
;*********************
  ;*** if turn on, cloud frac diff will also use this level ***
  ;if (varname.eq."CLDLOW" .or. varname.eq."CLDMED" .or. varname.eq."CLDHGH" .or. varname.eq."CLDTOT") then
  ;if (isatt(varin,"units") .and. varin@units .eq. "fraction") then
  ;  res@cnLevelSelectionMode = "ExplicitLevels"
  ;  res@cnLevels = (/0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
  ;end if

;*************************
; change variables string
;*************************
;
  if (isatt(varin,"long_name")) then
    res@gsnLeftString = varin@long_name
  end if

  var_strings = (/  \
                  (/"CLDLOW"		, "Low cloud fraction"/),  \ 
                  (/"CLDMED"		, "Middle cloud fraction"/),  \ 
                  (/"CLDHGH"		, "High cloud fraction"/),  \ 
                  (/"CLDTOT"		, "Total cloud fraction"/),  \ 
                  (/"EVAP"		, "Evaporation rate"/),  \ 
                  (/"FLDS"		, "Dn LW flux at the surf"/), \
                  (/"FLDSC"		, "clrsky surf Dn LW flux"/), \
                  (/"FLNS"		, "Net Up LW flux at the surface"/), \
                  (/"FLNSC"		, "clrsky Surf Net Up LW flux"/), \
                  (/"FLUS"		, "Up LW flux at the surf"/), \
                  (/"FLUT"		, "Outgoing LW radiation"/), \
                  (/"FLUTC"		, "clrsky OLR"/), \
                  (/"FSDS"		, "Dn SW flux at the surf"/), \
                  (/"FSDSC"		, "clr-sky Surf Dn SW flux"/), \
                  (/"FSNS"		, "Net Surf Dw SW flux"/), \
                  (/"FSNSC"		, "clrsky Net Surf Dw SW flux"/), \
                  (/"FSUS"		, "Up SW flux at the surf"/), \
                  (/"FSUSC"		, "clr-sky Surf Up SW flux"/), \
                  (/"FSNTOA"		, "Net Dw SW flux at the TOA"/), \
                  (/"FSNTOAC"		, "clrsky TOA Net Dw SW flux"/), \
                  (/"FSUTOA"            , "Up SW flux at the TOA"/), \
                  (/"FSUTOAC"           , "clrsky TOA Up SW flux"/), \
                  (/"ICEFRAC"           , "Sea ice fraction"/), \
                  (/"LHFLX"		, "Latent heat flux"/), \
                  (/"LWCF_SURF"         , "Surf LW Cloud Forcing"/), \
                  (/"LWCOOL_TOA_TO_SURF", "Column LW Cooling"/), \
                  (/"MET"		, "Meridional energy transport"/), \
                  (/"NETCOOL_TOA_TO_SURF","Column Net Cooling"/), \
                  (/"OMEGA500"		, "~F33~w~F~ at 500 hPa"/), \
                  (/"P_minus_E"         , "Precip - Evap"/), \
                  (/"PRECC"		, "Conv Prec"/), \
                  (/"PRECL"		, "Stratiform Prec"/), \
                  (/"PRECT"		, "Total Prec"/), \
                  (/"RADCOOL_TOA_TO_SURF","Column Net Rad Cooling"/), \
                  (/"RADSURF"		, "Surf net down energy flux"/), \
                  (/"SHFLX"		, "Sensible heat flux"/), \
                  (/"SWCF_SURF"         , "Surf SW Cloud Forcing"/), \
                  (/"SWWARM_TOA_TO_SURF", "Column SW Heating"/), \
                  (/"TMQ"		, "Precipitable water"/), \
                  (/"TS"		, "Surf Rad temperature"/), \
                  (/"TREFHT"		, "Surface air temperature"/), \
                  (/"TGCLDIWP"		, "Grid-box ice water path"/), \
                  (/"TGCLDLWP"		, "Grid-box liquid water path"/) \
                  ;(/""	, /), \
                /)
  nstring = dimsizes(var_strings(:,0))
  do i=0,nstring-1
    if (varname.eq.var_strings(i,0)) then
      res@gsnLeftString = var_strings(i,1)
      varin@long_name = var_strings(i,1)
    end if
  end do

end if ; end if of opt_chgatt

end

;### end of change_var_atts ###

;#############################################
undef("ignorant_diagnose_var")
function ignorant_diagnose_var (infilename[1]:string, varname[1]:string, opt[1]:logical)
begin

;-------------
; description
;-------------

  ; available varnames_diagnose, update 03/02/2018
  ;   EVAP                      : Evaporation rate, converted from latent heat flux
  ;   FLDSC                     : Clearsky Downwelling longwave flux at surface (assume blackbody surface)
  ;   FLUS                      : Upward LW flux at the surface
  ;   FSUS                      : Upward SW flux at the surface
  ;   FSDTOA                    : Downward SW flux at the TOA, i.e. SOLIN
  ;   LWCF_SURF                 : LW cloud forcing at the surface       
  ;   LWCOOL_MT_TO_SURF         : LW radiative cooling rate between model top and surface, +cooling, -warming
  ;   LWCOOL_TOA_TO_SURF        : LW radiative cooling rate between the TOA and surface, +cooling, -warming
  ;   OMEGA at a pressure level : OMEGA at a certain pressure level, e.g. OMEGA500
  ;   NETCOOL_TOA_TO_SURF       : Net column cooling rate (LW+SW+LH+SH) between then TOA and surface, +cooling, -warming
  ;   P_minus_E                 : Precipitaion minus Evaporation rate. Note that the evaporation rate is converted from the latent heat flux
  ;   PRECT                     : Total (conv+stratiform) precipitation rate
  ;   PW			: Column precipitable water
  ;   RADCOOL_MT_TO_SURF        : Net radiative cooling rate between model top and surface, +cooling, -warming
  ;   RADCOOL_TOA_TO_SURF       : Net radiative cooling rate between the TOA and surface, +cooling, -warming
  ;   RADMT                     : Radiation balance (SW+LW) at the model top, +cooling, -warming
  ;   RADSURF                   : Surface energy imbalance (SW-LW-SH-LH), +warming, -cooling
  ;   RADTOA                    : Radiation balance (SW+LW) at the TOA, +cooling, -warming
  ;   SWCF_SURF                 : SW cloud forcing at the surface, +:warming, -: cooling
  ;   SWWARM_TOA_TO_SURF        : SW radiative warming rate between the TOA and surface, +warming, -cooling
  ;   SWWARM_MT_TO_SURF         : SW radiative warming rate between model top and surface, +warming, -cooling
  ;   zonal_mpsi                : zonal mean meridional stream function
  ;   [var]" "[number]          : a variable at a certain pressure level, e.g. "T 850", "Q 500"
  ;   [var]"__"[number]         : a variable at a certain pressure level, e.g. "T__850", "Q__500"

;--------------
; initialize
;--------------
  if ( isatt(opt, "writeout") .and. opt@writeout ) then
    file_status = "w"
  else
    file_status = "r"
  end if

  plev_zmpsi = 100.*(/10., 20., 30., 50., 70., 100., 125., 150., 175., 200., 225., \
      250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 775., 800., 825., \
      850., 875., 900., 925., 950., 975., 1000./) ; 500 Pa < p(0) < p(1) < ... < 100500. Pa

  plev_zint = (/10., 20., 30., 50., 70., 100., 125., 150., 175., 200., 225., \
      250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 775., 800., 825., \
      850., 875., 900., 925., 950., 975., 1000./) ; 500 Pa < p(0) < p(1) < ... < 100500. Pa

  ;*** open file ***
  file_temp = addfile(infilename, file_status)
  varname_input = varname

  option_diagnose = True

;---------------------------------
; variable needed interpolation
;---------------------------------

  ;*** variables for interpolation ***
  option_interp = False
  varname0 = ""

  if (.not.ismissing(str_match(varname," "))) then
    string11 = str_split(varname," ")
    option_interp = True
  else if (.not.ismissing(str_match(varname,"__"))) then
    string11 = str_split(varname,"__")
    option_interp = True
  end if
  end if
    
  if (option_interp) then
    var_interp    = string11(0)
    plev_interp   = tofloat(string11(1))
    varname = string11(0)+string11(1)

    if (var_interp .eq. "qv_zint") then
      varname0 = "qv_zint"
      varname  = "qv_zint_suf_to_"+string11(1)+"hPa"
      option_interp = False
    end if
    if (var_interp .eq. "qu_zint") then
      varname0 = "qu_zint"
      varname  = "qu_zint_suf_to_"+string11(1)+"hPa"
      option_interp = False
    end if
    if (var_interp .eq. "qwind_zint") then
      varname0 = "qwind_zint"
      varname  = "qwind_zint_suf_to_"+string11(1)+"hPa"
      option_interp = False
    end if

    if (var_interp .eq. "qudiv") then
      varname0 = "qudiv"
      option_interp = False
    end if

  end if

  ;print(var_interp)
  ;print(plev_interp)
  ;print(varname)
  ;print(option_interp)
  ;exit

  ;*** if variables need to be interpolated ***
  if (option_interp) then
    varcesm_name      = var_interp
    plev              = plev_interp
    plev@long_name    = "pressure"
    plev@units        = "hPa"

    hyai_temp = file_temp->hyai      ; hybrid A coefficient at layer interface
    hybi_temp = file_temp->hybi      ; hybrid B coefficient at layer interface
    hyam_temp = file_temp->hyam      ; hybrid A coefficient at layer midpoints
    hybm_temp = file_temp->hybm      ; hybrid B coefficient at layer midpoints
    ps_temp   = file_temp->PS        ; surface pressures  in Pa
    p0_temp   = file_temp->P0        ; reference pressure in Pa
    varcesm_temp = file_temp->$varcesm_name$
  
    ;*** interpolate ***
    varcesm_plev_temp = vinth2p(varcesm_temp, hyam_temp, hybm_temp, plev, ps_temp, 1, p0_temp/100., 1, False)
    copy_VarAtts(varcesm_temp, varcesm_plev_temp)
    varcesm_plev_temp0  = varcesm_plev_temp(:,0,:,:)  ; reduce dimesion from (time,plev,lat,lon) to (time,lat,lon)
    ;printVarSummary(varcesm_plev_temp0)

    ;*** change_var_atts ***
    res11 = True
    opt11 = True
    opt11@chg_cncolors = False
    change_var_atts(varcesm_name,varcesm_plev_temp0,res11,opt11)
    delete(res11)
    delete(opt11)

    var_temp  = varcesm_plev_temp0
    var_temp@long_name = varcesm_name+" at "+plev+" hPa"
    option_diagnose = False
  end if  ; end if of option_interp

;-------------------------------
; variables already in the file
;-------------------------------
  varname1 = varname
  if (.not.isfilevar(file_temp, varname)) then
    varname1 = ignorant_convert_varname(varname, opt)
  end if

  ;*** if varname is in read-in files ***
  if (isfilevar(file_temp, varname1)) then
    var_temp  = file_temp->$varname1$

    res11 = True
    opt11 = True
    opt11@chg_cncolors = False    
    change_var_atts(varname,var_temp,res11,opt11)
    delete(res11)
    delete(opt11)
    option_diagnose = False
  end if

;-------------------------
; diagnose each variable
;-------------------------
if (option_diagnose) then

  if ( varname .eq. "PHIS_m" ) then
    PHIS_temp = file_temp->PHIS
    PHIS_m_temp = divide(PHIS_temp, 9.8)

    var_temp  = PHIS_m_temp
    var_temp@long_name = "Surface height"
    var_temp@units     = "m"
  end if

  if ( varname .eq. "SWCF_SURF" ) then
    FSDS_temp = file_temp->FSDS
    FSDSC_temp = file_temp->FSDSC
    SWCF_SURF_temp = diff(FSDS_temp, FSDSC_temp)

    var_temp  = SWCF_SURF_temp
    var_temp@long_name = "SW cloud forcing at the surface, +:warming, -: cooling"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "FLDSC" ) then
    FLNS_temp   = file_temp->FLNS
    FLDS_temp   = file_temp->FLDS
    FLNSC_temp  = file_temp->FLNSC
    FLDSC_temp0 = plus(FLNS_temp,FLDS_temp)
    FLDSC_temp  = diff(FLDSC_temp0,FLNSC_temp)

    var_temp  = FLDSC_temp
    var_temp@long_name = "Clearsky Downwelling longwave flux at surface (assume blackbody surface)"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "P_minus_E" ) then
    if (isfilevar(file_temp,"PRECT")) then
      PRECT_temp = file_temp->PRECT
    else
      PRECC_temp = file_temp->PRECC
      PRECL_temp = file_temp->PRECL
      PRECT_temp = plus(PRECC_temp,PRECL_temp)
    end if

    PRECT_temp  = unit_convert(PRECT_temp,"m/s","mm/day")

    ; evaporation rate is converted by the latent heat flux
    LHFLX_temp = file_temp->LHFLX
    EVAP_temp  = unit_convert(LHFLX_temp , "W/m2", "mm/day")

    PmE_temp  = diff(PRECT_temp,EVAP_temp)
    var_temp  = PmE_temp
    var_temp@long_name = "Precipitaion minus Evaporation rate"
    var_temp@units     = "mm/day"
  end if

  if ( varname .eq. "EVAP" ) then
    LHFLX_temp = file_temp->LHFLX
    EVAP_temp  = unit_convert(LHFLX_temp , "W/m2", "mm/day")

    var_temp  = EVAP_temp
    var_temp@long_name = "Evaporation rate (converted from latent heat flux)"
    var_temp@units     = "mm/day"
  end if

  if ( varname .eq. "FLDSC_over_FLDS" ) then
    FLDS_temp = file_temp->FLDS
    FLDSC_temp = file_temp->FLDSC

    FLDSC_temp_over_FLDS_temp = divide(FLDSC_temp, FLDS_temp)
    var_temp  = FLDSC_temp_over_FLDS_temp
    var_temp@long_name = "frac of clrsky Dn LW flux in all-sky flux"
    var_temp@units     = "fraction"
  end if

  if ( varname .eq. "LWCF_over_FLDS" ) then
    FLDS_temp = file_temp->FLDS
    FLDSC_temp = file_temp->FLDSC
    LWCFSURF_temp = diff(FLDS_temp, FLDSC_temp)

    LWCFSURF_temp_over_FLDS_temp = divide(LWCFSURF_temp, FLDS_temp)
    var_temp  = LWCFSURF_temp_over_FLDS_temp
    var_temp@long_name = "frac of LW cloud forcing in all-sky flux"
    var_temp@units     = "fraction"
  end if
 
  if ( varname .eq. "LWCF_SURF" ) then
    FLDS_temp = file_temp->FLDS
    FLDSC_temp = file_temp->FLDSC

    LWCF_SURF_temp = diff(FLDS_temp, FLDSC_temp)
    var_temp  = LWCF_SURF_temp
    var_temp@long_name = "LW cloud forcing at the surface, +:warming, -: cooling"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "FLUS" ) then
    FLDS_temp  = file_temp->FLDS
    FLNS_temp  = file_temp->FLNS

    FLUS_temp  = plus(FLDS_temp,FLNS_temp)
    var_temp  = FLUS_temp
    var_temp@long_name = "Upward LW flux at the surface"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "FLUSC" ) then
    FLDSC_temp  = file_temp->FLDSC
    FLNSC_temp  = file_temp->FLNSC

    FLUSC_temp  = plus(FLDSC_temp,FLNSC_temp)
    var_temp  = FLUSC_temp
    var_temp@long_name = "Upward LW flux at the surface"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "FSUS" ) then
    FSDS_temp = file_temp->FSDS
    FSNS_temp = file_temp->FSNS
  
    FSUS_temp = diff(FSDS_temp,FSNS_temp)
    var_temp = FSUS_temp
    var_temp@long_name = "Upward SW flux at the surface"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "FSUSC" ) then
    FSDSC_temp = file_temp->FSDSC
    FSNSC_temp = file_temp->FSNSC
  
    FSUSC_temp = diff(FSDSC_temp,FSNSC_temp)
    var_temp = FSUSC_temp
    var_temp@long_name = "clr-sky Upward SW flux at the surface"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "FSDTOA" ) then
    FSNTOA_temp = file_temp->FSNTOA
    FSUTOA_temp = file_temp->FSUTOA

    FSDTOA_temp = plus(FSUTOA_temp,FSNTOA_temp)
    var_temp = FSDTOA_temp
    var_temp@long_name = "Downward SW flux at the TOA, i.e. SOLIN"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "RADMT" ) then
    FSNT_temp = file_temp->FSNT
    FLNT_temp = file_temp->FLNT

    RADMT_temp = diff(FSNT_temp,FLNT_temp)
    var_temp = RADMT_temp
    var_temp@long_name = "Radiation balance (SW+LW) at the model top, +:gain energy, -:lose energy"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "RADTOA" ) then
    FSNTOA_temp = file_temp->FSNTOA
    FLUT_temp = file_temp->FLUT

    RADTOA_temp = diff(FSNTOA_temp,FLUT_temp)
    var_temp = RADTOA_temp
    ;var_temp@long_name = "Radiation balance (SW+LW) at the TOA, +:gain energy, -:lose energy"
    var_temp@long_name = "TOA net downward flux"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "PRECT" ) then
    PRECC_temp  = file_temp->PRECC
    PRECL_temp  = file_temp->PRECL

    PRECT_temp  = plus(PRECC_temp,PRECL_temp)
    var_temp  = PRECT_temp
    var_temp  = unit_convert(var_temp,"m/s","mm/day")
    var_temp@long_name = "Total (conv+stratiform) precipitation rate"
  end if

  if ( varname .eq. "RADCOOL_MT_TO_SURF" ) then
    FLNT_temp  = file_temp->FLNT
    FLNS_temp  = file_temp->FLNS
    FSNT_temp  = file_temp->FSNT
    FSNS_temp  = file_temp->FSNS
  
    LWCOOL_temp  = diff(FLNT_temp,FLNS_temp)
    SWWARM_temp  = diff(FSNT_temp,FSNS_temp)
    RADCOOL_temp = diff(LWCOOL_temp,SWWARM_temp)

    var_temp = RADCOOL_temp
    var_temp@long_name = "Net radiative cooling rate between model top and surface, +:cooling, -:warming"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "LWCOOL_MT_TO_SURF" ) then
    FLNT_temp  = file_temp->FLNT
    FLNS_temp  = file_temp->FLNS
  
    LWCOOL_temp  = diff(FLNT_temp,FLNS_temp)
    var_temp = LWCOOL_temp
    var_temp@long_name = "LW radiative cooling rate between model top and surface, +:cooling, -:warming"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "SWWARM_MT_TO_SURF" ) then
    FSNT_temp  = file_temp->FSNT
    FSNS_temp  = file_temp->FSNS
  
    SWWARM_temp  = diff(FSNT_temp,FSNS_temp)
    var_temp = SWWARM_temp
    var_temp@long_name = "SW radiative warming rate between model top and surface, +:warming, -:cooling"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "RADCOOL_TOA_TO_SURF" ) then
    FLUT_temp  = file_temp->FLUT
    FLNS_temp  = file_temp->FLNS
    FSNTOA_temp  = file_temp->FSNTOA
    FSNS_temp  = file_temp->FSNS
  
    LWCOOL_temp  = diff(FLUT_temp,FLNS_temp)
    SWWARM_temp  = diff(FSNTOA_temp,FSNS_temp)
    RADCOOL_temp = diff(LWCOOL_temp,SWWARM_temp)

    var_temp = RADCOOL_temp
    var_temp@long_name = "Net radiative cooling rate between the TOA and surface, +:cooling, -:warming"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "NETCOOL_TOA_TO_SURF" ) then
    FLUT_temp  = file_temp->FLUT
    FLNS_temp  = file_temp->FLNS
    FSNTOA_temp  = file_temp->FSNTOA
    FSNS_temp  = file_temp->FSNS

    LHFLX_temp = file_temp->LHFLX  
    SHFLX_temp = file_temp->SHFLX  

    LWCOOL_temp  = diff(FLUT_temp,FLNS_temp)
    SWWARM_temp  = diff(FSNTOA_temp,FSNS_temp)
    RADCOOL_temp = diff(LWCOOL_temp,SWWARM_temp)

    SURF_temp    = plus(LHFLX_temp,SHFLX_temp)
    NEWCOOL_temp = diff(RADCOOL_temp,SURF_temp)

    var_temp = NEWCOOL_temp
    var_temp@long_name = "Column Net cooling rate between the TOA and surface, +:cooling, -:warming"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "LWCOOL_TOA_TO_SURF" ) then
    FLUT_temp  = file_temp->FLUT
    FLNS_temp  = file_temp->FLNS

    LWCOOL_temp  = diff(FLUT_temp,FLNS_temp)
    var_temp = LWCOOL_temp
    var_temp@long_name = "LW radiative cooling rate between the TOA and surface, +:cooling, -:warming"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "SWWARM_TOA_TO_SURF" ) then
    FSNTOA_temp  = file_temp->FSNTOA
    FSNS_temp  = file_temp->FSNS
  
    SWWARM_temp  = diff(FSNTOA_temp,FSNS_temp)
    var_temp = SWWARM_temp
    var_temp@long_name = "SW radiative warming rate between the TOA and surface, +:warming, -:cooling"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "RADSURF" ) then
    SHFLX_temp  = file_temp->SHFLX
    FLNS_temp   = file_temp->FLNS
  
    LHFLX_temp  = file_temp->LHFLX
    FSNS_temp   = file_temp->FSNS
  
    RADSURF_temp  = LHFLX_temp
    RADSURF_temp = FSNS_temp - FLNS_temp - SHFLX_temp - LHFLX_temp
  
    var_temp = RADSURF_temp
    ;var_temp@long_name = "Surface energy imbalance (SW-LW-SH-LH), +:gain energy, -:lose energy"
    var_temp@long_name = "Surface net downward flux (rad & turbulence)"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "ENGSURF" ) then
    SHFLX_temp  = file_temp->SHFLX
    FLNS_temp   = file_temp->FLNS
  
    LHFLX_temp  = file_temp->LHFLX
    FSNS_temp   = file_temp->FSNS
  
    ENGSURF_temp  = LHFLX_temp
    ENGSURF_temp =  FLNS_temp + SHFLX_temp + LHFLX_temp - FSNS_temp
  
    var_temp = ENGSURF_temp
    var_temp@long_name = "Surface net upward energy flux (LW+SH+LH-SW), +:upward, -:downward"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "zonal_mpsi" ) then
    varcesm_name      = "V"
    plev              = plev_zmpsi
    plev@long_name    = "pressure"
    plev@units        = "Pa"

    hyai_temp = file_temp->hyai      ; hybrid A coefficient at layer interface
    hybi_temp = file_temp->hybi      ; hybrid B coefficient at layer interface
    hyam_temp = file_temp->hyam      ; hybrid A coefficient at layer midpoints
    hybm_temp = file_temp->hybm      ; hybrid B coefficient at layer midpoints
    ps_temp   = file_temp->PS        ; surface pressures  in Pa
    p0_temp   = file_temp->P0        ; reference pressure in Pa
    varcesm_temp = file_temp->$varcesm_name$
    lat_temp     =varcesm_temp&lat

    varcesm_plev_temp = vinth2p(varcesm_temp, hyam_temp, hybm_temp, plev/100., ps_temp, 1, p0_temp/100., 1, False)
    copy_VarAtts(varcesm_temp, varcesm_plev_temp)
    ;printVarSummary(varcesm_plev_temp)

    zmpsi_temp = zonal_mpsi_Wrap(varcesm_plev_temp, lat_temp, plev, ps_temp)   
    var_temp  = zmpsi_temp
    var_temp@long_name = "zonal mean meridional stream function"
    var_temp@units     = "kg/s"
  end if

  ;*** newcase ***
  if ( varname .eq. "PW" ) then 
    Q_temp = file_temp->Q

    hyai_temp = file_temp->hyai      ; hybrid A coefficient at layer interface
    hybi_temp = file_temp->hybi      ; hybrid B coefficient at layer interface
    ps_temp   = file_temp->PS        ; surface pressures  in Pa
    p0_temp   = file_temp->P0        ; reference pressure in Pa

    dp_temp = dpres_hybrid_ccm(ps_temp, p0_temp, hyai_temp, hybi_temp)
    dp_temp!0 = "time"
    dp_temp!1 = "lev"
    dp_temp!2 = "lat"
    dp_temp!3 = "lon"

    pw_temp = prcwater_dp(Q_temp(time|:,lat|:,lon|:,lev|:),dp_temp(time|:,lat|:,lon|:,lev|:))    
    copy_VarCoords(Q_temp(:,0,:,:),pw_temp)
   
    pw_temp@units = "kg/m2"
    pw_temp = unit_convert(pw_temp, "kg/m2", "mm")
    var_temp = pw_temp
    var_temp@long_name = "total column precipitable water"
    ;printVarSummary(var_temp)
  end if

  ; RRTMG_LW band widths
  ;   band 01: 10-350 cm~S~-1
  ;   band 02: 350-500 cm~S~-1
  ;   band 03: 500-630 cm~S~-1
  ;   band 04: 630-700 cm~S~-1
  ;   band 05: 700-820 cm~S~-1
  ;   band 06: 820-980 cm~S~-1
  ;   band 07: 980-1080 cm~S~-1
  ;   band 08: 1080-1180 cm~S~-1
  ;   band 09: 1180-1390 cm~S~-1
  ;   band 10: 1390-1480 cm~S~-1
  ;   band 11: 1480-1800 cm~S~-1
  ;   band 12: 1800-2080 cm~S~-1
  ;   band 13: 2080-2250 cm~S~-1
  ;   band 14: 2250-2380 cm~S~-1
  ;   band 15: 2380-2600 cm~S~-1
  ;   band 16: 2600-3250 cm~S~-1

  if ( varname .eq. "FLUS_RRTMG_10-630" ) then
    wn_start = 10
    wn_end   = 630
    vars1  = (/"FLSU01","FLSU02","FLSU03"/)
    nvar1  = dimsizes(vars1)
    FLSUXX = file_temp->$vars1(0)$
    FLSUXX = 0.

    do vv=0,nvar1-1
      var1 = vars1(vv)
      FLSUX1  = file_temp->$var1$
      FLSUXX  = plus(FLSUXX,FLSUX1)
      delete(FLSUX1)
    end do

    kk_surf = dimsizes(FLSUXX&ilev)
    var_temp = FLSUXX(:,kk_surf-1,:,:)
    var_temp@long_name = "Upward LW flux at the surface from "+wn_start+" to "+wn_end+" cm-1"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "FLUS_RRTMG_820-1180" ) then
    wn_start = 820
    wn_end   = 1180
    vars1  = (/"FLSU06","FLSU07","FLSU08"/)
    nvar1  = dimsizes(vars1)
    FLSUXX = file_temp->$vars1(0)$
    FLSUXX = 0.

    do vv=0,nvar1-1
      var1 = vars1(vv)
      FLSUX1  = file_temp->$var1$
      FLSUXX  = plus(FLSUXX,FLSUX1)
      delete(FLSUX1)
    end do

    kk_surf = dimsizes(FLSUXX&ilev)
    var_temp = FLSUXX(:,kk_surf-1,:,:)
    var_temp@long_name = "Upward LW flux at the surface from "+wn_start+" to "+wn_end+" cm-1"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "FLUT_RRTMG_10-630" ) then
    wn_start = 10
    wn_end   = 630
    vars1  = (/"FLSU01","FLSU02","FLSU03"/)
    nvar1  = dimsizes(vars1)
    FLSUXX = file_temp->$vars1(0)$
    FLSUXX = 0.

    do vv=0,nvar1-1
      var1 = vars1(vv)
      FLSUX1  = file_temp->$var1$
      FLSUXX  = plus(FLSUXX,FLSUX1)
      delete(FLSUX1)
    end do

    kk_toa = 0
    var_temp = FLSUXX(:,kk_toa,:,:)
    var_temp@long_name = "Upward LW flux at the TOA from "+wn_start+" to "+wn_end+" cm-1"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "FLUT_RRTMG_820-1180" ) then
    wn_start = 820
    wn_end   = 1180
    vars1  = (/"FLSU06","FLSU07","FLSU08"/)
    nvar1  = dimsizes(vars1)
    FLSUXX = file_temp->$vars1(0)$
    FLSUXX = 0.

    do vv=0,nvar1-1
      var1 = vars1(vv)
      FLSUX1  = file_temp->$var1$
      FLSUXX  = plus(FLSUXX,FLSUX1)
      delete(FLSUX1)
    end do

    kk_toa = 0
    var_temp = FLSUXX(:,kk_toa,:,:)
    var_temp@long_name = "Upward LW flux at the TOA from "+wn_start+" to "+wn_end+" cm-1"
    var_temp@units     = "W/m2"
  end if

  if ( varname .eq. "MET" ) then

    ; Meridional energy transport is based on Eq (2), Frierson and Hwang (2011, J. Climate)
    pi = acos(-1.)
    a  =  6371. * 1000. ; earth radius (m)

    FSNTOA_temp = file_temp->FSNTOA
    FLUT_temp   = file_temp->FLUT

    FLNS_temp   = file_temp->FLNS
    FSNS_temp   = file_temp->FSNS
    LHFLX_temp   = file_temp->LHFLX
    SHFLX_temp   = file_temp->SHFLX

    ;*** energy dicifit in a column ***
    Fa_temp = FSNTOA_temp
    Fa_temp = (FSNTOA_temp-FLUT_temp) - (FSNS_temp-FLNS_temp-LHFLX_temp-SHFLX_temp)
    ;printVarSummary(Fa_temp)

    ;*** get lat/lon grid size, dlat & dlon ***
    lat = tofloat(Fa_temp&lat)  ; latitude must in -90 to 90
    lon = tofloat(Fa_temp&lon)  ; latitude must in -90 to 90
    nlat = dimsizes(lat)
    coslat = abs(cos(lat/180.*pi))

    dlon_rad = (lon(1)-lon(0))/180.*pi
    dlat_rad = (lat(1)-lat(0))/180.*pi

    ;*** change unit from W/m2 to W ***
    Fa_temp_W = Fa_temp
    do j=0,nlat-1
      Fa_temp_W(:,j,:) = Fa_temp(:,j,:) * coslat(j) * dlon_rad * dlat_rad * a^2
    end do

    ;*** summation over longitude ***
    Fa_temp_W_isum = dim_sum_n_Wrap(Fa_temp_W,2)

    ;*** calculate meridional energy transport ***
    Fa_temp_W_isum_tran = Fa_temp_W_isum
    Fa_temp_W_isum_tran(:,0) = 0. ;Fa_temp_W_isum

    do j=1,nlat-1
      Fa_temp_W_isum_tran(:,j) = Fa_temp_W_isum_tran(:,j-1) + Fa_temp_W_isum(:,j)
    end do

    ;*** change unit from W to PW, i.e. 10^15 W ***
    Fa_temp_W_isum_tran = Fa_temp_W_isum_tran / 10.^15

    ;*** return variable ***
    var_temp = Fa_temp_W_isum_tran
    var_temp@long_name = "Meridional energy transport"
    var_temp@units     = "PW"
    ;printVarSummary(var_temp)
  end if

  if ( varname0 .eq. "qudiv") then
    varname_q         = "Q"
    varname_u         = "U"
    varname_v         = "V"

    hyai_temp = file_temp->hyai      ; hybrid A coefficient at layer interface
    hybi_temp = file_temp->hybi      ; hybrid B coefficient at layer interface
    hyam_temp = file_temp->hyam      ; hybrid A coefficient at layer midpoints
    hybm_temp = file_temp->hybm      ; hybrid B coefficient at layer midpoints
    ps_temp   = file_temp->PS        ; surface pressures  in Pa
    p0_temp   = file_temp->P0        ; reference pressure in Pa
    varcesm_q = file_temp->$varname_q$
    varcesm_v = file_temp->$varname_v$
    varcesm_u = file_temp->$varname_u$

    lat = varcesm_q&lat
    lon = varcesm_q&lon

;printVarSummary(varcesm_q)
;printVarSummary(varcesm_u)
;printVarSummary(varcesm_v)
;print(plev_interp)

    varcesm_q_plev = vinth2p(varcesm_q, hyam_temp, hybm_temp, plev_interp, ps_temp, 1, p0_temp/100., 1, False)
    copy_VarAtts(varcesm_q, varcesm_q_plev)

    varcesm_v_plev = vinth2p(varcesm_v, hyam_temp, hybm_temp, plev_interp, ps_temp, 1, p0_temp/100., 1, False)
    copy_VarAtts(varcesm_v, varcesm_v_plev)

    varcesm_u_plev = vinth2p(varcesm_u, hyam_temp, hybm_temp, plev_interp, ps_temp, 1, p0_temp/100., 1, False)
    copy_VarAtts(varcesm_u, varcesm_u_plev)

;printVarSummary(varcesm_q_plev)
;printVarSummary(varcesm_u_plev)
;printVarSummary(varcesm_v_plev)

    ntime = dimsizes(varcesm_q_plev&time)
    nlev  = dimsizes(varcesm_q_plev&lev_p)
    nlat  = dimsizes(varcesm_q_plev&lat)
    nlon  = dimsizes(varcesm_q_plev&lon)

    varcesm_qv_plev = varcesm_q_plev
    varcesm_qv_plev = 0.
    varcesm_qu_plev = varcesm_qv_plev

    do tt=0,ntime-1
    do pp=0,nlev-1
    do jj=0,nlat-1
    do ii=0,nlon-1
      varcesm_qv_plev(tt,pp,jj,ii) = varcesm_q_plev(tt,pp,jj,ii) * varcesm_v_plev(tt,pp,jj,ii)
      varcesm_qu_plev(tt,pp,jj,ii) = varcesm_q_plev(tt,pp,jj,ii) * varcesm_u_plev(tt,pp,jj,ii)
      ;print("q="+varcesm_q_plev(tt,pp,jj,ii)+", v="+varcesm_wind_plev(tt,pp,jj,ii)+", q*v="+varcesm_qv_plev(tt,pp,jj,ii))
    end do
    end do
    end do
    end do

    boundOpt = 0 ; Boundary points are set to the missing value.
    qudiv = uv2dv_cfd(varcesm_qu_plev, varcesm_qv_plev, lat, lon, boundOpt)
    ;qudiv = uv2dvF_Wrap(varcesm_qu_plev, varcesm_qv_plev) ; for global data and non-missing values only

    copy_VarCoords(varcesm_q_plev, qudiv)
    ;printVarSummary(qudiv)

;print(qudiv(10,0,5,5))
    qudiv = -multiply(qudiv, 10.^8)
;print(qudiv(10,0,5,5))
    var_temp = qudiv(:,0,:,:)
    ;var_temp@long_name = "-~F34~QW~F~(qv) at "+plev_interp+" hPa"
    var_temp@long_name = "-~F34~QW~F~(qU) at "+plev_interp+"hPa"
    ;var_temp@long_name = "-~F34~QW~F~(qv) & wind at "+plev_interp+" hPa"
    var_temp@units     = "10~S~-5~N~g kg~S~-1~N~s~S~-1"
    ;printVarSummary(var_temp)

  end if


  if ( varname0 .eq. "qv_zint" .or. varname0 .eq. "qu_zint" .or. varname0 .eq. "qwind_zint") then
    if (varname0 .eq. "qv_zint") then
      ;varname_v         = "V"
      long_name = "q*v, int from surf to "+plev_interp+" hPa"
    end if
    if (varname0 .eq. "qu_zint") then
      ;varname_v         = "U"
      long_name = "q*u, int from surf to "+plev_interp+" hPa"
    end if
    if (varname0 .eq. "qwind_zint") then
      long_name = "q*|u|, int from surf to "+plev_interp+" hPa"
    end if

    varname_q         = "Q"
    plev              = plev_zint
    plev@long_name    = "pressure"
    plev@units        = "hPa"

    varname_u         = "U"
    varname_v         = "V"

    hyai_temp = file_temp->hyai      ; hybrid A coefficient at layer interface
    hybi_temp = file_temp->hybi      ; hybrid B coefficient at layer interface
    hyam_temp = file_temp->hyam      ; hybrid A coefficient at layer midpoints
    hybm_temp = file_temp->hybm      ; hybrid B coefficient at layer midpoints
    ps_temp   = file_temp->PS        ; surface pressures  in Pa
    p0_temp   = file_temp->P0        ; reference pressure in Pa
    varcesm_q = file_temp->$varname_q$
    varcesm_v = file_temp->$varname_v$
    varcesm_u = file_temp->$varname_u$

    varcesm_q_plev = vinth2p(varcesm_q, hyam_temp, hybm_temp, plev, ps_temp, 1, p0_temp/100., 1, False)
    copy_VarAtts(varcesm_q, varcesm_q_plev)

    varcesm_v_plev = vinth2p(varcesm_v, hyam_temp, hybm_temp, plev, ps_temp, 1, p0_temp/100., 1, False)
    copy_VarAtts(varcesm_v, varcesm_v_plev)

    varcesm_u_plev = vinth2p(varcesm_u, hyam_temp, hybm_temp, plev, ps_temp, 1, p0_temp/100., 1, False)
    copy_VarAtts(varcesm_u, varcesm_u_plev)

    if (varname0 .eq. "qv_zint") then
      varcesm_wind_plev = varcesm_v_plev
    end if
    if (varname0 .eq. "qu_zint") then
      varcesm_wind_plev = varcesm_u_plev
    end if
    if (varname0 .eq. "qwind_zint") then
      varcesm_wind_plev = wind_speed(varcesm_u_plev,varcesm_u_plev)
    end if

    ntime = dimsizes(varcesm_q_plev&time)
    nlev  = dimsizes(varcesm_q_plev&lev_p)
    nlat  = dimsizes(varcesm_q_plev&lat)
    nlon  = dimsizes(varcesm_q_plev&lon)

    varcesm_qv_plev = varcesm_q_plev
    varcesm_qv_plev = 0.

    do tt=0,ntime-1
    do pp=0,nlev-1
    do jj=0,nlat-1
    do ii=0,nlon-1
      varcesm_qv_plev(tt,pp,jj,ii) = varcesm_q_plev(tt,pp,jj,ii) * varcesm_wind_plev(tt,pp,jj,ii)
      ;print("q="+varcesm_q_plev(tt,pp,jj,ii)+", v="+varcesm_wind_plev(tt,pp,jj,ii)+", q*v="+varcesm_qv_plev(tt,pp,jj,ii))
    end do
    end do
    end do
    end do

    ; vertical integration of moisture divergence via NCL function vibeta (http://www.ncl.ucar.edu/Document/Functions/Built-in/vibeta.shtml)
    ; note: 1. ptop MUST be the minimum of p, 2. psfc uses surface pressure, 3. vertical integration from psfc to ptop
    linlog = 1         ; 1 = linear interpolation, 2 = log interpolation
    pbot = 1100.       ; unit: hPa
    ptop = plev_interp ; unit: hPa    
    varcesm_qv_plev_reverse = varcesm_qv_plev(:,::-1,:,:)
    varcesm_qv_plev_zint= vibeta (plev(::-1),varcesm_qv_plev_reverse(time|:,lat|:,lon|:,lev_p|:),linlog,ps_temp/100.,pbot,ptop) ; in unit: kg/kg * m/s * hPa
    varcesm_qv_plev_zint = divide(varcesm_qv_plev_zint, 0.01*9.8) ; change unit to kg/m/s

    copy_VarCoords(ps_temp, varcesm_qv_plev_zint)
    ;printVarSummary(varcesm_qv_plev_zint)

    ;printVarSummary(varcesm_qv_plev)
    ;printVarSummary(varcesm_q_plev)
    ;printVarSummary(varcesm_wind_plev)

    ;tt = 0
    ;pp = 5
    ;ii = 10
    ;jj = 7
    ;print(varcesm_qv_plev_zint(tt,jj,ii))
    ;print("q="+varcesm_q_plev(tt,:,jj,ii)+", v="+varcesm_wind_plev(tt,:,jj,ii)+", qv="+varcesm_qv_plev(tt,:,jj,ii)+", qv_zint="+qv_zint)

    var_temp = varcesm_qv_plev_zint
    ;var_temp@long_name = "vertically-integrated q*v, from surf to "+plev_interp+" hPa"
    var_temp@long_name = long_name
    var_temp@units     = "kg/m/s"
    var_temp@z_integral_range  = "surface to "+plev_interp+" hPa"
    ;printVarSummary(var_temp)

  end if

;  ;*** newcase ***
;  if ( varname .eq. "" ) then 
;    = file_temp->
;
;    var_temp =
;    var_temp@long_name = ""
;    var_temp@units     = ""
;    printVarSummary(var_temp)
;  end if
end if ; end if of option_diagnose

;------------------------
; write out to the file
;------------------------
  if ( isatt(opt, "writeout") .and. opt@writeout ) then
    print("")
    print("*****************************")
    print("write ["+varname+"] into the file ["+infilename+"] ...")
    file_temp->$varname$ = var_temp
    print("Done. write ["+varname+"] into the file ["+infilename+"]")
    print("*****************************")
  end if

;------------------------------------------------
; modify return variables according to land mask
;------------------------------------------------
  ;if ( isatt(opt, "select_landfrac") .and. opt@select_landfrac ) then
  if ( isatt(opt, "landfrac")) then
    select_landfrac(var_temp, opt)
  end if

;-----------------
; return variable
;-----------------

  if (isdefined("var_temp")) then
    varname = varname_input
    return(var_temp)
  else
    varname = varname_input
    print("")
    print("ERROR: (ignorant_diagnose_var) input varialbe ["+varname+"] is not defined in this function")
    print("")
    return(0)
  end if 
end

;### end of ignorant_diagnose_var ###

;#############################################
undef("ignorant_get_var_files")
function ignorant_get_var_files(infile_names:string, read_var_name[1]:string, opt)

begin

;----------------------------------------------
; determinie input files $ varible properties
;----------------------------------------------

; number of read-in files/varible
  file_num = dimsizes(infile_names)
  var_num  = dimsizes(read_var_name)

; properties of read-in varibles
  opt_diag = opt
  var = ignorant_diagnose_var(infile_names(0),read_var_name,opt_diag)

  ;tmp_infile   = addfile( infile_names(0),"r" )
  ;if ( isfilevar(tmp_infile,read_var_name) ) then
  ;  var          = tmp_infile->$read_var_name$
  ;else
  ;  opt_diag = True
  ;  var = ignorant_diagnose_var(infile_names(0),read_var_name,opt_diag)
  ;end if

  ;var_type     = getfilevartypes(tmp_infile, read_var_name)
  ;var_dimnames = getfilevardims(tmp_infile, read_var_name)
  ;var_dimsizes = getfilevardimsizes( tmp_infile, read_var_name )
  ;var_atts     = getfilevaratts(tmp_infile, read_var_name)

  var_type     = typeof(var) ;getfilevartypes(tmp_infile, read_var_name)
  var_dimnames = getvardims(var) ;(tmp_infile, read_var_name)
  var_dimsizes = dimsizes(var) ;( tmp_infile, read_var_name )
  var_atts     = getvaratts(var) ;(tmp_infile, read_var_name)

  var_dimnames_sizes = dimsizes(var_dimnames)
  var_atts_sizes     = dimsizes(var_atts) 

;------------------------------------
; determine output array properties
;------------------------------------

; define outpuy array size
  aa = 1 + var_dimnames_sizes
  work_dimsizes = new( aa, "integer" )

  work_dimsizes(0) = file_num
  do i = 1, aa-1
    work_dimsizes(i) = var_dimsizes(i-1)
  end do

  work_var = new( work_dimsizes, var_type )

; define output array named dimensions
  work_var!0 = "file"
  do i = 1, var_dimnames_sizes
    work_var!i = var_dimnames(i-1)

    ; copy coordinate
    if (iscoord(var,var_dimnames(i-1))) then
      work_var&$var_dimnames(i-1)$ = var&$var_dimnames(i-1)$
    end if
  end do

; define output array attributes
  file_att = new(1, "string")
  file_att = " "
  do i = 0, file_num-1
    aptr = tostring(i)
    file_att = file_att + aptr + ":[" + infile_names(i) + "], "
  end do

  work_var@file = file_att

  do i = 0, var_atts_sizes-1
    work_var@$var_atts(i)$ = var@$var_atts(i)$
  end do

;------------------------------------------------------
; read-in varibles from files & store into output array
;------------------------------------------------------ 

  do ff = 0, file_num-1

    opt_diag = opt 
    tmp_infile_name = infile_names(ff)
    ;print(opt_diag)
    tmp_var = ignorant_diagnose_var(tmp_infile_name,read_var_name,opt_diag)

    ;tmp_infile_name = infile_names(ff)
    ;tmp_infile  = addfile( tmp_infile_name, "r" )   
    ;if ( isfilevar(tmp_infile,read_var_name) ) then
    ;  tmp_var     = tmp_infile->$read_var_name$
    ;else
    ;  opt_diag = True
    ;  tmp_var = ignorant_diagnose_var(tmp_infile_name,read_var_name,opt_diag)
    ;end if

    if ( var_dimnames_sizes .eq. 1 ) then
      work_var(ff,:) = (/tmp_var/)
    end if
    if ( var_dimnames_sizes .eq. 2 ) then
      work_var(ff,:,:) = (/tmp_var/)
    end if

    if ( var_dimnames_sizes .eq. 3 ) then
      work_var(ff,:,:,:) = (/tmp_var/)
    end if

    if ( var_dimnames_sizes .eq. 4 ) then
      work_var(ff,:,:,:,:) = (/tmp_var/)
    end if 

  end do

;------------------------
; return to main program
;------------------------

  return(work_var)

end

;### end of ignorant_get_var_files ###

;#############################################
undef("ignorant_get_vars_file")
function ignorant_get_vars_file (infilenames[1]:string, varvars[*]:string, opt:logical)
begin

;----------------------------------------------
; determinie input files $ varible properties
;----------------------------------------------

; number of read-in files/varible
  nfile = dimsizes(infilenames)
  nvar  = dimsizes(varvars)

; properties of read-in varibles
  file1 = infilenames(0)
  var1  = varvars(0)

  opt_diag = opt
  var = ignorant_diagnose_var(file1,var1,opt_diag)

  ;tmp_infile   = addfile( file1,"r" )
  ;if ( isfilevar(tmp_infile,var1) ) then
  ;  var          = tmp_infile->$var1$
  ;else
  ;  opt_diag = True
  ;  var = ignorant_diagnose_var(file1,var1,opt_diag)
  ;end if

  ;var_type     = getfilevartypes(tmp_infile, var1)
  ;var_dimnames = getfilevardims(tmp_infile, var1)
  ;var_dimsizes = getfilevardimsizes( tmp_infile, var1 )
  ;var_atts     = getfilevaratts(tmp_infile, var1)

  var_type     = typeof(var) ;getfilevartypes(tmp_infile, read_var_name)
  var_dimnames = getvardims(var) ;(tmp_infile, read_var_name)
  var_dimsizes = dimsizes(var) ;( tmp_infile, read_var_name )
  var_atts     = getvaratts(var) ;(tmp_infile, read_var_name)

  var_dimnames_sizes = dimsizes(var_dimnames)
  var_atts_sizes     = dimsizes(var_atts)

;------------------------------------
; determine output array properties
;------------------------------------
;
; define outpuy array size
  aa = 1 + var_dimnames_sizes
  work_dimsizes = new( aa, "integer" )
;
  work_dimsizes(0) = nvar
  do i = 1, aa-1
    work_dimsizes(i) = var_dimsizes(i-1)
  end do
;
  work_var = new( work_dimsizes, var_type )
;
; define output array named dimensions
;
  work_var!0 = "variable"
  do i = 1, var_dimnames_sizes
    work_var!i = var_dimnames(i-1)
;
    ; copy coordinate
    if (iscoord(var,var_dimnames(i-1))) then
      work_var&$var_dimnames(i-1)$ = var&$var_dimnames(i-1)$
    end if
  end do
;
; define output array attributes
;
  file_att = new(nvar, "string")
  name_att = (/"var_name","long_name","units","Sampling_Sequence","cell_methods"/)
  natt     = dimsizes(name_att)

  do i=0,natt-1
    att1 = name_att(i)
    work_var@$att1$ = file_att
  end do

  do vv = 0,nvar-1
    work_var@var_name(vv) = varvars(vv)
  end do

  ;printVarSummary(work_var)

;---------------------------------------------------------------------
; read-in all varibles from the input file & store into output array
;---------------------------------------------------------------------
  ff = 0

  do vv = 0, nvar-1
    tmp_infile_name = infilenames(ff)
    varname     = varvars(vv)
    opt_diag = opt
    tmp_var = ignorant_diagnose_var(tmp_infile_name,varname,opt_diag)

    ;tmp_infile  = addfile( tmp_infile_name, "r" )
    ;if ( isfilevar(tmp_infile,varname) ) then
    ;  tmp_var          = tmp_infile->$varname$
   ;   res1 = True
   ;   opt_chgatt1 = True
   ;   opt_chgatt1@chg_cncolors = False
   ;   change_var_atts(varname, tmp_var, res1, opt_chgatt1)
   ; else
   ;   opt_diag = True
   ;   tmp_var = ignorant_diagnose_var(tmp_infile_name,varname,opt_diag)
   ; end if
    ;
    ; read and save variable values
    ;
    if ( var_dimnames_sizes .eq. 1 ) then
      work_var(vv,:) = (/tmp_var/)
    end if
    if ( var_dimnames_sizes .eq. 2 ) then
      work_var(vv,:,:) = (/tmp_var/)
    end if

    if ( var_dimnames_sizes .eq. 3 ) then
      work_var(vv,:,:,:) = (/tmp_var/)
    end if

    if ( var_dimnames_sizes .eq. 4 ) then
      work_var(vv,:,:,:,:) = (/tmp_var/)
    end if

    ;
    ; read and save variable attributes
    ;
    do aa=0,natt-1
      att1 = name_att(aa)
      if ( isatt(tmp_var,att1) ) then
        work_var@$att1$(vv) = tmp_var@$att1$
      end if
    end do

    delete(tmp_var)
  end do

;----------
; return
;----------
  return(work_var)
  ;printVarSummary(work_var)

end
;### end of ignorant_get_vars_file ###

;#############################################
undef("ignorant_convert_2Ddata")
function ignorant_convert_2Ddata (infilename[1]:string, varvar[1]:string, opt:logical)
begin

;---------------
; default setup
;---------------
  option_execute = False

  case_names = (/ \
                "tlji2pt_avgji", \
                "tlji2tji_intp", \
                "tlji2ji_intp_avgt", \
                "tlji2pn_crossji", \
                "tlji2pj_avgti", \
                "tlji2tpj_avgi", \
                "tji2jt_avgi", \
                "tji2ji_avgt", \
                "tpj2pj_avgt", \
                "monthly_series_tlji2ypm_avgji", \
                "monthly_series_tji2yjm_avgi", \
                "monthly_series_tji2ymji", \
                "sdlkfjsldjflsjdflsdjfl" \  
              /)

  plev_default = (/1., 2., 3., 5., 7., 10., 20., 30., 50., 70., 100., 125., 150., 175., 200., 225., \
      250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 775., 800., 825., \
      850., 875., 900., 925., 950., 975., 1000./)

  if (isatt(opt, "return_type")) then
    if (opt@return_type .eq. "pressure_time-latlon_avg") then
      opt@return_type = "tlji2pt_avgji"
    end if

    if (opt@return_type .eq. "pressure_lat-timelon_avg") then
      opt@return_type = "tji2jt_avgi"
    end if

    if (opt@return_type .eq. "lat_time-lon_avg") then
      opt@return_type = "tji2jt_avgi"
    end if

    if (opt@return_type .eq. "time_lat_lon-plev") then
      opt@return_type = "tlji2tji_intp"
    end if
  end if

;----------------
; read variable
;----------------
  filetemp = addfile(infilename, "r")

  opt_diag = opt
  varvar_tmp = ignorant_diagnose_var(infilename,varvar,opt_diag)

  ;if ( isfilevar(filetemp,varvar) ) then
  ;  varvar_tmp = filetemp->$varvar$
  ;else
  ;  opt_diag = opt
  ;  varvar_tmp = ignorant_diagnose_var(infilename,varvar,opt_diag)
  ;end if

;--------------------------------------------
; convert (time,lev/ilev,lat,lon) to (plev,time)
;   1. interpolate CESM hybrid level to the given pressure levels
;   2. do lat & lon average
;--------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tlji2pt_avgji") then

    ;*** read the variable ***
    hyai_tmp   = filetemp->hyai      ; hybrid A coefficient at layer interface
    hybi_tmp   = filetemp->hybi      ; hybrid B coefficient at layer interface
    hyam_tmp   = filetemp->hyam      ; hybrid A coefficient at layer midpoints
    hybm_tmp   = filetemp->hybm      ; hybrid B coefficient at layer midpoints
    PS_tmp     = filetemp->PS        ; surface pressures  in Pa
    P0_tmp     = filetemp->P0        ; reference pressure in Pa

    if (iscoord(varvar_tmp,"lev")) then
      hyaw_tmp  = hyam_tmp
      hybw_tmp  = hybm_tmp
    else if (iscoord(varvar_tmp,"ilev")) then
      hyaw_tmp  = hyai_tmp
      hybw_tmp  = hybi_tmp
    end if
    end if

    ;*** set pressure levels for interpolation ***
    if ( isatt(opt, "plev") ) then
      plev_work = opt@plev
    else
      plev_work = plev_default
    end if

    ; do interpolation
    varvar_tmp_plev = vinth2p(varvar_tmp, hyaw_tmp, hybw_tmp, plev_work, PS_tmp, 1, P0_tmp/100., 1, False)
    copy_VarAtts(varvar_tmp, varvar_tmp_plev)

    varvar_tmp_plev_ijavg = lat_wgt_avg(varvar_tmp_plev)

    ;*** return variable ***
    varvar_return = varvar_tmp_plev_ijavg(lev_p|:,time|:)

    varvar_return&lev_p@long_name = "Pressure"
    varvar_return&lev_p@units     = "hPa"

    ntime = dimsizes(varvar_return&time)
    time_coord = fspan(0.,tofloat(ntime-1),ntime)
    ;time_coord = ispan(0,ntime-1,1)
    delete(varvar_return&time)
    varvar_return&time = time_coord

    ;printVarSummary(varvar_tmp_plev)
    ;printVarSummary(varvar_tmp_plev_ijavg)
    ;printVarSummary(varvar_return)

    option_execute = True
  end if  ; end if of "tlji2pt_avgji"

;--------------------------------------------
; convert (time,lev/ilev,lat,lon) to (lat,lon) at a pressure level & do time avg
;   1. interpolate CESM hybrid level to the given pressure level
;   2. do time average
;--------------------------------------------

  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tlji2ji_intp_avgt") then

    hyai_tmp   = filetemp->hyai      ; hybrid A coefficient at layer interface
    hybi_tmp   = filetemp->hybi      ; hybrid B coefficient at layer interface
    hyam_tmp   = filetemp->hyam      ; hybrid A coefficient at layer midpoints
    hybm_tmp   = filetemp->hybm      ; hybrid B coefficient at layer midpoints
    PS_tmp     = filetemp->PS        ; surface pressures  in Pa
    P0_tmp     = filetemp->P0        ; reference pressure in Pa

    if (iscoord(varvar_tmp,"lev")) then
      hyaw_tmp  = hyam_tmp
      hybw_tmp  = hybm_tmp
    else if (iscoord(varvar_tmp,"ilev")) then
      hyaw_tmp  = hyai_tmp
      hybw_tmp  = hybi_tmp
    end if
    end if

    if (isatt(opt,"time_step") .and. opt@time_step .eq. "avg") then
      tt = -999
    end if

    ;*** set pressure levels for interpolation ***
    if ( isatt(opt, "plev") .and. dimsizes(opt@plev) .eq. 1) then
      plev_work = opt@plev
    else
      print("WARNING: (ignorant_convert_2Ddata) case [tlji2ji_intp_avgt] only allow a pressure level or opt@plev is not given. Set to 500 hPa")
      plev_work = 500.
    end if

    ; do interpolation
    varvar_tmp_plev = vinth2p(varvar_tmp, hyaw_tmp, hybw_tmp, plev_work, PS_tmp, 1, P0_tmp/100., 1, False)
    copy_VarAtts(varvar_tmp, varvar_tmp_plev)

    if (isatt(opt,"time_step") .and. opt@time_step .eq. "avg") then
      varvar_return = dim_avg_n_Wrap(varvar_tmp_plev(:,0,:,:),0) ; return (lat,lon) averaged over time
    else
      varvar_return0 = varvar_tmp_plev(:,0,:,:)  ; return (time,lat,lon)
    end if

    varvar_return = dim_avg_n_Wrap(varvar_return0,0)

    delete(varvar_return@lev_p)
    varvar_return@pressure_level = plev_work+" hPa"
    varvar_return@long_name = varvar_return@long_name+" at "+varvar_return@pressure_level

    option_execute = True
  end if  ; end if of "tlji2ji_intp_avgt"

;--------------------------------------------
; convert (time,lat,lon) to (lat,time)
;   1. do lon average
;--------------------------------------------

  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tji2jt_avgi") then

    varvar_tmp_lonavg = dim_avg_n_Wrap(varvar_tmp,2)
    ;printVarSummary(varvar_tmp_lonavg)
  
    latname = varvar_tmp!1    ; get the name of latitude dimension
    timename = varvar_tmp!0   ; get the name of time     dimension

    varvar_return = varvar_tmp_lonavg($latname$|:,$timename$|:) 

    ntime = dimsizes(varvar_return&$timename$)
    time_coord = fspan(0.,tofloat(ntime-1),ntime)
    delete(varvar_return&time)
    varvar_return&time = time_coord    

    option_execute = True
    ;printVarSummary(varvar_return)
  end if  ; end if of "tji2jt_avgi"

;--------------------------------------------
; convert (time,lat,lon) to (lat,lon)
;   1. do time average
;--------------------------------------------

  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tji2ji_avgt") then

    varvar_tmp_tavg = dim_avg_n_Wrap(varvar_tmp,0)
    ;printVarSummary(varvar_tmp_lonavg)
  
    varvar_return = varvar_tmp_tavg

    option_execute = True
    ;printVarSummary(varvar_return)
  end if  ; end if of "tji2jt_avgi"

;--------------------------------------------
; OUTDATED: convert (time,lat,lon1) to (time,lat,lon2), lon1 is not monotonic while lon2 is. e.g. lon1 = (/350,0,10/) and lon2 is (/-10,0,10/).
;--------------------------------------------
;  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "time_lat_lon-lon_mono") then
;
;    ;*** read the variable ***
;    filetemp = addfile(infilename, "r")
;      varvar_tmp = filetemp->$varvar$
;
;    ;if ( isatt(opt, "time_step") .and. isint(opt@time_step)) then
;    ;  tt = opt@time_step
;    ;  varvar_tmp = filetemp->$varvar$(tt,:,:)
;    ;else
;    ;  varvar_tmp = filetemp->$varvar$
;    ;end if
;
;    ;*** get lon ***
;    if ( iscoord(varvar_tmp,"lon")) then
;      lonname = "lon"
;      lon = filetemp&$lonname$
;      option_execute = True
;
;    else if ( iscoord(varvar_tmp,"longitude")) then
;      lonname = "longitude"
;      lon = filetemp&$lonname$
;      option_execute = True
;    else
;      print("ERROR: (ignorant_convert_2Ddata) coordinate lon/longitude is not found in return type: [time_lat_lon-lon_mono]")
;      print("")
;      option_execute = False
;    end if
;    end if
;
;    ;*** make lon monotonic ***
;    if (option_execute) then
;
;      nlon = dimsizes(lon)
;      lon1 = lon
;      do i=0,nlon-1
;        if (lon(i).gt.180.) then
;          lon1(i) = lon(i) - 360.
;        end if
;      end do
;
;    ;*** return variable ***
;    varvar_return = varvar_tmp
;    varvar_return&$lonname$ = (/lon1/)
;
;    end if
;  end if ; end if of "time_lat_lon-lon_mono"

;--------------------------------------------
; convert (time,lev/ilev,lat,lon) to (time,lat,lon) at a pressure level
;   1. interpolate CESM hybrid level to the given pressure level
;--------------------------------------------

  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tlji2tji_intp") then

    hyai_tmp   = filetemp->hyai      ; hybrid A coefficient at layer interface
    hybi_tmp   = filetemp->hybi      ; hybrid B coefficient at layer interface
    hyam_tmp   = filetemp->hyam      ; hybrid A coefficient at layer midpoints
    hybm_tmp   = filetemp->hybm      ; hybrid B coefficient at layer midpoints
    PS_tmp     = filetemp->PS        ; surface pressures  in Pa
    P0_tmp     = filetemp->P0        ; reference pressure in Pa

    if (iscoord(varvar_tmp,"lev")) then
      hyaw_tmp  = hyam_tmp
      hybw_tmp  = hybm_tmp
    else if (iscoord(varvar_tmp,"ilev")) then
      hyaw_tmp  = hyai_tmp
      hybw_tmp  = hybi_tmp
    end if
    end if

    if (isatt(opt,"time_step") .and. opt@time_step .eq. "avg") then
      tt = -999
    end if

    ;*** set pressure levels for interpolation ***
    if ( isatt(opt, "plev") .and. dimsizes(opt@plev) .eq. 1) then
      plev_work = opt@plev
    else
      print("WARNING: (ignorant_convert_2Ddata) case [tlji2tji_intp] only allow a pressure level or opt@plev is not given. Set to 500 hPa")
      plev_work = 500.
    end if

    ; do interpolation
    varvar_tmp_plev = vinth2p(varvar_tmp, hyaw_tmp, hybw_tmp, plev_work, PS_tmp, 1, P0_tmp/100., 1, False)
    copy_VarAtts(varvar_tmp, varvar_tmp_plev)

    if (isatt(opt,"time_step") .and. opt@time_step .eq. "avg") then
      varvar_return = dim_avg_n_Wrap(varvar_tmp_plev(:,0,:,:),0) ; return (lat,lon) averaged over time
    else
      varvar_return = varvar_tmp_plev(:,0,:,:)  ; return (time,lat,lon)
    end if

    delete(varvar_return@lev_p)
    varvar_return@pressure_level = plev_work+" hPa"
    varvar_return@long_name = varvar_return@long_name+" at "+varvar_return@pressure_level

    option_execute = True
  end if  ; end if of "tlji2tji_intp"

;--------------------------------------------
; convert (time,lev/ilev,lat,lon) to (time,plev,lat)
;   1. interpolate CESM hybrid level to the given pressure levels
;   2. do time & lon average
;--------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tlji2tpj_avgi") then

    hyai_tmp   = filetemp->hyai      ; hybrid A coefficient at layer interface
    hybi_tmp   = filetemp->hybi      ; hybrid B coefficient at layer interface
    hyam_tmp   = filetemp->hyam      ; hybrid A coefficient at layer midpoints
    hybm_tmp   = filetemp->hybm      ; hybrid B coefficient at layer midpoints
    PS_tmp     = filetemp->PS        ; surface pressures  in Pa
    P0_tmp     = filetemp->P0        ; reference pressure in Pa

    if (iscoord(varvar_tmp,"lev")) then
      hyaw_tmp  = hyam_tmp
      hybw_tmp  = hybm_tmp
    else if (iscoord(varvar_tmp,"ilev")) then
      hyaw_tmp  = hyai_tmp
      hybw_tmp  = hybi_tmp
    end if
    end if

    ;*** set pressure levels for interpolation ***
    if ( isatt(opt, "plev") ) then
      plev_work = opt@plev
    else
      plev_work = plev_default
    end if

    ; do interpolation
    varvar_tmp_plev = vinth2p(varvar_tmp, hyaw_tmp, hybw_tmp, plev_work, PS_tmp, 1, P0_tmp/100., 1, False)
    copy_VarAtts(varvar_tmp, varvar_tmp_plev)
    ;printVarSummary(varvar_tmp_plev)

    varvar_tmp_plev_iavg = dim_avg_n_Wrap(varvar_tmp_plev,3)
    ;printVarSummary(varvar_tmp_plev_itavg)

    varvar_return = varvar_tmp_plev_iavg
    option_execute = True

  end if ; end if of "tlji2tpj_avgi"

;--------------------------------------------
; convert (time,lev/ilev,lat,lon) to (plev,lat)
;   1. interpolate CESM hybrid level to the given pressure levels
;   2. do time & lon average
;--------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tlji2pj_avgti") then

    hyai_tmp   = filetemp->hyai      ; hybrid A coefficient at layer interface
    hybi_tmp   = filetemp->hybi      ; hybrid B coefficient at layer interface
    hyam_tmp   = filetemp->hyam      ; hybrid A coefficient at layer midpoints
    hybm_tmp   = filetemp->hybm      ; hybrid B coefficient at layer midpoints
    PS_tmp     = filetemp->PS        ; surface pressures  in Pa
    P0_tmp     = filetemp->P0        ; reference pressure in Pa

    if (iscoord(varvar_tmp,"lev")) then
      hyaw_tmp  = hyam_tmp
      hybw_tmp  = hybm_tmp
    else if (iscoord(varvar_tmp,"ilev")) then
      hyaw_tmp  = hyai_tmp
      hybw_tmp  = hybi_tmp
    end if
    end if

    ;*** set pressure levels for interpolation ***
    if ( isatt(opt, "plev") ) then
      plev_work = opt@plev
    else
      plev_work = plev_default
    end if

    ; do interpolation
    varvar_tmp_plev = vinth2p(varvar_tmp, hyaw_tmp, hybw_tmp, plev_work, PS_tmp, 1, P0_tmp/100., 1, False)
    copy_VarAtts(varvar_tmp, varvar_tmp_plev)
    ;printVarSummary(varvar_tmp_plev)

    varvar_tmp_plev_tavg = dim_avg_n_Wrap(varvar_tmp_plev,0)
    varvar_tmp_plev_itavg = dim_avg_n_Wrap(varvar_tmp_plev_tavg,2)
    ;printVarSummary(varvar_tmp_plev_itavg)

    varvar_return = varvar_tmp_plev_itavg
    option_execute = True

  end if ; end if of "tlji2pj_avgti"

;--------------------------------------------
; convert (time,lev/ilev,lat,lon) to (time,plev,points)
;   1. interpolate CESM hybrid level to the given pressure levels
;   2. do vertical cross section between two given lat/lon points
;--------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tlji2pn_crossji") then

    ;*** set parameters ***
    leftlat   = opt@leftlat
    leftlon   = opt@leftlon 

    rightlat  = opt@rightlat 
    rightlon  = opt@rightlon 

    if (isatt(opt,"npts")) then
      npts = opt@npts
    else
      npts = 100 
    end if

    if ( isatt(opt, "plev") ) then
      plev_work = opt@plev
    else
      plev_work = plev_default
    end if

    hyai_tmp   = filetemp->hyai      ; hybrid A coefficient at layer interface
    hybi_tmp   = filetemp->hybi      ; hybrid B coefficient at layer interface
    hyam_tmp   = filetemp->hyam      ; hybrid A coefficient at layer midpoints
    hybm_tmp   = filetemp->hybm      ; hybrid B coefficient at layer midpoints
    PS_tmp     = filetemp->PS        ; surface pressures  in Pa
    P0_tmp     = filetemp->P0        ; reference pressure in Pa

    if (iscoord(varvar_tmp,"lev")) then
      hyaw_tmp  = hyam_tmp
      hybw_tmp  = hybm_tmp
    else if (iscoord(varvar_tmp,"ilev")) then
      hyaw_tmp  = hyai_tmp
      hybw_tmp  = hybi_tmp
    end if
    end if

    ;*** get lon/lat coordinate ***
    if (isdim(varvar_tmp,"lat")) then
      lat_var = varvar_tmp&lat
    end if

    if (isdim(varvar_tmp,"lon")) then
      lon_var = varvar_tmp&lon
      if (.not.opt@linint2_points_cyclic .and. isMonotonic(lon_var) .eq. 0) then ; lon is not monotonic
        nlon = dimsizes(lon_var)
        do i=0,nlon-1
          if (lon_var(i).gt.180.) then
            lon_var(i) = lon_var(i) - 360.
          end if
        end do
        varvar_tmp&lon = (/lon_var/)
      end if
    end if

    ;*** set pressure levels for interpolation ***
    if ( isatt(opt, "plev") ) then
      plev_work = opt@plev
    else
      plev_work = plev_default
    end if

    ;*** do interpolation ***
    varvar_tmp_lon_work = varvar_tmp
    varvar_tmp_plev = vinth2p(varvar_tmp_lon_work, hyaw_tmp, hybw_tmp, plev_work, PS_tmp, 1, P0_tmp/100., 1, False)
    copy_VarAtts(varvar_tmp_lon_work, varvar_tmp_plev)
    ;printVarSummary(varvar_tmp_plev)

    ;*** use NCL function "gc_latlon" to finds the great circle distance (true surface distance) between two given points ***
    ;*** gc_latlon, https://www.ncl.ucar.edu/Document/Functions/Built-in/gc_latlon.shtml ***
    dist  = gc_latlon(leftlat,leftlon,rightlat,rightlon,npts,2) ; iu=2, return the distance in degrees

    ;*** get all lat/lon between two given points ***
    lon_cross = dist@gclon
    lat_cross = dist@gclat
    ;print(lon_cross+" , "+lat_cross)

    if (.not.opt@linint2_points_cyclic .and. isMonotonic(lon_cross) .eq. 0) then ; lon is not monotonic
      nlon = dimsizes(lon_cross)
      do i=0,nlon-1
        if (lon_cross(i).gt.180.) then
          lon_cross(i) = lon_cross(i) - 360.
        end if
      end do
    end if
    ;print(lon_cross+" , "+lat_cross)

    ;*** use NCL function "linint2_points" to interpolate data between two given points ***
    ;*** linint2_points, https://www.ncl.ucar.edu/Document/Functions/Built-in/linint2_points.shtml ***
    varvar_cross = linint2_points_Wrap(lon_var, lat_var, varvar_tmp_plev, opt@linint2_points_cyclic, lon_cross, lat_cross, 2) ; 2:Reserved for future use

    ;*** get points coordinate ***
    points   = fspan(0.,npts-1,npts)
    varvar_cross@pts = points

    varvar_cross@lat_points = lat_cross
    varvar_cross@lon_points = lon_cross

    if (isatt(opt,"time_step") .and. opt@time_step .eq. "avg") then
      varvar_return = dim_avg_n_Wrap(varvar_cross,0)
    else
      varvar_return = varvar_cross 
    end if

    ;*** set XB coordinates ***
    i0 = 0
    i1 = toint(npts/3)
    i2 = toint(npts/3*2)
    iN = npts-1

    opt@tmXBValues = (/points(i0),points(i1),points(i2),points(iN)/) 
    lat_XB = (/lat_cross(i0),lat_cross(i1), lat_cross(i2), lat_cross(iN)/) 
    lon_XB = (/lon_cross(i0),lon_cross(i1), lon_cross(i2), lon_cross(iN)/) 

    nlon_XB = dimsizes(lon_XB)
    do i=0,nlon_XB-1
      if (lon_XB(i) .gt. 180.) then
        lon_XB(i) = lon_XB(i) - 360.
      end if
    end do
;
    nstring = dimsizes(opt@tmXBValues)
    tmXBLabels = new(nstring, "string")
    do i=0,nstring-1
      tmXBLabels(i) = "lon = "+sprintf ("%3.0f" ,lon_XB(i)) +"~S~o~N~~C~lat = "+sprintf ("%2.0f" , lat_XB(i))+"~S~o"
    end do
    opt@tmXBLabels = tmXBLabels

    ;printVarSummary(varvar_return)
    option_execute = True

  end if ; end if of "tlji2pn_crossji"

;--------------------------------------------
; convert (time,plev,lat) to (plev,lat)
;   1. do time average
;--------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "monthly_series_tlji2ypm_avgji") then

    ntime = dimsizes(varvar_tmp(:,0,0,0))
    if (mod(ntime,12).eq.0) then

      nyears = ntime/12

      ;*** read the variable ***
      hyai_tmp   = filetemp->hyai      ; hybrid A coefficient at layer interface
      hybi_tmp   = filetemp->hybi      ; hybrid B coefficient at layer interface
      hyam_tmp   = filetemp->hyam      ; hybrid A coefficient at layer midpoints
      hybm_tmp   = filetemp->hybm      ; hybrid B coefficient at layer midpoints
      PS_tmp     = filetemp->PS        ; surface pressures  in Pa
      P0_tmp     = filetemp->P0        ; reference pressure in Pa
  
      if (iscoord(varvar_tmp,"lev")) then
        hyaw_tmp  = hyam_tmp
        hybw_tmp  = hybm_tmp
      else if (iscoord(varvar_tmp,"ilev")) then
        hyaw_tmp  = hyai_tmp
        hybw_tmp  = hybi_tmp
      end if
      end if
  
      ;*** set pressure levels for interpolation ***
      if ( isatt(opt, "plev") ) then
        plev_work = opt@plev
      else
        plev_work = plev_default
      end if
  
      ; do interpolation
      varvar_tmp_plev = vinth2p(varvar_tmp, hyaw_tmp, hybw_tmp, plev_work, PS_tmp, 1, P0_tmp/100., 1, False)
      copy_VarAtts(varvar_tmp, varvar_tmp_plev)

      varvar_tmp_plev_jiavg = lat_wgt_avg(varvar_tmp_plev)

      nplev = dimsizes(varvar_tmp_plev_jiavg(0,:))

      varvar_tmp_plev_jiavg_1d = ndtooned(varvar_tmp_plev_jiavg)
      varvar_tmp_ypm0           = onedtond(varvar_tmp_plev_jiavg_1d,(/nyears,12,nplev/))
      copy_VarAtts(varvar_tmp, varvar_tmp_ypm0)

      varvar_tmp_ypm0!0 = "year"
      varvar_tmp_ypm0!1 = "month"
      varvar_tmp_ypm0!2 = "lev_p"

      year_coord  = fspan(0.,tofloat(nyears-1),nyears)
      month_coord = fspan(0.,11.,12)
      lev_p_coord = varvar_tmp_plev&lev_p

      varvar_tmp_ypm0&year  = year_coord
      varvar_tmp_ypm0&month = month_coord
      varvar_tmp_ypm0&lev_p = lev_p_coord

      varvar_tmp_ypm = varvar_tmp_ypm0(year|:,lev_p|:,month|:)

      varvar_return = varvar_tmp_ypm
      ;printVarSummary(varvar_return)

      ;*** check ***
      ;ppp = 36
      ;mmm = 11
      ;print(varvar_tmp_plev_jiavg(mmm::12,ppp)+" , "+varvar_return(:,ppp,mmm))

    else
      print("ERROR: (ignorant_convert_2Ddata) use ["+opt@return_type+"] time steps MUST be 12x times")
      varvar_return = 0.
    end if

    option_execute = True    

  end if ; end if of "monthly_series_tlji2ypm_avgji"

;--------------------------------------------
; convert (year*month,lat,lon) to (year,month,lat)
;   1. do lon average
;--------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "monthly_series_tji2yjm_avgi") then
    ntime = dimsizes(varvar_tmp(:,0,0))

    if (mod(ntime,12).eq.0) then
      aa=0
    else if (ntime .gt. 12) then
      ntime_new = ntime - mod(ntime,12)
      varvar_tmp_new = varvar_tmp(:ntime_new-1,:,:)

      delete(varvar_tmp)
      varvar_tmp = varvar_tmp_new
      ntime = dimsizes(varvar_tmp(:,0,0))

      print("WARNING: (ignorant_convert_2Ddata) use ["+opt@return_type+"] time steps MUST be 12x times")
      print("  only use ["+ntime_new+"] months")

    else
      print("ERROR: (ignorant_convert_2Ddata) use ["+opt@return_type+"] time steps MUST be 12x times")
      varvar_return = 0.
    end if
    end if

      nyears = ntime/12
      nlat = dimsizes(varvar_tmp(0,:,0))

      varvar_tmp_iavg = dim_avg_n_Wrap(varvar_tmp,2)
      varvar_tmp_iavg_1d = ndtooned(varvar_tmp_iavg)
      varvar_tmp_yjm0    = onedtond(varvar_tmp_iavg_1d,(/nyears,12,nlat/))
      copy_VarAtts(varvar_tmp, varvar_tmp_yjm0)

      varvar_tmp_yjm0!0 = "year"
      varvar_tmp_yjm0!1 = "month"
      varvar_tmp_yjm0!2 = "lat"

      year_coord  = fspan(0.,tofloat(nyears-1),nyears)
      month_coord = fspan(0.,11.,12)
      lat_coord = varvar_tmp&lat

      varvar_tmp_yjm0&year  = year_coord
      varvar_tmp_yjm0&month = month_coord
      varvar_tmp_yjm0&lat = lat_coord

      varvar_tmp_yjm = varvar_tmp_yjm0(year|:,lat|:,month|:)

      varvar_return = varvar_tmp_yjm


    ;printVarSummary(varvar_return)

    ;*** check ***
    ;jjj = 10
    ;mmm = 11
    ;print(varvar_tmp_iavg(mmm::12,jjj)+" , "+varvar_return(:,jjj,mmm))

    option_execute = True    

  end if ; end if of "monthly_series_tji2yjm_avgi"

;--------------------------------------------
; convert (year*month,lat,lon) to (year,month,lat,lon)
;--------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "monthly_series_tji2ymji") then
    ntime = dimsizes(varvar_tmp(:,0,0))
    if (mod(ntime,12).eq.0) then

      nyears = ntime/12

      nlat = dimsizes(varvar_tmp(0,:,0))
      nlon = dimsizes(varvar_tmp(0,0,:))

      varvar_tmp_1d = ndtooned(varvar_tmp)
      varvar_tmp_ymji0    = onedtond(varvar_tmp_1d,(/nyears,12,nlat,nlon/))
      copy_VarAtts(varvar_tmp, varvar_tmp_ymji0)

      varvar_tmp_ymji0!0 = "year"
      varvar_tmp_ymji0!1 = "month"
      varvar_tmp_ymji0!2 = "lat"
      varvar_tmp_ymji0!3 = "lon"

      year_coord  = fspan(0.,tofloat(nyears-1),nyears)
      month_coord = fspan(0.,11.,12)
      lat_coord = varvar_tmp&lat
      lon_coord = varvar_tmp&lon

      varvar_tmp_ymji0&year  = year_coord
      varvar_tmp_ymji0&month = month_coord
      varvar_tmp_ymji0&lat = lat_coord
      varvar_tmp_ymji0&lon = lon_coord

      varvar_tmp_ymji = varvar_tmp_ymji0

      varvar_return = varvar_tmp_ymji

    else
      print("ERROR: (ignorant_convert_2Ddata) use ["+opt@return_type+"] time steps MUST be 12x times")
      varvar_return = 0.
    end if

    ;printVarSummary(varvar_return)

    ;*** check ***
    ;jjj = 10
    ;mmm = 11
    ;print(varvar_tmp_iavg(mmm::12,jjj)+" , "+varvar_return(:,jjj,mmm))

    option_execute = True    

  end if ; end if of "monthly_series_tji2ymji"
;--------------------------------------------
; convert (time,plev,lat) to (plev,lat)
;   1. do time average
;--------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tpj2pj_avgt") then

    varvar_temp_tavg = dim_avg_n_Wrap(varvar_temp,0)
    ;printVarSummary(varvar_temp_tavg)

    varvar_return = varvar_temp_tavg
    option_execute = True
  end if ; end if of "tpj2pj_avgt"

;---------
; default
;---------
  if (.not.option_execute) then

    if (isatt(opt, "return_type")) then 
      print("WARNING: (ignorant_convert_2Ddata) fail to execute the case ["+opt@return_type+"]")
      print("  avail case: "+ case_names)
      print("return default variable")
    end if

    varvar_return = varvar_tmp

  end if


;------------------------
; option: change monotonically lon
;------------------------

  option_execute1 = False

  ;*** get longitude array ***
  if ( iscoord(varvar_return,"lon")) then
    lonname = "lon"
    lon = filetemp&$lonname$

    if (isMonotonic(lon) .eq. 0) then ; lon is not monotonic
      option_execute1 = True
    end if

  else if ( iscoord(varvar_return,"longitude")) then
    lonname = "longitude"
    lon = filetemp&$lonname$

    if (isMonotonic(lon) .eq. 0) then ; lon is not monotonic
      option_execute1 = True
    end if

  end if
  end if

  ;*** make lon monotonic ***
  if (option_execute1) then
    nlon = dimsizes(lon)
    lon1 = lon
    do i=0,nlon-1
      if (lon(i).gt.180.) then
        lon1(i) = lon(i) - 360.
      end if
    end do
    varvar_return&$lonname$ = (/lon1/)
  end if

;----------------------------------------
; option: select the returned time step
;----------------------------------------

  if ( isatt(opt, "time_step") .and. isint(opt@time_step)) then
    if ( iscoord(varvar_return,"time") .and. dimsizes(dimsizes(varvar_return)).eq.3) then
      tt = opt@time_step
      varvar_return1 = varvar_return(tt,:,:)
      delete(varvar_return)
      varvar_return = varvar_return1
    ;else
    ; print("ERROR: (ignorant_convert_2Ddata) time coordinate is not found so it cannot be processed")
    ; print("")
    end if
  end if

;----------------------
; return the variable
;----------------------
  res11 = True
  opt_chgatt = True
  opt_chgatt@chg_cncolors = False
  change_var_atts(varvar, varvar_return, res11, opt_chgatt)

  varvar_return@infilename = infilename
  return(varvar_return)
end
;### end of ignorant_convert_2Ddata ###

;#############################################
undef("gsn_csm_contour_siglvl")
function gsn_csm_contour_siglvl(wks:graphic, varvar_any[*][*][*]:numeric, varvar_base[*][*][*]:numeric, res_var[1]:logical, res_sig[1]:logical, opt[1]:logical)
begin

;*************
; initialize
;*************

; get significance level
  if ( opt .and. isatt(opt,"siglvl")  ) then
    siglvl = opt@siglvl
  else
    siglvl = 0.05
  end if

; parameters for ttest
  siglev = 100.*(1.-siglvl)
  iflag= False

;****************
; program start
;****************

; take difference
  varvar_any_tavg = dim_avg_n_Wrap(varvar_any,0)
  varvar_base_tavg = dim_avg_n_Wrap(varvar_base,0)
  varvar_anyMbase_tavg = diff(varvar_any_tavg, varvar_base_tavg)

  ;varvar_anyMbase = diff(varvar_any,varvar_base)
  ;varvar_anyMbase_tavg = dim_avg_n_Wrap(varvar_anyMbase,0)

; res for plot
  res_var_work = res_var
    if (isatt(opt,"cnColors_PosNeg") .and. .not.opt@cnColors_PosNeg ) then
      opt1 = False
    else
      opt1 = True
      copy_VarAtts(opt,opt1)
      cnColors_PosNeg(varvar_anyMbase_tavg,res_var_work,opt1)
    end if
  copy_VarAtts(res_var,res_var_work)

; res for siglvl
;  res_sig_work = res_var
  res_sig_work = True
  res_sig_work@gsnFrame = False
  res_sig_work@gsnDraw = False
  
  res_sig_work@cnLevelSelectionMode = "ExplicitLevels"
  res_sig_work@cnLevels = (/siglev/)

  res_sig_work@cnFillOn = True
  res_sig_work@cnLinesOn = False
  res_sig_work@cnLineLabelsOn = False
  res_sig_work@cnInfoLabelOn = False

  res_sig_work@cnLineThicknessF = 3.

  res_sig_work@cnMonoFillColor = False
  res_sig_work@cnFillColors = (/-1,1/)

  res_sig_work@cnMonoFillPattern = False
  res_sig_work@cnFillPatterns = (/0,4/)   ; slash
  ;res_sig_work@cnFillPatterns = (/0,17/)  ; dotted

  res_sig_work@lbLabelBarOn = False
  res_sig_work@gsnLeftString = ""

  ;copy_VarAtts(res_var,res_sig_work)
  copy_VarAtts(res_sig,res_sig_work)
  ;print(res_sig_work)

  ;if (isatt(res_var,"gsnAddCyclic")) then
  ;  res_sig_work@gsnAddCyclic = res_var@gsnAddCyclic
  ;end if

;************ 
; do t-test
;************
  opt_ttest = True
  copy_VarAtts(res_sig,opt_ttest)
    
  prob = ignorant_ttest(varvar_any, varvar_base, opt_ttest)

  alpha = prob
  alpha = 100.*(1.-prob) ; A significance of 0.05 returned by ttest would yield 95% for alpha. This is often done for plotting.
                         ; ref: http://www.ncl.ucar.edu/Document/Functions/Built-in/ttest.shtml, Example 2

  ;**********************
  ; set plot variables
  ;**********************

  ;print(prob+" , "+alpha)
  ;print(res_sig_work)

  varvar_anyMbase_siglev = alpha

  plot_anyMbase = gsn_csm_contour(wks, varvar_anyMbase_tavg, res_var_work)
  plot_siglvl   = gsn_csm_contour(wks, alpha, res_sig_work)

  overlay(plot_anyMbase,plot_siglvl)

  return(plot_anyMbase)
  ;return(plot_siglvl)
end
;### end of gsn_csm_contour_siglvl ###

;#############################################
undef("gsn_csm_contour_map_siglvl")
function gsn_csm_contour_map_siglvl(wks:graphic, varvar_any0[*][*][*]:numeric, varvar_base0[*][*][*]:numeric, res_var[1]:logical, res_sig[1]:logical, opt[1]:logical)
begin

;*************
; initialize
;*************

;---------------------------------------
; read data_input in certain lat/lon
;---------------------------------------
varvar_any  = get_VarLatLon(varvar_any0 ,res_var)
varvar_base = get_VarLatLon(varvar_base0,res_var)

;-------------
; mean format
;-------------
mean_format = "%7.3f"
mean_prefix = "mean: "
LatWeightCheck = True

if ( isatt(opt, "MeanFormat")) then
  mean_format = opt@MeanFormat
end if
    
if ( isatt(opt, "MeanPrefix")) then
  mean_prefix = opt@MeanPrefix
end if
  
if ( isatt(opt, "LatWeighted")) then
  LatWeightCheck = opt@LatWeighted
;else
;  LatWeightCheck = LatWeightCheck_default
end if

; take difference
  varvar_any_avg = dim_avg_n_Wrap(varvar_any,0)
  varvar_base_avg = dim_avg_n_Wrap(varvar_base,0)
  varvar_anyMbase_avg = diff(varvar_any_avg, varvar_base_avg)
  varvar_anyMbase_tijavg = lat_wgt_avg(varvar_anyMbase_avg)

  ;varvar_anyMbase = diff(varvar_any,varvar_base)
  ;varvar_anyMbase_avg = dim_avg_n_Wrap(varvar_anyMbase,0)
  ;varvar_anyMbase_tijavg = lat_wgt_avg(varvar_anyMbase_avg)

; get lat/lon dimemsion
  nlat = -1
  nlon = -1

  ntime = dimsizes(varvar_base(:,0,0))

  if ( isdim(varvar_base,"lat") ) then
    nlat = dimsizes(varvar_base&lat)
  end if

  if ( isdim(varvar_base,"lon") ) then
    nlon = dimsizes(varvar_base&lon)
  end if

; get significance level
  if ( opt .and. isatt(opt,"siglvl")  ) then
    siglvl = opt@siglvl
  else
    siglvl = 0.05
  end if

; parameters for ttest
  siglev = 100.*(1.-siglvl)
  iflag= False

; res for plot
  res_var_work = res_var
    if (isatt(opt,"cnColors_PosNeg") .and. .not.opt@cnColors_PosNeg ) then
      opt1 = False
    else
      opt1 = True
      copy_VarAtts(opt,opt1)
      cnColors_PosNeg(varvar_anyMbase_avg,res_var_work,opt1)
    end if
  copy_VarAtts(res_var,res_var_work)

; res for siglvl
  res_sig_work = True
  res_sig_work@cnLevelSelectionMode = "ExplicitLevels"
  res_sig_work@cnLevels = (/siglev/)

  res_sig_work@cnFillOn = True
  res_sig_work@cnLinesOn = False
  res_sig_work@cnLineLabelsOn = False
  res_sig_work@cnInfoLabelOn = False

  res_sig_work@cnLineThicknessF = 3.

  res_sig_work@cnMonoFillColor = False
  res_sig_work@cnFillColors = (/-1,1/)

  res_sig_work@cnMonoFillPattern = False
  res_sig_work@cnFillPatterns = (/0,4/)   ; slash
  ;res_sig_work@cnFillPatterns = (/0,17/)  ; dotted

  copy_VarAtts(res_sig,res_sig_work)

  ;if (isatt(res_var,"gsnAddCyclic")) then
  ;  res_sig_work@gsnAddCyclic = res_var@gsnAddCyclic
  ;end if

  if ( .not.isatt(res_var_work, "tiMainString")) then
    res_var_work@tiMainString = "Significance level at "+siglev+"% for "+ntime+" steps"
  end if

  if ( .not.isatt(res_var_work, "gsnCenterString")) then
    res_var_work@gsnCenterString = (mean_prefix + sprintf(mean_format,varvar_anyMbase_tijavg))
  end if

;************ 
; do t-test
;************
  prob = varvar_any(0,:,:)
  prob@long_name = "significance level"
  prob@units = "NONE"

  do j=0,nlat-1
  do i=0,nlon-1
    xtmp = varvar_any(:,j,i)
    ytmp = varvar_base(:,j,i)

  ;*** compute time-average and variance ***
    xAve = dim_avg_Wrap (xtmp)              ; calculate means at each grid point 
    yAve = dim_avg_Wrap (ytmp)
    xVar = dim_variance_Wrap (xtmp)         ; calculate variances
    yVar = dim_variance_Wrap (ytmp)

  ;*** get size ***
    xN = dimsizes(xtmp)
    yN = dimsizes(ytmp)

    prob(j,i) = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False)
  end do
  end do

  alpha = prob
  alpha = 100.*(1.-prob) ; A significance of 0.05 returned by ttest would yield 95% for alpha. This is often done for plotting.
                         ; ref: http://www.ncl.ucar.edu/Document/Functions/Built-in/ttest.shtml, Example 2

  ;when all alpha.eq.siglev with polar projection, significance plot covers the entire plot area, leading to dots/slashes outside the map
  ;To avoid this, manually set two alpha values less than siglev
  if (all(alpha.ge.siglev)) then
    alpha(0,0) = siglev-5.
    alpha(0,1) = siglev-5.
  end if

  ;**********************
  ; set plot variables
  ;**********************

  varvar_anyMbase_siglev = alpha

  ;********
  ; plot
  ;********
  if (isatt(opt,"mapview") .and. opt@mapview .eq. "polar_NH") then
    
    res_var_work@mpProjection = "Stereographic"
    res_var_work@gsnPolar   = "NH"
   
    ;res_sig_work@mpProjection = res_var_work@mpProjection
    ;res_sig_work@gsnPolar     = res_var_work@gsnPolar
    plot_anyMbase_siglev = gsn_csm_contour_map_overlay(wks, varvar_anyMbase_avg, varvar_anyMbase_siglev, res_var_work, res_sig_work)
  
  else if (isatt(opt,"mapview") .and. opt@mapview .eq. "polar_SH") then
    res_var_work@mpProjection = "Stereographic"
    res_var_work@gsnPolar   = "SH"
   
    ;res_sig_work@mpProjection = res_var_work@mpProjection
    ;res_sig_work@gsnPolar     = res_var_work@gsnPolar
    plot_anyMbase_siglev = gsn_csm_contour_map_overlay(wks, varvar_anyMbase_avg, varvar_anyMbase_siglev, res_var_work, res_sig_work)
 
  else
    plot_anyMbase_siglev = gsn_csm_contour_map_overlay(wks, varvar_anyMbase_avg, varvar_anyMbase_siglev, res_var_work, res_sig_work)
  
  end if
  end if

  return(plot_anyMbase_siglev)
end
;### end of gsn_csm_contour_map_siglvl ###

;#############################################
undef("ignorant_convert_1Ddata")
procedure ignorant_convert_1Ddata (infilenames[*]:string, varvars[*]:string, res[1]:logical, opt[1]:logical)
begin

;---------------
; initialization
;---------------
  ;case_name = var(input dimension) + "2" + xx + "_" + yy + process , e.g. ftji2f_t_avgji
  case_names = (/ \
                "tji2tijOneD", \
                "tji2ijOneD_avgt", \
                "tji2tjOneD_avgi", \
                "tji2t_avgji", \
                "vtji2t_vt_avgji", \
                "vtlji2vp_p_avgtji", \
                "vtj2j_vj_avgt", \
                "vtji2j_vj_avgti", \
                "ftlji2fp_p_avgtji", \
                "ftji2t_ft_avgji", \
                "ftj2j_fj_avgt", \
                "ftji2j_fj_avgti", \
                "monthly_series_ftji2m_fm_avgyji", \
                "monthly_series_vtji2m_vm_avgyji", \
                "monthly_series_tji2tjOneD_season_avgi", \
                "monthly_series_tji2tjOneD__avgi", \
                "zzzlkjdflajsldfjlsdjf-no_use" \
                /)
;flag11

  ;*** set res_work ***
  res_work = True
  res_work@xyMonoMarkLineMode = True
  res_work@xyMarkLineMode     = "Lines"
  if (isatt(opt,"xy_colors1")) then
    xy_colors1 = opt@xy_colors1
  else
    xy_colors1 = (/ "red", \
                    "blue",  \
                    "orange", \
;                    "yellow",  \
                    "green", \
                    "purple",  \
                    "pink",  \
                    "brown",  \
                    "gray",  \
                    "purple",  \
                   "black"  \
                  /) 
  end if
  res_work@xyMonoLineColor    = False
  res_work@xyLineColors       = xy_colors1

  res_work@xyMonoDashPattern  = True
  res_work@xyDashPattern      = 0

  res_work@xyMonoLineThickness = True
  res_work@xyLineThicknessF    = 2.5

  res_work@xyMonoMarkerColor = False
  res_work@xyMarkerColors = xy_colors1

  ;*** set opt for change_var_atts ***
  opt_change = True

  ;*** set opt@update_res ***
  if ( isatt(opt, "update_res") .and. .not.opt@update_res ) then
    opt@update_res = False
  else
    opt@update_res = True
  end if

  ;*** set default pressure levels for interpolation ***
  plev_default = (/1., 2., 3., 5., 7., 10., 20., 30., 50., 70., 100., 125., 150., 175., 200., 225., \
      250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 775., 800., 825., \
      850., 875., 900., 925., 950., 975., 1000./)

;------------
; check case
;------------
  ncase  = dimsizes(case_names)
  check1 = False

  do i=0,ncase-1
    if ( isatt(opt, "return_type") .and. opt@return_type .eq. case_names(i)) then
      check1 = True
    end if    
  end do

  if (.not.check1) then
    print("WARNING: (ignorant_convert_1Ddata) fail to execute the case ["+opt@return_type+"]")
    print("  avail case: "+ case_names)
    print("Not return var_xx & var_yy")
  end if

;------------------
; check time step
;------------------
  if (isatt(opt,"time_step_start") .and. isatt(opt,"time_step_end") ) then
    tt1 = opt@time_step_start
    tt2 = opt@time_step_end
  else
    ftemp_name = infilenames(0)
    ftemp1     = addfile(ftemp_name, "r")
    ;printVarSummary(ftemp1)
    if ( isdim(ftemp1,"time") ) then
      ntime = dimsizes(ftemp1&time)
      tt1 = 0
      tt2 = ntime-1
    end if
  end if
  ;print(tt1+" , "+tt2)

;-----------------------------------------------------------------------------------------------------------------
; convert (var,year*month,lat,lon) to xx=(month), yy=(month,var)), vmy=(var,year,month) averaged over lon.lat
;-----------------------------------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "monthly_series_vtji2m_vm_avgyji") then
    varvar_all0 = ignorant_get_vars_file(infilenames,varvars,opt)
    ;printVarSummary(varvar_all0)

    varvar_all1 = lat_wgt_avg(varvar_all0) ; lat-weighted averafe, vars (files,time)
    ;printVarSummary(varvar_all1)
  
    nfile = dimsizes(varvar_all1(:,0))
    ntime = dimsizes(varvar_all1(0,:))

    if (mod(ntime,12).eq.0) then
      nyears = ntime/12

      varvar_all1_1d     = ndtooned(varvar_all1)
      varvar_all1_3d     = onedtond(varvar_all1_1d,(/nfile,nyears,12/))
      varvar_all1_3d!0 = "file"
      varvar_all1_3d!1 = "year"
      varvar_all1_3d!2 = "month"
      copy_VarAtts(varvar_all1,varvar_all1_3d)
      ;printVarSummary(varvar_all1_3d)

      var_xx = fspan(0.,11.,12)
      var_yy = dim_avg_n_Wrap(varvar_all1_3d,1)
      var_vym = varvar_all1_3d

      opt@var_vym = var_vym  ; return to opt
      ;printVarSummary(var_xx)
      ;printVarSummary(var_yy)
      ;printVarSummary(var_fym)

    else
      print("ERROR: (ignorant_convert_1Ddata) use ["+opt@return_type+"] time steps MUST be 12x times")
      var_xx = 0.
      var_yy = 0.
    end if

    unit1 = var_yy@units(0)
    delete(var_yy@units)
    var_yy@units = unit1

    if (isatt(var_yy,"long_name")) then
      delete(var_yy@long_name) ; b.c. both @long_name and @units are a multidimensional assry. gsn_csm_xy would fail.
    end if
  end if 

;-----------------------------------------------------------------------------------------------------------------
; convert (file,year*month,lat,lon) to xx=(month), yy=(month,file)), fmy=(file,year,month) averaged over lon.lat
;-----------------------------------------------------------------------------------------------------------------
;flag11
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "monthly_series_ftji2m_fm_avgyji") then
    varvar_all0 = ignorant_get_var_files(infilenames,varvars,opt)
    ;printVarSummary(varvar_all0)

    varvar_all1 = lat_wgt_avg(varvar_all0) ; lat-weighted averafe, vars (files,time)
    ;printVarSummary(varvar_all1)
  
    nfile = dimsizes(varvar_all1(:,0))
    ntime = dimsizes(varvar_all1(0,:))

    if (mod(ntime,12).eq.0) then
      nyears = ntime/12

      varvar_all1_1d     = ndtooned(varvar_all1)
      varvar_all1_3d     = onedtond(varvar_all1_1d,(/nfile,nyears,12/))
      varvar_all1_3d!0 = "file"
      varvar_all1_3d!1 = "year"
      varvar_all1_3d!2 = "month"
      copy_VarAtts(varvar_all1,varvar_all1_3d)
      ;printVarSummary(varvar_all1_3d)

      var_xx = fspan(0.,11.,12)
      var_yy = dim_avg_n_Wrap(varvar_all1_3d,1)
      var_fym = varvar_all1_3d

      opt@var_fym = var_fym  ; return to opt
      ;printVarSummary(var_xx)
      ;printVarSummary(var_yy)
      ;printVarSummary(var_fym)

    else
      print("ERROR: (ignorant_convert_1Ddata) use ["+opt@return_type+"] time steps MUST be 12x times")
      var_xx = 0.
      var_yy = 0.
    end if
  end if 

;------------------------------------------------------------------------------------
; convert (time,lat,lon) to xx=1D(lat,lon), yy=1D(lat,lon), averaged over time
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tji2tijOneD") then
    var_temp0 = ignorant_get_vars_file(infilenames,varvars,opt)
    var_temp  = var_temp0(0,tt1:tt2,:,:)
    ;printVarSummary(var_temp)

    var_temp_1D = ndtooned(var_temp)
    copy_VarAtts(var_temp,var_temp_1D)
    var_temp_1D@redimension = "convert (lat,lon) to 1D data using ndtooned"
    ;printVarSummary(var_temp_1D)

    var_xx = var_temp_1D
    var_yy = var_temp_1D

    ;atts1 = getvaratts(var_temp_tavg)
    ;natts = dimsizes(atts1)
    ;do i=0,natts-1
    ;  att = atts1(i)
    ;  opt@$att$ = var_temp_tavg@$att$
    ;end do
    ;copy_VarAtts(var_temp_tavg,opt) ; "this would have error message: totype: warning: cannot convert input variable type <float> to type: <logical>" 
  end if

;------------------------------------------------------------------------------------
; convert (time,lat,lon) to xx=1D(lat,lon), yy=1D(lat,lon), averaged over time
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tji2ijOneD_avgt") then
    var_temp0 = ignorant_get_vars_file(infilenames,varvars,opt)
    var_temp  = var_temp0(0,tt1:tt2,:,:)
    var_temp_tavg = dim_avg_n_Wrap(var_temp,0)
    ;printVarSummary(var_temp_tavg)

    var_temp_tavg_1D = ndtooned(var_temp_tavg)
    copy_VarAtts(var_temp_tavg,var_temp_tavg_1D)
    var_temp_tavg_1D@redimension = "convert (lat,lon) to 1D data using ndtooned"
    ;printVarSummary(var_temp_tavg_1D)

    var_xx = var_temp_tavg_1D
    var_yy = var_temp_tavg_1D

    ;atts1 = getvaratts(var_temp_tavg)
    ;natts = dimsizes(atts1)
    ;do i=0,natts-1
    ;  att = atts1(i)
    ;  opt@$att$ = var_temp_tavg@$att$
    ;end do
    ;copy_VarAtts(var_temp_tavg,opt) ; "this would have error message: totype: warning: cannot convert input variable type <float> to type: <logical>" 
  end if

;------------------------------------------------------------------------------------
; convert (time,lat,lon) to xx=1D(time,lat), yy=1D(time,lat), averaged over lon
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tji2tjOneD_avgi") then
    var_temp0 = ignorant_get_vars_file(infilenames,varvars,opt) ; var_temp0 = (var,time,lat,lon)
    ;printVarSummary(var_temp0)
    var_temp  = var_temp0(0,tt1:tt2,:,:)  ; var_temp = (time,lat,lon)
    var_temp_iavg = dim_avg_n_Wrap(var_temp,2)
    ;printVarSummary(var_temp_tavg)

    var_temp_iavg_1D = ndtooned(var_temp_iavg)
    copy_VarAtts(var_temp_iavg,var_temp_iavg_1D)
    var_temp_iavg_1D@redimension = "convert (time,lat) to 1D data using ndtooned"
    ;printVarSummary(var_temp_tavg_1D)

    var_xx = var_temp_iavg_1D
    var_yy = var_temp_iavg_1D

    var_atts = var_xx

    ;atts1 = getvaratts(var_temp_tavg)
    ;natts = dimsizes(atts1)
    ;do i=0,natts-1
    ;  att = atts1(i)
    ;  opt@$att$ = var_temp_tavg@$att$
    ;end do
    ;copy_VarAtts(var_temp_tavg,opt) ; "this would have error message: totype: warning: cannot convert input variable type <float> to type: <logical>" 
  end if

;------------------------------------------------------------------------------------
; convert (time,lat,lon) to xx=1D(time), yy=1D(time), averaged over lon, lat
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "tji2t_avgji") then
    var_temp0 = ignorant_get_vars_file(infilenames,varvars,opt) ; var_temp0 = (var,time,lat,lon)
    var_temp  = var_temp0(0,tt1:tt2,:,:)  ; var_temp = (time,lat,lon)
    var_temp_ijavg = lat_wgt_avg(var_temp)    

    ;printVarSummary(var_temp)
    ;printVarSummary(var_temp_ijavg)
    ;exit

    var_xx = var_temp_ijavg
    var_yy = var_temp_ijavg

    var_atts = var_xx

    ;atts1 = getvaratts(var_temp_tavg)
    ;natts = dimsizes(atts1)
    ;do i=0,natts-1
    ;  att = atts1(i)
    ;  opt@$att$ = var_temp_tavg@$att$
    ;end do
    ;copy_VarAtts(var_temp_tavg,opt) ; "this would have error message: totype: warning: cannot convert input variable type <float> to type: <logical>" 
  end if

;-----------------------------------------------------------------------------------------------------------------
; convert (var,year*month,lat,lon) to xx=1D(month of season,lat), yy=1D(month of season,lat), averaged over lon
;-----------------------------------------------------------------------------------------------------------------

;flag11

  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "monthly_series_tji2tjOneD_season_avgi") then
    var_temp0 = ignorant_get_vars_file(infilenames,varvars,opt) ; var_temp0 = (var,time,lat,lon)
    ;printVarSummary(var_temp0)

    var_temp  = var_temp0(0,tt1:tt2,:,:)  ; var_temp = (time,lat,lon)
    ntime = dimsizes(var_temp(:,0,0))

    if (mod(ntime,12).eq.0 .and. isatt(opt,"season")) then
      if (opt@season.eq."DJF") then
        mm = (/0,1,11/)
      end if
      if (opt@season.eq."MAM") then
        mm = (/2,3,4/)
      end if
      if (opt@season.eq."JJA") then
        mm = (/5,6,7/)
      end if
      if (opt@season.eq."SON") then
        mm = (/8,9,10/)
      end if

      var_temp_iavg = dim_avg_n_Wrap(var_temp,2) ; var_temp_iavg(time,lat)
      ;print(var_temp_iavg)
      ;printVarSummary(var_temp_iavg)

      var_temp_iavg_month1 = var_temp_iavg(mm(0)::12,:)
      var_temp_iavg_month2 = var_temp_iavg(mm(1)::12,:)
      var_temp_iavg_month3 = var_temp_iavg(mm(2)::12,:)

      var_temp_iavg_season = (/var_temp_iavg_month1,var_temp_iavg_month2,var_temp_iavg_month3/)
      var_temp_iavg_season_1D = ndtooned(var_temp_iavg_season)
      copy_VarAtts(var_temp_iavg,var_temp_iavg_season_1D)
      var_temp_iavg_season_1D@redimension = "convert (time,lat) to 1D data using ndtooned"
      ;print(var_temp_iavg_season_1D)

      ;print(opt)
      ;printVarSummary(var_temp_iavg_season_1D)

      var_xx = var_temp_iavg_season_1D
      var_yy = var_temp_iavg_season_1D

      var_atts = var_xx

    else
      print("ERROR: (ignorant_convert_1Ddata) use ["+opt@return_type+"] MUST set opt@season = DJF,MAM,JJA,SON")
      print("ERROR: (ignorant_convert_1Ddata) use ["+opt@return_type+"] time steps MUST be 12x times. current time step ["+ntime+"]")
      var_xx = 0.
      var_yy = 0.
    end if
  end if


;------------------------------------------------------------------------------------
; convert (var,time,lat,lon) to xx=lat, yy=(var,lat), averaged over time,longitude
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "vtji2j_vj_avgti") then
    varvar_all0 = ignorant_get_vars_file(infilenames,varvars,opt)
    varvar_all1 = varvar_all0(:,tt1:tt2,:,:)

    varvar_all1_iavg = dim_avg_n_Wrap(varvar_all1,3)
    varvar_all1_tiavg = dim_avg_n_Wrap(varvar_all1_iavg,1)

    var_yy = varvar_all1_tiavg

    if ( isdim(varvar_all1,"lat")) then
      var_xx = varvar_all1&lat

    else if ( isdim(varvar_all1,"latitude")) then
      var_xx = varvar_all1&latitude

    else
      nlat = dimsizes(varvar_all1(0,:))
      var_xx = fspan(0.,tofloat(nlat-1),nlat)
    end if
    end if

    var_vtj = varvar_all1_iavg
    opt@var_vtj = var_vtj

  end if

;------------------------------------------------------------------------------------
; convert (file,time,lat,lon) to xx=lat, yy=(file,lat), averaged over time,longitude
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "ftji2j_fj_avgti") then
    varvar_all0 = ignorant_get_var_files(infilenames,varvars,opt)
    varvar_all1 = varvar_all0(:,tt1:tt2,:,:)
    varvar_all1_iavg = dim_avg_n_Wrap(varvar_all1,3)
    varvar_all1_tiavg = dim_avg_n_Wrap(varvar_all1_iavg,1)

    var_yy = varvar_all1_tiavg

    if ( isdim(varvar_all1,"lat")) then
      var_xx = varvar_all1&lat

    else if ( isdim(varvar_all1,"latitude")) then
      var_xx = varvar_all1&latitude

    else
      nlat = dimsizes(varvar_all1(0,:))
      var_xx = fspan(0.,tofloat(nlat-1),nlat)
    end if
    end if

    var_ftj = varvar_all1_iavg
    opt@var_ftj = var_ftj
    ;printVarSummary(var_ftj)
  end if

;------------------------------------------------------------------------------------
; convert (file,time,lat) to xx=lat, yy=(file,lat), averaged over time
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "ftj2j_fj_avgt") then
    varvar_all0 = ignorant_get_var_files(infilenames,varvars,opt)
    ;printVarSummary(varvar_all0)
 
    varvar_all00 = varvar_all0(:,tt1:tt2,:)
    varvar_all1 = dim_avg_n_Wrap(varvar_all00,1)

    var_yy = varvar_all1

    if ( isdim(varvar_all1,"lat")) then
      var_xx = varvar_all1&lat

    else if ( isdim(varvar_all1,"latitude")) then
      var_xx = varvar_all1&latitude

    else
      nlat = dimsizes(varvar_all1(0,:))
      var_xx = fspan(0.,tofloat(nlat-1),nlat)
    end if
    end if

    var_ftj = varvar_all00
    opt@var_ftj = var_ftj

  end if  ; end if of "ftj2j_fj_avgt"

;------------------------------------------------------------------------------------
; convert (file,time,lat,lon) to xx=time, yy=(flile,time), averaged over lat and lon
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "ftji2t_ft_avgji") then
    varvar_all0 = ignorant_get_var_files(infilenames,varvars,opt)
    ;printVarSummary(varvar_all0)

    varvar_all1 = lat_wgt_avg(varvar_all0) ; lat-weighted averafe, vars (files,time)
    ;printVarSummary(varvar_all1)
 
    var_yy = varvar_all1(:,tt1:tt2)

    if ( isdim(varvar_all1,"time")) then
      var_xx = varvar_all1&time(tt1:tt2)
    else
      ntime = dimsizes(varvar_all1(0,tt1:tt2))
      var_xx = fspan(0.,tofloat(ntime-1),ntime)
    end if
  end if  ; end if of "ftji2t_ft_avgji"

;------------------------------------------------------------------------------------
; convert (var,time,lat,lon) to xx=time, yy=(var,time), averaged over lat and lon
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "vtji2t_vt_avgji") then
    varvar_all0 = ignorant_get_vars_file(infilenames,varvars,opt)
    ;printVarSummary(varvar_all0)

    varvar_all1 = lat_wgt_avg(varvar_all0) ; lat-weighted averafe, vars (variables,time)
    ;printVarSummary(varvar_all1)
 
    var_yy = varvar_all1(:,tt1:tt2)
    delete(var_yy@long_name)  ; b.c. both @long_name and @units are a multidimensional assry. gsn_csm_xy would fail.

    unit1 = var_yy@units(0)
    delete(var_yy@units)
    var_yy@units = unit1
 
    if ( isdim(varvar_all1,"time")) then
      var_xx = varvar_all1&time(tt1:tt2)
    else
      ntime = dimsizes(varvar_all1(0,tt1:tt2))
      var_xx = fspan(0.,tofloat(ntime-1),ntime)
    end if
  end if  ; end if of "vtji2t_vt_avgji

;------------------------------------------------------------------------------------
; convert (file,lev,time,lat,lon) to xx=(file,pres), yy=(pres), averaged over lat and lon
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "ftlji2fp_p_avgtji") then
    hyai_all   = ignorant_get_var_files(infilenames,"hyai",opt)
    hybi_all   = ignorant_get_var_files(infilenames,"hybi",opt)
    hyam_all   = ignorant_get_var_files(infilenames,"hyam",opt)
    hybm_all   = ignorant_get_var_files(infilenames,"hybm",opt)
    PS_all     = ignorant_get_var_files(infilenames,"PS",opt)
    P0_all     = ignorant_get_var_files(infilenames,"P0",opt)
    varvar_all = ignorant_get_var_files(infilenames,varvars,opt)

    if (iscoord(varvar_all,"lev")) then
      hyaw_all = hyam_all(0,:)
      hybw_all = hybm_all(0,:)
      P0       = P0_all(0,0)
    else if (iscoord(varvar_all,"ilev")) then
      hyaw_all = hyai_all(0,:)
      hybw_all = hybi_all(0,:)
      P0       = P0_all(0,0)
    end if
    end if

    if (isatt(opt,"plev")) then
      plev_work = opt@plev
    else
      plev_work = plev_default
    end if

    ; do interpolation
    varvar_all_plev = vinth2p(varvar_all, hyaw_all, hybw_all, plev_work, PS_all, 1, P0/100., 1, False)
    copy_VarAtts(varvar_all, varvar_all_plev)

    ;*** do lat/lon avg first then time. This is more reasonable ***
    varvar_all_plev_ijavg = lat_wgt_avg(varvar_all_plev)
    ;printVarSummary(varvar_all_plev_ijavg)

    varvar_all_plev_ijavg_tavg = dim_avg_n_Wrap(varvar_all_plev_ijavg(:,tt1:tt2,:),1)
    ;printVarSummary(varvar_all_plev_ijavg_tavg)
    ;*** do lat/lon first then time ***

    ;*** do time avg first then lat/lon ***
    ;varvar_all_plev_tavg = dim_avg_n_Wrap(varvar_all_plev(:,tt1:tt2,:,:,:),1)
    ;printVarSummary(varvar_all_plev_tavg)

    ;varvar_all_plev_tavg_ijavg = lat_wgt_avg(varvar_all_plev_tavg)
    ;;printVarSummary(varvar_all_plev_tavg_ijavg)
    ;*** do time avg first then lat/lon ***

    ;*** these two MUST be equal ***
    ;var111 = varvar_all_plev_ijavg_tavg
    ;var222 = varvar_all_plev_tavg_ijavg
    ;var11m22 = diff(var111,var222)
    ;print("qqq, "+var111+" , "+var222+" , "+var11m22)

    var_xx = varvar_all_plev_ijavg_tavg
    var_yy = plev_work

    var_ftp = varvar_all_plev_ijavg
    opt@var_ftp = var_ftp
    
  end if ; end if of "ftlji2fp_p_avgtji"    

;------------------------------------------------------------------------------------
; convert (var,lev,time,lat,lon) to xx=(var,pres), yy=(pres), averaged over lat and lon
;------------------------------------------------------------------------------------
  if ( isatt(opt, "return_type") .and. opt@return_type .eq. "vtlji2vp_p_avgtji") then
    ftemp1 = addfile(infilenames,"r")
    hyai_all   = ftemp1->hyai
    hybi_all   = ftemp1->hybi
    hyam_all   = ftemp1->hyam
    hybm_all   = ftemp1->hybm
    P0_all     = ftemp1->P0
    varvar_all = ignorant_get_vars_file(infilenames,varvars,opt)

    PS_name = varvars
    PS_name = "PS"
    PS_all = ignorant_get_vars_file(infilenames,PS_name,opt)

    if (iscoord(varvar_all,"lev")) then
      hyaw_all = hyam_all
      hybw_all = hybm_all
      P0       = P0_all
    else if (iscoord(varvar_all,"ilev")) then
      hyaw_all = hyai_all
      hybw_all = hybi_all
      P0       = P0_all
    end if
    end if

    if (isatt(opt,"plev")) then
      plev_work = opt@plev
    else
      plev_work = plev_default
    end if

    ; do interpolation
    varvar_all_plev = vinth2p(varvar_all, hyaw_all, hybw_all, plev_work, PS_all, 1, P0/100., 1, False)
    copy_VarAtts(varvar_all, varvar_all_plev)

    ; do time & lat & lon average
    varvar_all_plev_tavg = dim_avg_n_Wrap(varvar_all_plev,1)
    varvar_all_plev_tavg_ijavg = lat_wgt_avg(varvar_all_plev_tavg)

    var_xx = varvar_all_plev_tavg_ijavg
    var_yy = plev_work

    delete(var_xx@long_name)  ; b.c. both @long_name and @units are a multidimensional assry. gsn_csm_xy would fail.
    unit1 = var_xx@units(0)
    delete(var_xx@units)
    var_xx@units = unit1

  end if ; end if of "ftlji2fp_p_avgtji"    

  ;******************
  ; if var_xx = number
  ;******************
  if ( isatt(opt, "var_xx_number") .and. opt@var_xx_number) then
    nxx = dimsizes(var_xx)
    delete(var_xx)
    var_xx = fspan(0.,tofloat(nxx-1),nxx)
  end if

  ;******************
  ; if var_yy = number
  ;******************
  if ( isatt(opt, "var_yy_number") .and. opt@var_yy_number) then
    nyy = dimsizes(var_yy)
    delete(var_yy)
    var_yy = fspan(0.,tofloat(nyy-1),nyy)
  end if

;-------------
; modify res
;-------------

  if ( isatt(opt, "update_res") .and. opt@update_res ) then

    ;*************************
    ; variable name & unit
    ;*************************
    if (    isatt(opt, "return_type") .and. opt@return_type .eq. "ftji2t_ft_avgji" \
        .or.isatt(opt, "return_type") .and. opt@return_type .eq. "monthly_series_ftji2m_fm_avgyji" \
       ) then
      change_var_atts(varvars,var_yy,res_work,opt_change)
      if (isatt(var_xx,"long_name") .and. isatt(var_xx,"units")) then
        res_work@tiXAxisString = var_xx@long_name+" ("+var_xx@units+")"
      end if
      res_work@gsnRightString = var_yy@units
      res_work@tiYAxisString  = var_yy@units
    end if

    if ( isatt(opt, "return_type") .and. opt@return_type .eq. "ftlji2fp_p_avgtji") then
      res_work@trYReverse = True
      res_work@tiYAxisString = "Pressure (hPa)"

      change_var_atts(varvars,var_xx,res_work,opt_change)
      if (isatt(var_xx,"long_name") .and. isatt(var_xx,"units")) then
        res_work@tiXAxisString = var_xx@long_name+" ("+var_xx@units+")"
      end if
    end if

    if (    isatt(opt, "return_type") .and. opt@return_type .eq. "vtji2t_vt_avgji" \
        .or.isatt(opt, "return_type") .and. opt@return_type .eq. "monthly_series_vtji2m_vm_avgyji" \
       ) then

      if (isatt(var_xx,"long_name") .and. isatt(var_xx,"units")) then
        res_work@tiXAxisString = var_xx@long_name+" ("+var_xx@units+")"
      end if
      change_var_atts("ksjdfkjsdkfjksd",var_yy,res_work,opt_change)  ; just change units so ser varname to a strange string
      if (isatt(var_yy,"units")) then
        res_work@tiYAxisString  = var_yy@units
        res_work@gsnRightString = res_work@tiYAxisString 
      end if

      ;if (isatt(var_yy,"var_name")) then
      ;  string1 = str_join(var_yy@var_name,", ")
      ;  res_work@gsnLeftString = string1
      ;end if
    end if

    if ( isatt(opt, "return_type") .and. opt@return_type .eq. "vtlji2vp_p_avgtji") then
      res_work@trYReverse = True
      res_work@tiYAxisString = "Pressure (hPa)"
      change_var_atts(varvars(0),var_xx,res_work,opt_change)

      if (isatt(var_xx,"units")) then
        res_work@tiXAxisString = var_xx@units
      end if

      if (isatt(var_xx,"var_name")) then
        string1 = str_join(var_xx@var_name,", ")
        res_work@gsnLeftString = string1
      end if
    end if

    if ( isatt(opt, "return_type") .and. opt@return_type .eq. "ftj2j_fj_avgt" \
     .or.isatt(opt, "return_type") .and. opt@return_type .eq. "ftji2j_fj_avgti" ) then
      res_work@tiXAxisString = "Latitude"

      change_var_atts(varvars,var_yy,res_work,opt_change)  ; just change units so ser varname to a strange string
      res_work@tiYAxisString = ""
      ;if (isatt(var_yy,"long_name")) then
      if (isatt(res_work,"gsnLeftString")) then
        res_work@tiYAxisString = res_work@gsnLeftString
        ;res_work@gsnLeftString = var_yy@long_name
      end if
      if (isatt(var_yy,"units")) then
        res_work@tiYAxisString = res_work@tiYAxisString+" ("+var_yy@units+")"
      end if

    end if

    ;**********************************
    ; x_coordinate = "12months_number"
    ;**********************************
    if ( isatt(opt, "x_coordinate") .and. opt@x_coordinate .eq. "12months_number") then
      delete(var_xx)
      var_xx = fspan(0.,11.,12)
      res_work@tmXBMode = "Explicit"
      res_work@tmXBValues = var_xx
      res_work@trXMinF = 0.
      res_work@trXMaxF = 11.
  
      res_work@xyMarkLineMode  = "MarkLines"  ; "Lines" "Markers"
      res_work@xyMonoMarker = True
      res_work@xyMarker = 4                   ; empty circle

      res_work@tiXAxisString = "Month"

      ;res_work@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
      res_work@tmXBLabels = (/"1","2","3","4","5","6","7","8","9","10","11","12"/)
    end if

    ;**********************************
    ; x_coordinate = "12months_1characters"
    ;**********************************
    if ( isatt(opt, "x_coordinate") .and. opt@x_coordinate .eq. "12months_1characters") then
      delete(var_xx)
      var_xx = fspan(0.,11.,12)
      res_work@tmXBMode = "Explicit"
      res_work@tmXBValues = var_xx
      res_work@trXMinF = 0.
      res_work@trXMaxF = 11.
  
      res_work@xyMarkLineMode  = "MarkLines"  ; "Lines" "Markers"
      res_work@xyMonoMarker = True
      res_work@xyMarker = 4                   ; empty circle

      res_work@tiXAxisString = "Month"

      res_work@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
    end if

    ;**********************************
    ; x_coordinate = "12months_3characters"
    ;**********************************
    if ( isatt(opt, "x_coordinate") .and. opt@x_coordinate .eq. "12months_3characters") then
      delete(var_xx)
      var_xx = fspan(0.,11.,12)
      res_work@tmXBMode = "Explicit"
      res_work@tmXBValues = var_xx
      res_work@trXMinF = 0.
      res_work@trXMaxF = 11.
  
      res_work@xyMarkLineMode  = "MarkLines"  ; "Lines" "Markers"
      res_work@xyMonoMarker = True
      res_work@xyMarker = 4                   ; empty circle

      res_work@tiXAxisString = "Month"

      res_work@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
    end if

    ;**********************************
    ; x_coordinate = "MPACE_3hr"
    ;**********************************
    if ( isatt(opt, "x_coordinate") .and. opt@x_coordinate .eq. "MPACE_3hr") then
      delete(var_xx)
      var_xx = fspan(0.,138.,139)
      res_work@tmXBMode = "Explicit"
      tt_all = fspan(8-2./3.,136.-2./3.,17)
      tt_major = tt_all(::2)
      tt_minor = tt_all(1::2)
      a_chrtime = (/"06~C~Oct",    \
                    "08","10","12","14","16","18", \
                    "20"," 22~C~Oct"/)
    
      a_labeltime  = tt_major
      a_labelminor = tt_minor
      res_work@tmXBValues = a_labeltime
      res_work@tmXBLabels = a_chrtime
      res_work@tmXBMinorValues = a_labelminor
    
      res_work@tmXBMode = "Explicit"
      res_work@tmXBMajorThicknessF = 1.2
      res_work@tmXBMinorThicknessF = 1.
    end if

    ;**********************************
    ; x_coordinate = "MPACE_3hr"
    ;**********************************
    if ( isatt(opt, "x_coordinate") .and. opt@x_coordinate .eq. "MPACE_3hr_05OCT14Z_21OCT11Z") then
      delete(var_xx)
      var_xx = fspan(0.,127.,128)
      tt_all = fspan(3.+1./3.,123.+1./3.,16)
      tt_major = tt_all(::2)
      tt_minor = tt_all(1::2)
      a_chrtime = (/"06~C~Oct",    \
                    "08","10","12","14","16","18", \
                    " 20~C~Oct"/)
    
      a_labeltime  = tt_major
      a_labelminor = tt_minor
      res_work@tmXBValues = a_labeltime
      res_work@tmXBLabels = a_chrtime
      res_work@tmXBMinorValues = a_labelminor
    
      res_work@tmXBMode = "Explicit"
      res_work@tmXBMajorThicknessF = 1.2
      res_work@tmXBMinorThicknessF = 1.
    end if

    ;*************
    ; return res
    ;*************
    copy_VarAtts(res_work,res)

  end if  ; end if of update_res

;-----------------
; return variable
;-----------------
  opt@var_xx = var_xx
  opt@var_yy = var_yy

  if (\
            isatt(opt, "return_type") .and. opt@return_type .eq. "ftji2t_ft_avgji" \
       .or. isatt(opt, "return_type") .and. opt@return_type .eq. "vtji2t_vt_avgji" \
       .or. isatt(opt, "return_type") .and. opt@return_type .eq. "vtji2j_vj_avgti" \
       .or. isatt(opt, "return_type") .and. opt@return_type .eq. "ftj2j_fj_avgt" \
       .or. isatt(opt, "return_type") .and. opt@return_type .eq. "ftji2j_fj_avgti" \
       .or. isatt(opt, "return_type") .and. opt@return_type .eq. "monthly_series_ftji2m_fm_avgyji" \
       .or. isatt(opt, "return_type") .and. opt@return_type .eq. "monthly_series_vtji2m_vm_avgyji" \
       .or. isatt(opt, "return_type") .and. opt@return_type .eq. "tji2ijOneD_avgt" \
       .or. isatt(opt, "return_type") .and. opt@return_type .eq. "tji2tijOneD" \
     ) then
     var_atts = var_yy
  end if
  if (\
           isatt(opt, "return_type") .and. opt@return_type .eq. "ftlji2fp_p_avgtji" \
      .or. isatt(opt, "return_type") .and. opt@return_type .eq. "vtlji2vp_p_avgtji" \
     ) then
     var_atts = var_xx
  end if

  atts1 = getvaratts(var_atts)
  natts = dimsizes(atts1)
  do i=0,natts-1
    att = atts1(i)
    opt@$att$ = var_atts@$att$
    varvars@$att$ = var_atts@$att$
  end do

end
;### end of ignorant_convert_1Ddata ###

;#############################################
undef("dimsizes_compare")
function dimsizes_compare (array1, array2)
begin

  varvar_return = True

  varvar_return = True

  dimsizes_1 = dimsizes(array1)
  ndim_1 = dimsizes(dimsizes_1)

  dimsizes_2 = dimsizes(array2)
  ndim_2 = dimsizes(dimsizes_2)

  if (ndim_1 .ne. ndim_2) then
    varvar_return = False
  else
    do i=0,ndim_1-1
      if (dimsizes_1(i) .ne. dimsizes_2(i))
        varvar_return = False
      end if
    end do
  end if

  return(varvar_return)

end
;### end of dimsizes_compare ###

;#############################################
undef("ignorant_sameLB")
procedure ignorant_sameLB(infilenames[*]:string, varvar[1]:string, nclres[1]:logical, opt_sameLB[1]:logical)
begin

if (opt_sameLB) then

  ;*************
  ; initialize
  ;*************
  if (isatt(opt_sameLB,"factor") .and. isfloat(opt_sameLB@factor)) then
    factor = opt_sameLB@factor
  else
    factor = 1.
  end if

  if (isatt(opt_sameLB,"max_levels") .and. isint(opt_sameLB@max_levels)) then
    max_levels = opt_sameLB@max_levels
  else
    max_levels = 15
  end if

  ;filename_base = ""
  if (isatt(opt_sameLB,"filename_base") .and. isstring(opt_sameLB@filename_base)) then
    opt_sameLB@diff2base = True
  end if

  ;*********************************************************
  ; calculate the maximum and minimun in the input variables
  ;*********************************************************
  nfile = dimsizes(infilenames)
  cnmin = 9999999999999.
  cnmax = -9999999999999.

  if (isatt(opt_sameLB,"diff2base") .and. opt_sameLB@diff2base) then
    if (isatt(opt_sameLB,"filename_base") .and. isstring(opt_sameLB@filename_base)) then
      filename_base = opt_sameLB@filename_base
    else
      filename_base = infilenames(nfile-1)
      print("WARNING: (ignorant_sameLB) @filename_base is not set. Set to "+filename_base)
    end if

    opt_base = False
    if (dimsizes(filename_base).eq.1) then
      var_base0 = ignorant_convert_2Ddata(filename_base, varvar, opt_sameLB)
      opt_base = True
    end if
    ;printVarSummary(var_base)

    opt_sameLB@filename_base = filename_base
  end if

  do ff=0,nfile-1

    if (isatt(opt_sameLB,"diff2base") .and. opt_sameLB@diff2base) then
      if (opt_base) then
        var_base = var_base0
      else
        infilename_base = filename_base(ff)
        var_base = ignorant_convert_2Ddata(infilename_base, varvar, opt_sameLB)
      end if
    else
      var_base = 0
    end if
    ;printVarSummary(var_base)

    infilename_any = infilenames(ff)
    var_any = ignorant_convert_2Ddata(infilename_any, varvar, opt_sameLB)
    ;printVarSummary(var_any)

    if (isatt(opt_sameLB,"diff2base") .and. opt_sameLB@diff2base) then
      var_work = diff(var_any,var_base)
    else
      var_work = var_any
    end if

    ;*** get boundaries by Max/Min ***
    ;min1 = min(var_work)
    ;max1 = max(var_work)

    ;*** get boundaries by percentile ***
    stat_var_work = stat_dispersion(var_work,False)
    ;min1 = stat_var_work(3)   ; lower dectile
    ;max1 = stat_var_work(13)  ; upper dectile 
    min1 = stat_var_work(24)  ; lower 5%
    max1 = stat_var_work(25)  ; ipper 5%

    if (min1 .le. cnmin) then
      cnmin = min1
    end if

    if (max1 .ge. cnmax) then
      cnmax = max1
    end if

    delete(var_base)
    delete(infilename_any)
    delete(var_any)
  end do  ; end loop of files
  ;print("min: "+cnmin+", max: "+cnmax)

  ;******************************************
  ; get boundaries of contour levels - by Max/Min values
  ;******************************************

  ;*** modify the max & min
  if (cnmin .lt. 0.) then
    factor_min = factor
  else
    factor_min = 1. + 1.-factor
  end if

  if (cnmax .lt. 0.) then
    factor_max = 1. + 1.-factor
  else
    factor_max = factor
  end if

  cnmin1 = cnmin * factor_min
  cnmax1 = cnmax * factor_max
  if (cnmin1 .gt. cnmax1) then
    cnmin1 = cnmin
    cnmax1 = cnmax
  end if

  ;print("min: "+cnmin+", max: "+cnmax)
  ;print("factor min: "+factor_min+", factor max: "+factor_max)
  ;print("min work: "+cnmin1+", max work: "+cnmax1)

  ;********************** 
  ; get contour levels
  ;********************** 
  mnmxint = nice_mnmxintvl(cnmin1, cnmax1, max_levels, False)
  ;print(mnmxint)
  cnmin1 = mnmxint(0)
  cnmax1 = mnmxint(1)
  cnint1 = mnmxint(2)

  ;*********
  ; return
  ;*********

  nclres@cnLevelSelectionMode = "ManualLevels"
  nclres@cnMinLevelValF = cnmin1
  nclres@cnMaxLevelValF = cnmax1
  nclres@cnLevelSpacingF = cnint1

  opt_sameLB@infilenames   = infilenames
  ;opt_sameLB@filename_base = filename_base

  ;print("min: "+cnmin+", max: "+cnmax)

end if ; end if of opt_sameLB

end
;### end of ignorant_sameLB ###
;#############################################

;#############################################
undef("ignorant_sameScale")
procedure ignorant_sameScale(var_work:numeric, nclres[1]:logical, opt_sameLB[1]:logical)
begin

  ;*************
  ; initialize
  ;*************

  ;*** min/max returned will times factor ***
  if (isatt(opt_sameLB,"factor") .and. isfloat(opt_sameLB@factor)) then
    factor = opt_sameLB@factor
  else
    factor = 1.
  end if

  ;*** set levels for labelbar ***
  if (isatt(opt_sameLB,"max_levels") .and. isint(opt_sameLB@max_levels)) then
    max_levels = opt_sameLB@max_levels
  else
    max_levels = 15
  end if

  ;*** set option_scale ***    
  if (isatt(opt_sameLB,"OptionScale") .and. opt_sameLB@OptionScale .eq. "XScale") then
    option_scale = "XScale"
    if (.not.isatt(opt_sameLB,"percentile")) then
      opt_sameLB@percentile = 0
    end if

  else if (isatt(opt_sameLB,"OptionScale") .and. opt_sameLB@OptionScale .eq. "YScale") then
    option_scale = "YScale"
    if (.not.isatt(opt_sameLB,"percentile")) then
      opt_sameLB@percentile = 0
    end if

  else
    option_scale = "LBScale"
  end if
  end if

  ;*** Some stats values in stat_dispersion will set to FillValue if size < 1000
  ndim = dimsizes(ndtooned(var_work))
  ;print("ndim= "+ndim)
  if (ndim.lt.1000) then
    if (option_scale.eq."XScale" .or. option_scale.eq."YScale") then
      opt_sameLB@percentile = 0
    end if
    if (option_scale.eq."LBScale") then
      opt_sameLB@percentile = 10
    end if
  end if

  ;**************************
  ; get low/high boundaries
  ;**************************

  ;http://www.ncl.ucar.edu/Document/Functions/Contributed/stat_dispersion.shtml
  stat_var_work = stat_dispersion(var_work,False)

  if (isatt(opt_sameLB,"percentile") .and. opt_sameLB@percentile.eq.10) then
    cnmin = stat_var_work(3)   ; lower 10%
    cnmax = stat_var_work(13)  ; upper 10%

  else if (isatt(opt_sameLB,"percentile") .and. opt_sameLB@percentile.eq.1) then
    cnmin = stat_var_work(23)  ; lower 1%
    cnmax = stat_var_work(26)  ; ipper 1%

  else if (isatt(opt_sameLB,"percentile") .and. opt_sameLB@percentile.eq.0) then
    cnmin = stat_var_work(2)   ; minimum
    cnmax = stat_var_work(14)  ; maximum

  else
    cnmin = stat_var_work(24)  ; lower 5%
    cnmax = stat_var_work(25)  ; ipper 5%
  end if
  end if
  end if

  ;print(stat_var_work)
  ;print("min  = "+min(var_work)+", max  = "+max(var_work))
  ;print("cnmin= "+cnmin        +", cnmax= "+cnmax)

  ;*** modify the max & min
  if (cnmin .lt. 0.) then
    factor_min = factor
  else
    factor_min = 1. + 1.-factor
  end if

  if (cnmax .lt. 0.) then
    factor_max = 1. + 1.-factor
  else
    factor_max = factor
  end if

  cnmin1 = cnmin * factor_min
  cnmax1 = cnmax * factor_max
  if (cnmin1 .gt. cnmax1) then
    cnmin1 = cnmin
    cnmax1 = cnmax
  end if

  ;************************************* 
  ; get Scale and return - X, Y, or LB
  ;************************************* 

  ;*** XScale ***
  if (option_scale .eq. "XScale") then
    nclres@trXMinF = cnmin1
    nclres@trXMaxF = cnmax1
  end if

  ;*** YScale ***
  if (option_scale .eq. "YScale") then
    nclres@trYMinF = cnmin1
    nclres@trYMaxF = cnmax1
  end if

  ;*** LBScale ***
  if (option_scale .eq. "LBScale") then
    mnmxint = nice_mnmxintvl(cnmin1, cnmax1, max_levels, False)
    ;print(mnmxint)
    cnmin1 = mnmxint(0)
    cnmax1 = mnmxint(1)
    cnint1 = mnmxint(2)
  
    nclres@cnLevelSelectionMode = "ManualLevels"
    nclres@cnMinLevelValF = cnmin1
    nclres@cnMaxLevelValF = cnmax1
    nclres@cnLevelSpacingF = cnint1
  
    ;opt_sameLB@cnColors_PosNeg = True
    if (isatt(opt_sameLB,"cnColors_PosNeg") .and. opt_sameLB@cnColors_PosNeg) then
      opt111 = opt_sameLB
      opt111@cnlevels_predefine = True
      copy_VarAtts(nclres, opt111)
      ;print(opt111)
      cnColors_PosNeg(var_work, nclres, opt111)
    end if
  end if  ; end if of LBScale

end
;### end of ignorant_sameScale ###
;#############################################

;#############################################
undef("trim_fig")
procedure trim_fig(wks_type[1]:string, wks_name[1]:string, opt1[1]:logical)
begin

  if (opt1) then
    if (isatt(opt1,"fig_name")) then
      fig_name = opt1@fig_name
    else
      fig_type = wks_type
      fig_name = wks_name+"."+wks_type
    end if
    system("convert -trim "+fig_name+" "+fig_name)
    ;system("convert -trim "+fig_name+" "+fig_name+" && echo 'Done. trim figures'")
  end if
end
;### end of trim_fig ###
;#############################################

;#############################################
undef("convert_files_vars")
function convert_files_vars(filenames[*]:string, varnames[*]:string, opt[1]:logical)
begin
 
;----------------
; program start
;----------------
  opt_work_1D = False
  opt_work_2D = False
  res = True

  nfile = dimsizes(filenames)
  nvar  = dimsizes(varnames)

  names_atts = (/"long_name","units"/)  ; these attributes will be copied to the return variable
  natts = dimsizes(names_atts)

  string_atts      = new((/nvar, natts/), "string") ; attributes

  ;-----------------------
  ; get case name for 1D
  ;-----------------------
  case_names_1D = (/ \
                "tji2tijOneD", \
                "tji2ijOneD_avgt", \
                "tji2tjOneD_avgi", \
                "tji2t_avgji", \
                "vtji2t_vt_avgji", \
                "vtlji2vp_p_avgtji", \
                "vtj2j_vj_avgt", \
                "vtji2j_vj_avgti", \
                "ftlji2fp_p_avgtji", \
                "ftji2t_ft_avgji", \
                "ftj2j_fj_avgt", \
                "ftji2j_fj_avgti", \
                "monthly_series_ftji2m_fm_avgyji", \
                "monthly_series_vtji2m_vm_avgyji", \
                "monthly_series_tji2tjOneD_season_avgi", \
                "monthly_series_tji2tjOneD__avgi", \
                "zzzlkjdflajsldfjlsdjf-no_use" \
                /)
  n1Dcase = dimsizes(case_names_1D)

  do i=0,n1Dcase-1
    if ( opt@return_type .eq. case_names_1D(i)) then
      opt_work_1D = True
      copy_VarAtts(opt, opt_work_1D)
    end if
  end do

  ;-----------------------
  ; get case name for 2D
  ;-----------------------
  case_names_2D = (/ \
                "tlji2pt_avgji", \
                "tlji2tji_intp", \
                "tlji2ji_intp_avgt", \
                "tlji2pn_crossji", \
                "tlji2pj_avgti", \
                "tlji2tpj_avgi", \
                "tji2jt_avgi", \
                "tji2ji_avgt", \
                "tpj2pj_avgt", \
                "monthly_series_tlji2ypm_avgji", \
                "monthly_series_tji2yjm_avgi", \
                "monthly_series_tji2ymji", \
                "sdlkfjsldjflsjdflsdjfl" \
              /)
  n2Dcase = dimsizes(case_names_2D)

  do i=0,n2Dcase-1
    if ( opt@return_type .eq. case_names_2D(i)) then
      opt_work_2D = True
      copy_VarAtts(opt, opt_work_2D)
    end if
  end do

  ;--------------
  ; if a 2D case
  ;--------------
  if (opt_work_2D) then
    filename1 = filenames(0)
    varname1  = varnames(0)
    var1 = ignorant_convert_2Ddata(filename1,varname1,opt)
    dim = dimsizes(var1)  ; get the return dimensions

    files_vars_2D = new( (/nfile, nvar,dim(0), dim(1)/), typeof(var1))
   
    do ff=0, nfile-1
    do vv=0, nvar-1
      ;*** get the variable ***
      filename1 = filenames(ff)
      varname1  = varnames(vv)
      files_vars_2D(ff,vv,:,:) = ignorant_convert_2Ddata(filename1,varname1,opt)

      ;*** get attributes ***
      do aa=0,natts-1
        att1 = names_atts(aa)
        if (isatt(var1, att1)) then
          string_atts(vv,aa) = var1@$att1$
        end if
      end do ; end do of natts
    end do   ; end do of vars
    end do   ; end do of files

    ;*** delete unused attributes ***
    if (isatt(files_vars_2D,"infilename")) then
      delete(files_vars_2D@infilename)
    end if

  end if  ; end if of opt_work_2D

  ;--------------
  ; if a 1D case
  ;--------------
  if (opt_work_1D) then

    ;*********************
    ; know return x or y
    ;*********************
    coord = ""

    if (isatt(opt, "coord")) then
      coord = opt@coord
    end if

    ;*** read a variables in the first two files ***
    if (nfile.ge.2 .and. coord.eq."") then
      opt1 = opt_work_1D
      filename1 = filenames(0)

      opt2 = opt_work_1D
      filename2 = filenames(1)

      varname1  = varnames(0)

      ignorant_convert_1Ddata(filename1,varname1,res,opt1)
      ignorant_convert_1Ddata(filename2,varname1,res,opt2)

      if (all(opt1@var_xx.eq.opt2@var_xx)) then
        coord = "xx"
        ;printVarSummary(opt1)
        ;printVarSummary(opt2)
      else
        coord = "yy"
      end if 
    end if

    ;*** read the first two variabless in a file ***
    if (nvar.ge.2 .and. coord.eq."") then
      opt1 = opt_work_1D
      filename1 = filenames(0)

      varname1  = varnames(0)
      varname2  = varnames(1)

      ignorant_convert_1Ddata(filename1,varname1,res,opt1)
      ignorant_convert_1Ddata(filename1,varname2,res,opt2)

      if (all(opt1@var_xx.eq.opt2@var_xx)) then
        coord = "xx"
        ;printVarSummary(opt1)
        ;printVarSummary(opt2)
      else
        coord = "yy"
      end if 
    end if

    ;***********************
    ; set the return array
    ;***********************
    opt111 = opt_work_1D
    filename1 = filenames(0)
    varname1  = varnames(0)
    ignorant_convert_1Ddata(filename1,varname1,res,opt111)
    if (coord.eq."xx") then
      coord_values  = opt111@var_xx
      type = typeof(opt111@var_yy)
    else if (coord.eq."yy") then
      coord_values  = opt111@var_yy
      type = typeof(opt111@var_xx)
    else
      coord_values  = opt111@var_xx
      type = typeof(opt111@var_yy)
    end if
    end if

    ncoord = dimsizes(coord_values)

    files_vars_1D = new( (/nfile, nvar, ncoord/), type)
    files_vars_1D!2 = "coord"
    files_vars_1D&coord = coord_values
    ;printVarSummary(files_vars_1D)

    ;***********************
    ; read and save variables
    ;***********************
    delete(opt111)
    do ff=0, nfile-1
    do vv=0, nvar-1
      opt111 = opt_work_1D
      filename1 = filenames(ff)
      varname1  = varnames(vv)
      ignorant_convert_1Ddata(filename1,varname1,res,opt111)

      if (coord.eq."xx") then
        files_vars_1D(ff,vv,:) = opt111@var_yy
      else if (coord.eq."yy") then
        files_vars_1D(ff,vv,:) = opt111@var_xx
      else
        files_vars_1D(ff,vv,:) = opt111@var_yy
      end if
      end if

      do aa=0,natts-1
        att1 = names_atts(aa)
        if (isatt(opt111, att1)) then
          string_atts(vv,aa) = opt111@$att1$
        end if
      end do
      delete(opt111)
    end do
    end do

  end if  ; end if of opt_work_1D

;------------------------------
; set attributes & filenames
;------------------------------
  if (opt_work_1D) then
    varvar_return = files_vars_1D
  else if (opt_work_2D) then
    varvar_return = files_vars_2D
  else
    print("ERROR: ")
    varvar_return = 0
  end if
  end if

  varvar_return!0 = "files"
  varvar_return!1 = "vars"
 
  varvar_return@return_type = opt@return_type

  if (isatt(opt, "write2file") .and. opt@write2file) then

    do aa=0,natts-1
      att1 = names_atts(aa)
      if (isatt(varvar_return, att1)) then
        delete(varvar_return@$att1$)
      end if

      varvar_return@$att1$ = string_1Dto0D(string_atts(:,aa))
      ;end if
    end do
    varvar_return@varnames    = string_1Dto0D(varnames)
    varvar_return@infilenames = string_1Dto0D(filenames)

  else

    varvar_return&files = filenames
    varvar_return&vars  = varnames

    do aa=0,natts-1
      att1 = names_atts(aa)
      if (isatt(varvar_return, att1)) then
        delete(varvar_return@$att1$)
      end if

      varvar_return@$att1$ = string_atts(:,aa)
      ;end if
    end do
    varvar_return@varnames    = varnames
    varvar_return@infilenames = filenames
  end if

  return(varvar_return)
  ;printVarSummary(varvar_return)

end
;### end of convert_files_vars ###
;#############################################

;#############################################
undef("ignorant_get_vars_files")
function ignorant_get_vars_files(filenames[*]:string, varnames[*]:string, opt[1]:logical)
begin

;----------------
; program start 
;----------------

  opt111 = opt
  nfile = dimsizes(filenames)
  nvar  = dimsizes(varnames)

  ;***************************
  ; get variable dimensions
  ;***************************
  filename1 = filenames(0)
  varname1  = varnames(0)

  option_var = ""

  if (isatt(opt111,"var_choice") .and. opt111@var_choice.eq."2Ddata") then
    var01 = ignorant_convert_2Ddata(filename1, varnames, opt111)
    option_var = "2Ddata"

  ;if (isatt(opt111,"var_choice" .and. opt111@var_choice.eq."1Ddata")) then
  else  
    var01 = ignorant_get_vars_file(filename1, varnames, opt111)
    option_var = "normal"
  end if

  var_dimsizes = dimsizes(var01) ;( tmp_infile, read_var_name )  
  var_dimnames = getvardims(var01) ;(tmp_infile, read_var_name)
  var_dimnames_sizes = dimsizes(var_dimnames)

  ;***************************
  ; using ignorant_2Ddata
  ;***************************
  if (option_var .eq. "2Ddata")

    ;*** create return variable ***
    if ( var_dimnames_sizes .eq. 2 ) then
      n1 = var_dimsizes(0)  
      n2 = var_dimsizes(1)
      work_dimsizes = (/nfile,nvar,n1,n2/) 
    end if    
    if ( var_dimnames_sizes .eq. 3 ) then
      n1 = var_dimsizes(0)  
      n2 = var_dimsizes(1)
      n3 = var_dimsizes(2)
      work_dimsizes = (/nfile,nvar,n1,n2,n3/) 
    end if    
  
    work_var = new( work_dimsizes, typeof(var01))

    ;*** read variables ***
    do ff=0, nfile-1
    do vv=0, nvar-1
      filename1 = filenames(ff)
      varname1  = varnames(vv)
      tmp_var = ignorant_convert_2Ddata(filename1, varname1, opt111)

      if ( var_dimnames_sizes .eq. 2 ) then
        work_var(ff,vv,:,:) = tmp_var
      end if
  
      if ( var_dimnames_sizes .eq. 3 ) then
        work_var(ff,vv,:,:,:) = tmp_var
      end if
      delete(tmp_var)
    end do
    end do
  end if ; end if of option_var = 2Ddata

  ;***************************
  ; using ignorant_get_vars_file
  ;***************************
  if (option_var .eq. "normal")
    ;*** define outpuy array size ***
    aa = 1 + var_dimnames_sizes
    work_dimsizes = new( aa, "integer" )
  
    work_dimsizes(0) = nfile
    do i = 1, aa-1
      work_dimsizes(i) = var_dimsizes(i-1)
    end do

    work_var = new( work_dimsizes, typeof(var01))
    ;*** read variables ***

    do ff=0, nfile-1
      filename1 = filenames(ff)
      tmp_var = ignorant_get_vars_file(filename1, varnames, opt111)
  
      if ( var_dimnames_sizes .eq. 1 ) then
        work_var(ff,:) = tmp_var
      end if
      if ( var_dimnames_sizes .eq. 2 ) then
        work_var(ff,:,:) = tmp_var
      end if
  
      if ( var_dimnames_sizes .eq. 3 ) then
        work_var(ff,:,:,:) = tmp_var
      end if
  
      if ( var_dimnames_sizes .eq. 4 ) then
        work_var(ff,:,:,:,:) = tmp_var
      end if
    end do
  end if ; end if of option_var=normal

  ;*** set return variable ***
  work_var!0 = "files"
  work_var&files = filenames
  work_var@filenames = filenames

  work_var!1 = "variable"
  work_var&variable = varnames
  ;printVarSummary(work_var)  

  return(work_var)

end
;### end of ignorant_get_vars_files ###
;#############################################

;#############################################
undef("regline_print")
function regline_print(x[*]:numeric, y[*]:numeric, string_input[1]:string, opt[1]:logical)
begin

  ;*** determine attributes ****
  if (isatt(x,"long_name")) then
    VARx_name = x@long_name
  else
    VARx_name = "??"
  end if

  if (isatt(y,"long_name")) then
    VARy_name = y@long_name
  else
    VARy_name = "??"
  end if

  if (isatt(x,"units")) then
    VARx_units = x@units
  else
    VARx_units = "??"
  end if

  if (isatt(y,"units")) then
    VARy_units = y@units
  else
    VARy_units = "??"
  end if

  if (isatt(opt,"format")) then
    format1 = opt@format
  else
    format1 = "%3.2f"
  end if

  if (isatt(opt,"printout")) then
    printout = opt@printout
  else
    printout = True
  end if

  nx = dimsizes(x)

  ; use NCL regline_stats, http://www.ncl.ucar.edu/Document/Functions/Contributed/regline_stats.shtml
  stats = regline_stats(x,y)

  ; compute the mean of slope and 95% interval
  bmean = 0.5*(stats@b95(1)+stats@b95(0))
  bint = 0.5*(stats@b95(1)-stats@b95(0))

  yint = 0.5*(stats@y95(1)-stats@y95(0))

  ;*** set string to print on screen ***
  string_print1 = string_input
  string_print2 = "#points = "+nx+", VARx= "+VARx_name+" ("+VARx_units+"), VARy= "+VARy_name+" ("+VARy_units+"), Y=bx+y"

  ;*** print ***
  if (printout) then
    print("")
    print("----------------------------------")
    print(string_print1+"")
    print(string_print2+"")
    print("----------------------------------")
    print("")
    print("b: "+stats@b(1))
    print("b95 interval: "+bint)
    print("b95: "+stats@b95)
    print("")
    print("y: "+stats@b(0))
    print("y95 interval: "+yint)
    print("y95: "+stats@y95)
    print("")
  end if

  ;*** return user-define variable ***
  stats@slope = stats@b(1)
  stats@slope_95IC = bint
  stats@intercept = stats@b(0)
  stats@intercept_95IC = yint

  y_fit = stats@slope*x + stats@intercept
  stats@y_fit = y_fit
  ;print(x+", "+y_fit)

  stats@line_xx = (/min(x),max(x)/)
  stats@line_yy = (/min(y_fit),max(y_fit)/)

  string_slope = sprintf(format1, stats@slope)
  string_slope_95IC = sprintf(format1, bint)
  string_intercept = sprintf(format1, abs(stats@intercept))
  string_intercept_95IC = sprintf(format1, yint)

  if (stats@intercept.lt.0.) then
    string_sign = " - "
  else
    string_sign = " + "
  end if

  string_r2 = sprintf(format1, stats@r2)
  string_fit = "y = "+ \
               string_slope+"(~F34~1~F21~"+string_slope_95IC+")*x" + string_sign + \
               string_intercept+"(~F34~1~F21~"+string_intercept_95IC+")"+  \
               ", r~S~2~N~= "+string_r2
  stats@string_fit = string_fit

  return(stats)

end
;### end of regline_print ##
;#############################################

;#############################################
undef("overlay_siglvl_dots")
procedure overlay_siglvl_dots(wks[1]:graphic, plot[1]:graphic, xx:numeric, yy:numeric, pp:numeric, res[1]:logical, opt_plot[1]:logical)
begin
    ;opt_plot@siglvl = 0.5  ; significance level. p-value less than this will be plotted
    ;opt_plot@Draw = True   ; draw and frame the plot

;****************
; set plot res
;****************
  if (.not.isatt(opt_plot, "siglvl")) then
    opt_plot@siglvl = 0.5  ; significance level. p-value less than this will be plotted
  end if

  xy_colors = (/ "red", \
                 "blue",  \
                 "orange", \
                 "yellow",  \
                 "green", \
                 "purple",  \
                 "pink",  \
                 "brown",  \
                 "gray",  \
                 "purple",  \
                 "black"  \
              /)  ; set colors in ignorant_xy

  polyres = True
  polyres@gsMarkerIndex = 16   ; solid dot
  polyres@gsMarkerSizeF = 5.
  ;polyres@gsMarkerColor = 1

  if (isatt(res, "xyLineColors")) then
    ;polyres@gsMarkerColor = res@xyLineColors
    xy_colors1 = res@xyLineColors
    option_xycolors = True
  else
   xy_colors1 = xy_colors
  end if

  copy_VarAtts(res, polyres)
  ;print(polyres)

  if (isatt(polyres, "gsnDraw")) then
    delete(polyres@gsnDraw)
  end if

  if (isatt(polyres, "gsnFrame")) then
    delete(polyres@gsnFrame)
  end if

;************************
; check dimension sizes
;************************

  nxx = dimsizes(dimsizes(xx))
  nyy = dimsizes(dimsizes(yy))
  npp = dimsizes(dimsizes(pp))
  ;print("nxx: "+nxx+", nyy: "+nyy+", npp: "+npp)
 
  if (nxx.gt.2 .or. nyy.gt.2 .or. npp.gt.2) then
    print("ERROR: (zzzz) only support 1D or 2D array")
    exit
  end if

  ;*** xx, yy, pp have same dimensions ***
  if (npp.eq.nxx .and. nxx.eq.nyy) then
    pp_work = pp
    xx_work = xx
    yy_work = yy

  ;*** xx, pp have 2 dimensions but yy has 1 dimension ***
  else if (npp.eq.2 .and. nxx.eq.npp .and. nyy.ne.npp) then
    npp1 = dimsizes(pp(:,0))
    npp2 = dimsizes(pp(0,:))
    nxx1 = dimsizes(xx(:,0))
    nxx2 = dimsizes(xx(0,:))

    if (npp1.ne.nxx1 .or. npp2.ne.nxx2) then
      print("ERROR: dimension is not consistent")
      print("pp[*,*]: "+npp1+","+npp2+" xx[*,*]: "+nxx1+","+nxx2)
      exit
    end if

    pp_work = pp
    xx_work = xx

    yy_work = xx
    yy_work = 0.
    do i=0,npp1-1
      yy_work(i,:) = yy
    end do

  ;*** yy, pp have 2 dimensions but xx has 1 dimension ***
  else if (npp.eq.2 .and. nxx.ne.npp .and. nyy.eq.npp) then
    npp1 = dimsizes(pp(:,0))
    npp2 = dimsizes(pp(0,:))
    nyy1 = dimsizes(yy(:,0))
    nyy2 = dimsizes(yy(0,:))

    if (npp1.ne.nyy1 .or. npp2.ne.nyy2) then
      print("ERROR: dimension is not consistent")
      print("pp[*,*]: "+npp1+","+npp2+" yy[*,*]: "+nyy1+","+nyy2)
      exit
    end if
    pp_work = pp
    yy_work = yy

    xx_work = yy
    xx_work = 0.
    do i=0,npp1-1
      xx_work(i,:) = xx
    end do
    
  else
    print("ERROR: (zzzz) dimension of input are NOT consistent")
    print("x: "+nxx+", y: "+nyy+", p-value: "+npp)
    exit
  end if
  end if
  end if

;***************************
; plot significant symbols
;***************************
  plot_polymarker = new(1000, graphic)  ; id for each p-values (polymarker)
  iip = 0  ; id number of each polymarker

  ;*** input dimension is 1 ***
  if (npp.eq.1) then
    npp1 = dimsizes(pp_work)

    do i=0,npp1-1
      polyres1 = polyres 
      if (.not.isatt(polyres1, "gsMarkerColor")) then
        polyres1@gsMarkerColor = xy_colors1(0)
      end if

      if (.not.ismissing(pp_work(i)) .and. pp_work(i) .le. opt_plot@siglvl) then
        plot_polymarker(iip) = gsn_add_polymarker(wks, plot, \
                                       xx_work(i),yy_work(i),polyres1)
        iip = iip+1
      end if
      delete(polyres1)
    end do
  end if

  ;*** input dimension is 2 ***
  if (npp.eq.2) then
    npp1 = dimsizes(pp_work(:,0))
    npp2 = dimsizes(pp_work(0,:))
    ;print(npp1+","+npp2)

    do i=0,npp1-1
    do j=0,npp2-1
      polyres1 = polyres 
      if (.not.isatt(polyres1, "gsMarkerColor")) then
        polyres1@gsMarkerColor = xy_colors1(i)
      end if

      if (.not.ismissing(pp_work(i,j)) .and. pp_work(i,j) .le. opt_plot@siglvl) then
        plot_polymarker(iip) = gsn_add_polymarker(wks, plot, \
                                       xx_work(i,j),yy_work(i,j),polyres1)
        iip = iip+1
      end if
      delete(polyres1)
    end do
    end do
  end if

;************
; draw plot
;************
  if (isatt(opt_plot, "Draw") .and. opt_plot@Draw) then
    draw(plot)
    frame(wks)
  end if

end
;### end of overlay_siglvl_dots ###
;#############################################


;#############################################
;undef("")
;function ()
;begin
;
;end
;### end of ###
;#############################################

