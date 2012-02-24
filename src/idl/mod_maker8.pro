; mod_maker7.pro plus subroutines
;*************************************************************
; Interpolates in time, lat/longitude into the NCEP re-analyses to
; produce a model file for a particular site for a series of dates.
;
; INPUTS:
;  mod_maker.input   Contains the names of the NetCDF files  
;                     (only first 5 lines are read)
;  ~/ggg/ncdf/Site_xxxxxxxxxx_AT.nc    Air Temperature
;  ~/ggg/ncdf/Site_xxxxxxxxxx_GH.nc    Geopotential Height
;  ~/ggg/ncdf/Site_xxxxxxxxxx_TP.nc    Tropopause pressure
;  ~/ggg/ncdf/Site_xxxxxxxxxx_SH.nc    Specific Humidity

; OUTPUTS:
;  ~/ggg/models/gnd/??yyyymmdd.mod  for all dates present in NCEP files
;
;
; Requires that the user has acquired the four necessary .nc files
; (Air Temperature, Geopotential Height, Tropopause Pressure, Specific Humidity)
; from the following four NOAA website:
; http://www.cdc.noaa.gov/cgi-bin/db_search/DBSearch.pl?&Dataset=NCEP+Reanalysis+Pressure+Level&Dataset=NCEP+Reanalysis+Daily+Averages+&Pressure+Level&Variable=Air+Temperature
; http://www.cdc.noaa.gov/cgi-bin/db_search/DBSearch.pl?&Dataset=NCEP+Reanalysis+Pressure+Level&Dataset=NCEP+Reanalysis+Daily+Averages+&Pressure+Level&Variable=Specific+Humidity
; http://www.cdc.noaa.gov/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis+Tropopause+Level&Dataset=NCEP+Reanalysis+Daily+Averages+Tropopause+Level&Variable=Pressure
; http://www.cdc.noaa.gov/cgi-bin/db_search/DBSearch.pl?&Dataset=NCEP+Reanalysis+Pressure+Level&Dataset=NCEP+Reanalysis+Daily+Averages+Pressure+Level&Variable=Geopotential+height
;
;  IMPORTANT: Enter longitude in deg E (eg, Park Falls: 240-280 E)
;             Use 6-hourly data (4 times daily)
;             Select all available pressure levels (1000 - 10 mb)
;             Select the time period
;             Select "Create subset without making a plot"
; 
; FTP the four resulting files from the NOAA website to ~/ggg/ncdf/ 
; giving them descriptive names, e.g.,
; ~/ggg/ncdf/NCEP_ParkFalls_40N_50N_265E_275E_20040501_20061130_6hourly_AT.nc
; ~/ggg/ncdf/NCEP_ParkFalls_40N_50N_265E_275E_20040501_20061130_6hourly_GH.nc
; ~/ggg/ncdf/NCEP_ParkFalls_40N_50N_265E_275E_20040501_20061130_6hourly_TP.nc
; ~/ggg/ncdf/NCEP_ParkFalls_40N_50N_265E_275E_20040501_20061130_6hourly_SH.nc
;
;  You must then edit mod_maker.input and enter at the top of the file
;  the lat/long of your site and the names of the four input files.
;  
;  You are now ready to run mod_maker8
;
;*************************************************************
;    MOD_MAKER7.PRO     13 Oct 2010
; Had to change the way that VARGET was called 
; since varid can no longer be hardwired to the variable #
; So, for example, instead of:
;    NCDF_VARGET, ncid, 1, Lat_AT       ; Read in variable 'lat'
; we now have:
;    varid=ncdf_varid(ncid,'lat')
;    NCDF_VARGET, ncid, varid, Lat_AT       ; Read in variable 'lat'
; Thanks to Hisako for reporting and fixing this problem.
;
;************************************************************
;    MOD_MAKER6.PRO     17 Feb 2009
; We were having memory problems when running with
; an entire year of global data. The larger input
; files were over 500MB each and so several GB of
; memory was needed by mod_maker5.  So I reduced the
; memory footprint by a factor ~4 by reorganizing
; the main program to re-use global data arrays
; (Global_Data) rather than having separate
; arrays for each parameter (AT, GH, TP, SH).
; This required that the interpolation to the lat/long
; of the site happen immediately after reading in the
; global data. Also removed some of the intermediate
; global arrays (e.g. Geopotential_Hgt, Specific_Humidity)
;*************************************************************
  
; Major changes introduced by GCT in Jan/Mar 2007:
;
; 1) Got rid of the need for inlog files for selecting the list
;  of dates for which to create models. MOD_MAKER now produces
;  a .mod file for every day that is covered by the NCEP files,
;  irrespective of whether there were spectra measured or not.
;  Of course, this will lead to models being created for days
;  having no spectra, but who cares? Previously, everyone using
;  mod_maker was creating inlog files containing every day of the
;  year anyway, because it was too tedious to delete the lines
;  representing completely cloudy days. Now mod_maker, uses the
;  IDL CALDAT intrinsic function to generate the dates of
;  the models from the julian day numbers in the NCEP files.
;
; 2) Replaced hard-wired values for the NCEP model latitude,
;  longitude, and time grid-steps with the actual values inside
;  the NCEP files. This provides flexibility to use models with
;  different lat/long grids or, more importantly, different
;  time-steps. For example, mod_maker can now use the 6-hourly
;  or the daily NCDF models as input. 
;
; 3) Replaced the use of the jplsummr.mod file for altitudes
;  above 10 mbar, with the US Standard atmosphere at levels
;  of 5, 2, 1, 0.1, 0.01, 0.001, 0.0001 mbar. This greatly
;  reduces the number of model levels above 10 mbar from
;  ~120 previously to 7 now.  For ground-based spectra, the
;  coarseness of the new models above 10 mbar should make
;  a negligible difference to the retrievals.
; 
; 4) Implemented piece-wise linear time-interpolation.
;  It now interpolates linearly in time between the NetCDF
;  models provided to obtain the T/P/H values for local noon.
;  (to within +-16 minutes).
;  For Park Falls this doesn't make much improvement because
;  local noon is very close to 18UT, which is a standard
;  NCEP model time.  But for Darwin (131E) local noon is
;  03:16 UT, which is almost half way between OO UT and 06 UT.
;
; 5) Implemented a more sophisticated scheme for calculating
;  the H2O vmr above 300 mbar (the top of the NCEP H2O profiles).
;  Previously, the H2O vmr was immediately dropped to the
;  stratospheric H2O vmr at the next level (250 mbar). This was
;  okay for mid- and high-latitude profiles when the tropopause
;  level was 250 mbar or more. But for low latitude sites (e.g.
;  Darwin) where the tropopause pressure is ~100 mbar this old
;  approach introduced a large discontinuity in the H2O profile.
;  Now, the H2O is reduced in a more gradual fashion until it
;  reaches the stratospheric value.
;
; 6) Implemented a more sophisticated scheme for dealing with
;  bad H2O vmrs in the NCEP file. Previously, -ve H2O vmrs
;  were replaced with the previous +ve value. But the NCEP
;  files contain many instances of H2O_vmr=0.0, or H2O_vmr<1ppm,
;  which were not caught previously. Now, all vmrs less than
;  3 ppm are replaced by SQRT(3ppm * Previous value).
;
; 7) Changed nearly all the variable names to make them more
;  consistent with each other. All variables with the suffix: 
;  AT denote Air Temperature
;  GH denote Geometric Height
;  TP denote Tropopause Pressure
;  SH denote Specific Humidity (mass mixing ratio of H2O).
;
; 8) Fixed a bug in the old mod_maker by which the index
;  of the site latitude and longitude tropopause pressure
;  (e.g. Site_Lat_IndexP) was tested againts the index limit
;  of the geometric height file (nH_Lat or nH_lon), i.e.
;  if((Site_Lat_IndexP gt nH_Lat-1) or $
;     (Site_Lat_IndexP lt 0) or $
;     (Site_Lon_IndexP gt nH_Lon-1) or $
;     (Site_Lon_IndexP lt 0)) then begin
; Fortunately, since all the NCEP files have the same Lat/Lon grid,
; this mistake didn't matter.
;
; 9) All input now comes from the file mod_maker.input.
;   This includes the lat/long of the site and the names
;   of the files containing the NCDF data. Thus there is
;   no longer any user input required for mod_maker6.
;   This facilitates automation, and improves tracability
;   of the results.
;
;10) Broke out the functionality of writing the .mod files into
;  a subroutine (write_mod.pro). This simplifies mod_maker.pro.
;  Wrote linear_interpN and bilinear_interpN subroutines to
;  handle the time and lat/long interpolations, respectively.
;
;
; Note: IDL is an interpretive language, so subroutines must
; be defined *before* they are called by the main program.
; So the subroutines are found at the beginning of the file
; and then the main program is found near the end.

; Subroutine svp_wv_over_ice
; Uses the Goff-Gratch equation to calculate the saturation vapor
; pressure of water vapor over ice at a user-specified temperature.
; Input:  temp (K)
; Output: svp (mbar)
pro svp_wv_over_ice, temp, svp
  t0=273.16   ; triple point temperature (K)
  tr=t0/temp
  yy=-9.09718*(tr-1)-3.56654*alog10(tr)+0.876793*(1-1/tr)
  svp=6.1173*10^yy   ; saturation vapor pressure over ice (mbar)
end

; SUBROUTINE WRITE_MOD:
pro write_mod, mod_path, Site_Lat, Lev_AT, sat, sgh, stp, ssh
; Creates a GGG-format .mod file when called by mod_maker.pro.
; INPUTS:
;     Site_Lat    ; The latitude of the site
;     Lev_AT      ; The pressure levels on which the data are tabulated
;     sat         ; Site Noon Atmospheric Temperature profile (vector)
;     sgh         ; Site Noon Geometric Height profile in km (vector)
;     ssh         ; Site Noon Specific Humidity profile (vector)
;     stp         ; Site Noon Tropopause Pressure (scalar)

;*************************************************************
; Define US Standard Atmosphere (USSA) for use above 10 mbar
  p_ussa=[10.0,  5.0,   2.0,   1.0,   0.1,   0.01,  0.001, 0.0001]
  t_ussa=[227.7, 239.2, 257.9, 270.6, 231.6, 198.0, 189.8, 235.0]
  z_ussa=[31.1,  36.8,  42.4,  47.8,  64.9,  79.3,  92.0,  106.3]


    print,Mod_Path,stp/100
; Export the head of .mod
    openw,lunw,Mod_Path,/get_lun
    printf,lunw,'4  5'
    printf,lunw,6378.137,6.000E-05,Site_Lat,9.81, $
      sgh[0],1013.25, stp/10^2, $
      format='(f8.3,1x,e11.4,1x,f7.3,1(1x,f5.3),3(1x,f8.3))'
    printf,lunw, ' mbar        Kelvin         km      g/mole      vmr'
    printf,lunw, 'Pressure  Temperature     Height     MMW        H2O'
  
;; Determine tropopause altitude and hence pressure.
;; Search for the first instance of lapse_rate > -2 K/km
;; above Level 6 (500 mbar). Will not find a trop below 500 mbar.
;    lr_thresh=-2.0  ;  Lapse Rate Threshold
;    k=6    ; Start at layer 6 (500 mbar) to avoid any PBL T-inversions
;    lapse_was=(sat(k-1)-sat(k-2))/(sgh(k-1)-sgh(k-2))  ; GCT 2010-02-13
;    lapse_rate=(sat(k)-sat(k-1))/(sgh(k)-sgh(k-1))
;    zlrwas=0.5*(sgh(k-1)+sgh(k-2))  ; GCT 2010-02-11
;    zlr=0.5*(sgh(k)+sgh(k-1))
;    while ( lapse_rate lt lr_thresh ) do begin
;       k=k+1
;       lapse_was=lapse_rate
;       zlrwas=zlr
;       lapse_rate=(sat(k)-sat(k-1))/(sgh(k)-sgh(k-1))
;       zlr=0.5*(sgh(k)+sgh(k-1)) 
;    endwhile
;    ztrop=zlrwas+(zlr-zlrwas)*(lr_thresh-lapse_was)/(lapse_rate-lapse_was)
;    sh=(sgh(k)-sgh(k-2))/alog(Lev_AT(k-2)/Lev_AT(k))
;    ptrop=Lev_AT(k-1)*exp(-(ztrop-sgh(k-1))/sh)
;    print,ztrop,ptrop,stp/100

; Export the Pressure, Temp and SHum for lower levels (1000 to 300 mbar)
    for k=0,N_elements(SSH)-1 do begin
      if (ssh[k] lt 4.0e-06 and k gt 0) then begin
        print, 'Replacing insufficient SHum ',mod_path, Lev_AT[k],SSH[k]
        SSH[k]=sqrt(4.0e-06*SSH[k-1])
      endif
      printf,lunw,Lev_AT[k], sat[k], sgh[k],28.9640,ssh[k], $
        format='(e9.3,4x,f7.3,4x,f7.3,4x,f7.4,4x,1e9.3)'
    endfor

; Export Pressure and Temp for middle levels (250 to 10 mbar)
;  which have no SHum reanalysis.
    dsh=SSH[N_elements(SSH)-1]
    lstp=alog10(stp)
    for k=N_elements(SSH),n_elements(Lev_AT)-1 do begin
        zz=alog10(Lev_AT[k])  ; log10[pressure]
        strat_h2o=7.5E-06*exp(-3.0*(zz/lstp)^2)
        dsh=dsh*(Lev_AT[k]/Lev_AT[k-1])^(5.5-Lev_AT[k]/100)
        svp_wv_over_ice, sat[k], svp
        satvmr=svp/Lev_AT[k]   ; vmr of saturated WV at T/P
;        print,Lev_AT[k],stp/100,sat[k],svp,dsh,dsh/satvmr
        if ( dsh gt satvmr ) then dsh=satvmr
        if (Lev_AT(k) lt stp/100) then dsh=sqrt(dsh*strat_h2o)  ; above the trop
      printf,lunw,Lev_AT[k], sat[k], $
        sgh[k],28.9640,max([dsh,strat_h2o]), $
        format='(e9.3,4x,f7.3,4x,f7.3,4x,f7.4,4x,1e9.3)'
    endfor

; Get the difference between the USSA and given site temperature at 10 mbar,
    Delta_T=sat[16]-t_ussa[0]

; Export the P-T profile above 10mbar
    for k=1,7 do begin
        Delta_T=Delta_T/2
        zz=alog10(p_ussa[k])  ; log10[pressure]
        strat_h2o=7.5E-06*exp(-0.25*zz^2)
;        print, p_ussa[k],strat_h2o
      printf,lunw,p_ussa[k],t_ussa[k]+Delta_T,z_ussa[k],28.9640,strat_h2o, $
        format='(e9.3,4x,f7.3,4x,f7.3,4x,f7.4,4x,1e9.3)'
    endfor
    free_lun,lunw
end

; SUBROUTINE LINEAR_INTERP0
; Evaluates  fout = fin(xx) for fout of dimension 0 (scalar)
pro linear_interp0, fin, xx, fout
   index_xx = fix(xx)
   if index_xx lt 0 then index_xx = 0
   fr_xx=xx-index_xx
   fout=(1-fr_xx)*fin[index_xx]+fr_xx*fin[index_xx+1]
end

; SUBROUTINE LINEAR_INTERP1
; Evaluates  fout = fin(xx) for fout of dimension 1 (vector)
pro linear_interp1, fin, xx, fout
   index_xx = fix(xx)
   if index_xx lt 0 then index_xx = 0
   fr_xx=xx-index_xx
   fout=(1-fr_xx)*fin[*,index_xx]+fr_xx*fin[*,index_xx+1]
end

; SUBROUTINE BILINEAR_INTERP1
; Evaluates  fout = fin(xx,yy) for fout of dimension 1 (vector)
pro bilinear_interp1, fin, xx, yy, fout
  index_xx=fix(xx)
  fr_xx=xx-index_xx
  index_yy=fix(yy)
  fr_yy=yy-index_yy
  fout= $
    (fin[index_xx,index_yy,*]*(1.0-fr_xx)+ $
     fin[index_xx+1,index_yy,*]*fr_xx)*(1.0-fr_yy)+ $
    (fin[index_xx,index_yy+1,*]*(1.0-fr_xx)+ $
     fin[index_xx+1,index_yy+1,*]*fr_xx)*fr_yy 
   fout=reform(fout,/overwrite)
end

; SUBROUTINE BILINEAR_INTERP2
; Evaluates  fout = fin(xx,yy) for fout of dimension 2 (array)
pro bilinear_interp2, fin, xx, yy, fout
  index_xx=fix(xx)
  fr_xx=xx-index_xx
  index_yy=fix(yy)
  fr_yy=yy-index_yy
  fout= $
    (fin[index_xx,index_yy,*,*]*(1.0-fr_xx)+ $
     fin[index_xx+1,index_yy,*,*]*fr_xx)*(1.0-fr_yy)+ $
    (fin[index_xx,index_yy+1,*,*]*(1.0-fr_xx)+ $
     fin[index_xx+1,index_yy+1,*,*]*fr_xx)*fr_yy 
  fout=reform(fout,/overwrite)
end

; MAIN PROGRAM:
pro mod_maker8
;*************************************************************
; Read site abbreviation, lat/long & names of NCDF files from "mod_maker.input"
; Note that it only reads the first 6 lines of this file.
; So you can keep other stuff further down for re-use later.
;*************************************************************
 print,'modmaker8.0  09-Sep-2011  GCT'
 Site_Abbrev='xx'
 ncdf_AT_file=' '
 ncdf_GH_file=' '
 ncdf_TP_file=' '
 ncdf_SH_file=' '
 home_path=string(getenv('GGGPATH'))
 openr, unitr,home_path+'/src/idl/mod_maker.input',/get_lun
 readf, unitr,Site_Abbrev
 readf, unitr,Site_Lat, Site_Lon
 readf, unitr,ncdf_AT_file
 readf, unitr,ncdf_GH_file
 readf, unitr,ncdf_TP_file
 readf, unitr,ncdf_SH_file
 close, unitr
 free_lun, unitr
if(Site_Lon lt 0.0) then Site_Lon=Site_Lon+360.0
if(Site_Lon gt 180.0 ) then Site_Lon_180=Site_Lon-360.0 else Site_Lon_180=Site_Lon
;*************************************************************
; Read Air Temperature file and interpolate to Lat/Long of site:
;*************************************************************
  ncid = NCDF_OPEN(ncdf_AT_file)       ; Open The NetCDF file
  varid=ncdf_varid(ncid,'level')
  NCDF_VARGET, ncid,  varid, Lev_AT       ; Read in variable 'level'
  varid=ncdf_varid(ncid,'lat')
  NCDF_VARGET, ncid,  varid, Lat_AT       ; Read in variable 'lat'
  varid=ncdf_varid(ncid,'lon')
  NCDF_VARGET, ncid,  varid, Lon_AT       ; Read in variable 'lon'
  varid=ncdf_varid(ncid,'time')
  NCDF_VARGET, ncid,  varid, Tim_AT       ; Read in variable 'time'
  varid=ncdf_varid(ncid,'air')
  NCDF_VARGET, ncid,  varid, Global_Data  ; Read in variable 'Global_Data'
  NCDF_ATTGET, ncid,  varid, 'add_offset', Global_Data_add_offset
  NCDF_ATTGET, ncid,  varid, 'scale_factor', Global_Data_scale_factor
  NCDF_CLOSE, ncid      ; Close the NetCDF file
  Global_Data = Global_Data*Float(Global_Data_scale_factor) + Float(Global_Data_add_offset)

  dlon=Lon_AT[N_Elements(Lon_AT)-1]-Lon_AT[0]
 if (Site_Lat-Lat_AT[N_Elements(Lat_AT)-1])*(Site_Lat-Lat_AT[0]) gt 0 or $
    (Site_Lon-Lon_AT[N_Elements(Lon_AT)-1])*(Site_Lon-Lon_AT[0])*dlon gt 0 then begin
       print,'The data file '+ncdf_AT_file+' does not cover this Site'
       print,'Latitudes: ',Lat_AT
       print,'Longitudes:',Lon_AT
       stop
 endif
  bilinear_interp2,Global_Data,(Site_Lon-Lon_AT[0])/(Lon_AT[1]-Lon_AT[0]),(Site_Lat-Lat_AT[0])/(Lat_AT[1]-Lat_AT[0]), Site_AT

;*************************************************************
; Read Geopotential Height file and interpolate to Lat/Long of site:
;*************************************************************
  ncid = NCDF_OPEN(ncdf_GH_file)       ; Open The NetCDF file
  varid=ncdf_varid(ncid,'level')
  NCDF_VARGET, ncid,  varid, Lev_GH       ; Read in variable 'level'
  varid=ncdf_varid(ncid,'lat')
  NCDF_VARGET, ncid,  varid, Lat_GH       ; Read in variable 'lat'
  varid=ncdf_varid(ncid,'lon')
  NCDF_VARGET, ncid,  varid, Lon_GH       ; Read in variable 'lon'
  varid=ncdf_varid(ncid,'time')
  NCDF_VARGET, ncid,  varid, Tim_GH       ; Read in variable 'time'
  varid=ncdf_varid(ncid,'hgt')
  NCDF_VARGET, ncid,  varid, Global_Data          ; Read in variable 'Global_Data'
  NCDF_ATTGET, ncid,  varid, 'add_offset', Global_Data_add_offset
  NCDF_ATTGET, ncid,  varid, 'scale_factor', Global_Data_scale_factor
  NCDF_CLOSE, ncid     ; Close the NetCDF file

  Global_Data = Global_Data*Float(Global_Data_scale_factor) + Float(Global_Data_add_offset)

  Earth_Radius=6378.137*1000.0

; Convert Geopotential height (H) to Geometric Height (Z).
; The latter is always larger than the former (i.e., Z > H)
;  Global_Data = Global_Data /(1-Global_Data/Earth_Radius)
;  z = h/(1-h/r)
;  h = z/(1+z/r) 

 dlon=Lon_GH[N_Elements(Lon_GH)-1]-Lon_GH[0]
 if (Site_Lat-Lat_GH[N_Elements(Lat_GH)-1])*(Site_Lat-Lat_GH[0]) gt 0 or $
    (Site_Lon-lon_GH[N_Elements(Lon_GH)-1])*(Site_Lon-Lon_GH[0])*dlon gt 0 then begin
       print,'The data file '+ncdf_GH_file+' does not cover this Site'
       print,'Latitudes: ',Lat_GH
       print,'Longitudes:',Lon_GH
       stop
 endif
  bilinear_interp2,Global_Data,(Site_Lon-Lon_GH[0])/(Lon_GH[1]-Lon_GH[0]),(Site_Lat-Lat_GH[0])/(Lat_GH[1]-Lat_GH[0]), Site_GH

;*************************************************************
; Read Tropopause Pressure file and interpolate to Lat/Long of site:
;*************************************************************
  ncid = NCDF_OPEN(ncdf_TP_file)       ; Open The NetCDF file
  varid=ncdf_varid(ncid,'lat')
  NCDF_VARGET, ncid,  varid, Lat_TP       ; Read in variable 'lat'
  varid=ncdf_varid(ncid,'lon')
  NCDF_VARGET, ncid,  varid, Lon_TP       ; Read in variable 'lon'
  varid=ncdf_varid(ncid,'time')
  NCDF_VARGET, ncid,  varid, Tim_TP       ; Read in variable 'time'
  varid=ncdf_varid(ncid,'pres')
  NCDF_VARGET, ncid,  varid, Global_Data         ; Read in variable 'Global_Data'
  NCDF_ATTGET, ncid,  varid, 'add_offset', Global_Data_add_offset
  NCDF_ATTGET, ncid,  varid, 'scale_factor', Global_Data_scale_factor
  NCDF_CLOSE, ncid     ; Close the NetCDF file

  Global_Data = Global_Data*Float(Global_Data_scale_factor) + Float(Global_Data_add_offset)

 dlon=Lon_TP[N_Elements(Lon_TP)-1]-Lon_TP[0]
 if (Site_Lat-Lat_TP[N_Elements(Lat_TP)-1])*(Site_Lat-Lat_TP[0]) gt 0 or $
    (Site_Lon-Lon_TP[N_Elements(Lon_TP)-1])*(Site_Lon-Lon_TP[0])*dlon gt 0 then begin
       print,'The data file '+ncdf_TP_file+' does not cover this Site'
       print,'Latitudes: ',Lat_TP
       print,'Longitudes:',Lon_TP
       stop
 endif

  bilinear_interp1,Global_Data,(Site_Lon-Lon_TP[0])/(Lon_TP[1]-Lon_TP[0]),(Site_Lat-Lat_TP[0])/(Lat_TP[1]-Lat_TP[0]), Site_TP
;*************************************************************
; Read Specific Humidity file and interpolate to Lat/Long of site:
;*************************************************************
  ncid = NCDF_OPEN(ncdf_SH_file)
  varid=ncdf_varid(ncid,'level')
  NCDF_VARGET, ncid,  varid, Lev_SH       ; Read in variable 'level'
  varid=ncdf_varid(ncid,'lat')
  NCDF_VARGET, ncid,  varid, Lat_SH       ; Read in variable 'lat'
  varid=ncdf_varid(ncid,'lon')
  NCDF_VARGET, ncid,  varid, Lon_SH       ; Read in variable 'lon'
  varid=ncdf_varid(ncid,'time')
  NCDF_VARGET, ncid,  varid, TIM_SH       ; Read in variable 'tim'
  varid=ncdf_varid(ncid,'shum')
  NCDF_VARGET, ncid,  varid, Global_Data         ; Read in variable 'Global_Data'
  NCDF_ATTGET, ncid, varid, 'add_offset', Global_Data_add_offset
  NCDF_ATTGET, ncid, varid, 'scale_factor', Global_Data_scale_factor
  NCDF_CLOSE, ncid     ; Close the NetCDF file

  Global_Data = (Global_Data*Float(Global_Data_scale_factor)+Float(Global_Data_add_offset))*28.964/18.02

 dlon=Lon_SH[N_Elements(Lon_SH)-1]-Lon_SH[0]
 if (Site_Lat-Lat_SH[N_Elements(Lat_SH)-1])*(Site_Lat-Lat_SH[0]) gt 0 or $
    (Site_Lon-Lon_SH[N_Elements(Lon_SH)-1])*(Site_Lon-Lon_SH[0])*dlon gt 0 then begin
       print,'The data file '+ncdf_SH_file+' does not cover this Site'
       print,'Latitudes: ',Lat_SH
       print,'Longitudes:',Lon_SH
       stop
 endif

  bilinear_interp2,Global_Data,(Site_Lon-Lon_SH[0])/(Lon_SH[1]-Lon_SH[0]),(Site_Lat-Lat_SH[0])/(Lat_SH[1]-Lat_SH[0]), Site_SH

;*************************************************************
; Check that the full number of levels have been downloaded
; Geopotential & Air Temperature files
;*************************************************************
  If(N_Elements(Lev_AT) lt 17) Then Begin
    print,'You must download all 17 levels of AT data from the'
    print,'NCEP/NCAR website. Please redo your AT download'
  EndIf
  If(N_Elements(Lev_GH) lt 17) Then Begin
    print,'You must download all 17 levels of GH data from the'
    print,'NCEP/NCAR website. Please redo your GH download'
  EndIf
  If(N_Elements(Lev_SH) lt 8) Then Begin
    print,'You must download all 8 levels of SH data from the'
    print,'NCEP/NCAR website. Please redo your SH download'
  EndIf

;*************************************************************************
; Define the name of the models
; Naming convention: 14 characters
;   1-2: Site name e.g. 'ca'-Caltech, 'pk'-Park Falls, 'tm'-Table Mountain
;   3-10: year+month+day, in the form yyyymmdd
;   11-14: suffix '.mod'
;*************************************************************************
  julian=Tim_AT[0]/24.0+julday(1,1,1,0)
  caldat,julian,month,day,year
  Site_Noon_Tim=(julday(month,day,year,12)-julday(1,1,1,0))*24.0 $
    -Site_Lon_180/15

; Main Loop over model dates
  while ( Site_Noon_Tim le Tim_AT[N_Elements(Tim_AT)-1] $
    and   Site_Noon_Tim le Tim_GH[N_Elements(Tim_GH)-1] $
    and   Site_Noon_Tim le Tim_TP[N_Elements(Tim_TP)-1] $
    and   Site_Noon_Tim le Tim_SH[N_Elements(Tim_SH)-1] ) do begin
  
; Interpolate in time to local noon
    linear_interp1,Site_AT,(Site_Noon_Tim-Tim_AT[0])/(Tim_AT[1]-Tim_AT[0]),Site_Noon_AT
    linear_interp1,Site_GH,(Site_Noon_Tim-Tim_GH[0])/(Tim_GH[1]-Tim_GH[0]),Site_Noon_GH
    linear_interp0,Site_TP,(Site_Noon_Tim-Tim_TP[0])/(Tim_TP[1]-Tim_TP[0]),Site_Noon_TP
    linear_interp1,Site_SH,(Site_Noon_Tim-Tim_SH[0])/(Tim_SH[1]-Tim_SH[0]),Site_Noon_SH
    Site_Noon_GH=Site_Noon_GH/1000   ; convert m to km

;  Write the .mod file 
    caldat,julian,month,day,year
    Mod_Name=string(Site_Abbrev,year,month,day,'.mod',FORMAT='(A2,I4,I2.2,I2.2,A4)')
    write_mod,home_path+'/models/gnd/'+mod_name, Site_Lat,Lev_AT,Site_Noon_AT, $
    Site_Noon_GH,Site_Noon_TP,Site_Noon_SH
    julian=julian+1
    Site_Noon_Tim=Site_Noon_Tim+24
  endwhile
end
