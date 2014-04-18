; mod_maker9.pro plus subroutines
;*************************************************************
; Interpolates in time, lat/longitude into the NCEP re-analyses to
; produce a model file for a particular site for a series of dates.
;
; There are two modes of operation:
; Mode 1) Traditional method driven by the mod_maker.input file
;  which contains the site abbreviation, lat/longitude, and
; the names of the NCEP NCDF-format files.
; Mode 2) New runlog-driven mode in which the time/lat/long
; info comes from the user-selected runlog. This provides
; flexibility for the runlogs to contain multiple sites or a
; moving instrument (e.g. balloon-borne, air-borne, ship-borne).  
;
; INPUTS:
;  runlog.grl 
; OR
;  mod_maker.input   Contains the names of the NetCDF files  
;                     (only first 5 lines are read)
;  ~/ggg/ncdf/Site_xxxxxxxxxx_AT.nc    Air Temperature
;  ~/ggg/ncdf/Site_xxxxxxxxxx_GH.nc    Geopotential Height
; (~/ggg/ncdf/Site_xxxxxxxxxx_TP.nc)   Tropopause pressure (Mode 1 only)
;  ~/ggg/ncdf/Site_xxxxxxxxxx_SH.nc    Specific Humidity
;
; OUTPUTS:
; Mode 1) ~/ggg/models/gnd/ssyyyymmdd.mod  for all dates present in NCEP files
; Mode 2) ~/ggg/models/gnd/NCEP_yyyymmdd_xxN_xxE.mod  for all lat/long/dates in runlog
;
;
; Requires that the user already downloaded to $GGGPATH/ncdf/ the necessary .nc files:
; Air Temperature, Geopotential Height, (Tropopause Pressure), Specific Humidity
; from the NOAA website. For example, to download the entire global dataset:
; wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/pressure/air.2013.nc
; wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/pressure/shum.2013.nc
; wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/pressure/hgt.2013.nc
; wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/tropopause/pres.tropp.2013.nc

; To download data for a smaller geographic region interact with the website:
; http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html
;
;  IMPORTANT: Enter longitude in deg E (eg, Park Falls: 240-280 E)
;             Use 6-hourly NCEP data (4 times daily)
;             Select all available pressure levels (1000 - 10 mb)
;             Select the time period
;             Select "Create subset without making a plot"
; 
;  You must then edit mod_maker.input and enter at the top of the file
;  the lat/long of your site and the names of the four input files.
;  
;  You are now ready to run mod_maker9
;----------------------------------------------------------------
; Version History Log (extracted from ggg.history)
;
; mod_maker_9.5  08-Feb-2013  GCT
; Changed the equation for converting the NCEP Specific Humidity (MMR)
; into a mole fraction (VMR) from
;   rmm=28.964/18.02             ; Ratio of Moleculas Masses (Dry_Air/H2O)
;   ssh=ssh*rmm                  ; Convert Mass Mixing Ratio to VMR
; to
;   rmm=28.964/18.02             ; Ratio of Moleculas Masses (Dry_Air/H2O)
;   ssh=ssh*rmm/(1-ssh*(1-rmm))  ; Convert Mass Mixing Ratio to VMR
;
; This makes little difference. Differences of ~1% under warm humid conditions.
; The old equation produced VMR = 1.607 when MMR = 1.0 which is nonsense.
; The new equation produces VMR = 1.0 when MMR = 1.0
;
; Fixed bug that caused large time extrapolation of NCEP data when day=0
; was entered into the runlog (because the CALDAT-derived year became
; different from the runlog year (IYR)
;
;---------------------------------------------------------------
; mod_maker_9.4  03-Aug-2012  GCT
;  Added code to calculate tropopause altitude internally from NCEP P/T profiles. 
;  instead of using latitude/season to generate a climatological tropopause.
;
;--------------------------------------------------------------
; mod_maker_9.3  29-Jul-2012  GCT
;
; I realized that the longitudes were confusing. There are three
; different flavors:
;  - NCEP requires longitudes that range 0-360 deg
; - The time equation requires longitudes that range -180 to +180
; - The runlog longitudes generally conform to the latter definition,
;   but there's no guarantee.
; The new version of the code therefore has three different variables
; representing the three different flavors of the site longitude:
;
; Resulting from this confusion, there was a mistake in the code that
; calculates the local date from the runlog info (year,doy,zpdtim,oblon)
; To fix this changed
;   julian=julday(1,idoy,iyr,zpdtim)-Site_Lon_360/360
; to
;   julian=julday(1,idoy,iyr,zpdtim)+Site_Lon_180/360
;
; For a site over the US (oblon=-120) this makes no difference
; because  -Site_Lon_360/360 = -0.66
; and      +Site_Lon_180/360 = -0.33
; so the same julian date is chosen.
; For a site at 120E, however, the former gives -0.33 whereas the
; latter gives +0.33, resulting in different dates.
;
;-----------------------------------------------------------
; mod_maker9.pro         GCT  2012-06-01
; Delayed the application of gains and offsets to the raw I*2 NCDF data
; until after tri-linear interpolation. This halves the memory footprint.
;
; The raw NCDF data are 16-bit integers (I*2). In order
; to maintain adequate precision and dynamic range, gains
; and offsets have been derived and removed by NCEP and
; must be re-applied by the user. Unfortunately, application
; of the gains and offsets doubles the memory footprint
; because the data are floated from I*2 to R*4.
; Previously, this was done immediately following reading
; the NCDF file, which means that the gains and offsets
; can then be discarded -- they are not needed again.
;
; In the new version, the gains and offsets are applied *after*
; the raw data have been interpolated. This has two advantages:
; 1) The large arrays containing the raw data stay as I*2
; 2) The gains and offsets are only applied to the data to be used,
;    not the entire dataset
;
; The disadvantage is that the vectors containing the gains and offsets
; for each parameter must be passed from read_global_data.pro, through
; the main program and down into the trilinear_interp.pro subroutines.
; Thanks to Debra and Chris for this suggestion.
;
;-------------------------------------------------------------------
; mod_maker9.pro         GCT  2012-05-30
; Added additional argument (NT) to call to trilinear.  This allows
; to check for array bound-violations, such as can occur on Dec 31
; when observations occur after 18:00 UT. Added two if-statements
; to trilinear to prevent array bound-violations by extrapolating
; in time rather than trying to interpolate.
;
;----------------------------------------------------------------
; mod_maker9.pro         GCT  2012-05-17
;  Moved the code that checks the latitude and longitude ranges of
; the NCEP data into the read_global_data_file subroutine. This way
; the code exists only in one place instead of 4 places. This required
; Site_Lat and Site_Lon to be added to the input arguments.
;
;--------------------------------------------------------------------
; mod_maker9.pro         GCT  2012-05-07
; Converted mod_maker to make it driven either:
; 1) from a user-specified runlog, or
; 2) the ncdf files listed in mod_maker.input.
; So the user can either enter the name of a runlog
; or hit CR and mod_maker will operate in the old way.
;
; The main advantage of the new option (1) is that the
; model files are created based on the lat, long, and dates
; in the runlog. This provides much greater flexibility, for
; example, for instruments that move around (e.g., MkIV) you
; can create all the nescssary models from a single input file
; (the runlog). Ditto for a runlog containing spectra from
; different TCCON sites (e.g., an aircraft overpass runlog).
; This also makes it simpler to generate models on a per
; spectrum basis, if we decide to go this way in the future.
;
; Another advantage is that you can run the new mod_maker
; over multiple years of data, and it will automatically
; switch from one year's NCEP data file to the next year's.
;
; A less important advantage is that model files are created
; only for the days for which spectra exist in the runlog,
; whereas previously they were created for every day of the year.
;
; I've also converted the code that reads the NCDF input
; data into a subroutine. So instead of having four blocks
; of nearly identical code in the main program, there are
; now just four subroutine calls.
; Also the "Dec 31" problem is fixed/kludged in the new version,
; in both operating modes. This problem occurs for observations
; (local solar noon) between 18 UT on Dec 31 and 00 UT on Jan 1.
; This happens in the Western USA when measurements made at
; local solar noon on Dec 31 are at ~20UT.  Theoretically, there
; could also be a Jan 1 problem on an island in the Western Pacific,
; 12+ hours East of the Greenwich meridian, and this possibility is
; covered too.
;
; Anyway, since these NCEP data typically exist in different annual
; files, it is not easy to interpolate between the last (18UT) entry
; on Dec 31 2004 and the first (00UT) entry on Jan 1 2005. This would
; require holding data from both years in memory simultaneously.
; The not-entirely-satisfactory solution to this problem is to
; extrapolate in time by inserting the two following if-statements
; into the code (trilinear_interp) that does the time interpolation
;   index_tt=fix(tt)
;   if index_tt+1 gt nt-1 then index_tt=nt-2  ; Prevent Dec 31 problem
;   if index_tt lt 0 then index_tt=0          ; Prevent Jan 1 problem
;   fr_tt=tt-index_tt
;
; I also decided to stop using the NCDF tropopause pressure
; (TP) file as input. For several years, the GGG code has
; calculated its own tropopause altitudes, which seems just
; as good (if not better) than the values provided by NCEP.
; Omitting the TP files helps simplify the subroutines that
; reads the input data and do the interpolation because the
; site TP is a scalar, whereas the other quantities (AT, GH,
; SH) are vectors (vertical profiles). The only down-side
; of this is that the write_mod subroutine uses the NCEP
; tropopause pressure to compute the H2O profile above the
; highest NCEP level (300 mbar). This is now done depending
; the latitude:
;   efflat=oblat+15*cos(2*3.14159265*(idoy-50)/365.25)
;   xx=(efflat/32.5)^2
;   Site_Noon_TP = 100*(95.0+175.0*xx/(1+xx)) ; Pa
; This changes the water profile for p<300 mbar especially
; when the tropopause pressure is far from climatology.
; But the impact on ground-based H2O retrievals is minimal.
;
; The down side of this more flexible approach is that the
; memory footprint increases because you can no longer
; immediately interpolate the global data to the lat/long
; of the observation site, because there might be multiple
; sites in the runlog.  This means that you can't re-use
; the same global data array for the three data input files
; (AT, GH, SH).  You have to store *all* the global data
; in memory all the time. This more than doubles the total
; memory usage.
;
; If you are running mod_maker on a 32-bit computer and are
; using annual, global, NCDF files as input, this will be a
; problem because IDL has a 2GB memory limit on 32-bit machines.
; [Since the annual, global, NCDF files total only 1.3GB, I'm
; not sure why the memory limit of 2GB is exceeded.  Possibly
; some inefficiency/duplication]. Anyway, on a 64-bit machine,
; there is no memory problem with the new version.
;  521829976  Aug  5  2010  air.2006.nc
;   30708976  Aug  5  2010  pres.tropp.2006.nc
;  245574588  Aug  5  2010  shum.2006.nc
;  521829984  Aug  5  2010  hgt.2006.nc
;
; There is no speed penalty with the new version. It still
; has to do the same interpolation. But whereas previously
; it did it in two steps (immediate bi-linear interpolation
; in lat/long followed by linear interpolation in time),
; it now does it all in a single tri-linear interpolation
; (lat/long/time)
;
; Although the contents of the .mod files will be virtually
; identical (tropopause altitude is missing), the file names
; will no longer have a two-letter site prefix identifying the
; site. Instead, the lat/longitude is built into the name of
; the model file, e.g., instead of pa20041222.mod we now have
; NCEP_20041222_46N_090W.mod
; Separating the model from a specific site will give us more
; flexibility in the future and less chance of an error (e.g.,
; putting an incorrect lat/long in the mod_maker.input file).
; The mod_maker.input file will no longer even be necessary
; since the runlog will be the driving input file.
;
;------------------------------------------------------------
;  mod_maker8.pro    13-Sep-2011    GCT
; Changed mod_maker8.pro so that it does not convert the NCEP geopotential
; altitudes into geometric altitudes. This was wasting memory.  Instead,
; the conversion is done in read_model_fc.f.  In fact, this only affects
; the tropopause altitude whose computation is the only use for the
; altitudes in the .mod file.
;
; Also, I resolved a small earth radius discrepancy in mod_maker8.pro.
; In the write_mod.pro subroutine, 6378.00 was used, whereas in the main
; program 6378.137 km was used.
;
; Finally, changed the parameterization of the stratospheric H2O in
; mod_maker8.pro. The new equation depends on the tropopause pressure
; (the old one didn't).
;
;--------------------------------------------------------------
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
;-------------------------------------------------------------
; mod_maker6  16-Nov-2009  GCT
;  The conversion from geopotential to geometric altitudes
;  was being done incorrectly. Also, I was using inconsistent
;  values for the earth radius between mod_maker6 and readmodFC.f
;
;  So I changed the statement:
;     Earth_Radius=6356.772*1000.0
;     Global_Data = Global_Data /(1+Global_Data/Earth_Radius)
;  to
;     Earth_Radius=6378.137*1000.0
;     Global_Data = Global_Data /(1-Global_Data/Earth_Radius)
;
;  Thanks to Markus Rettinger for pointing out these mistakes.
;
;--------------------------------------------------------------
; mod_maker6.1  03-Oct-2009  GCT
;  Re-organized the if-statements that check that the .nc files
;  downloaded from the NCEP/NCAR website contain sufficient levels.
;  The error messages are now more specific in stating which NCEP
;  file contains insufficient levels.
;
;  Also, I was noticing from fits to ground-based MkIV spectra
;  acquired during MOHAVE2009 that there was often too much H2O
;  in the 10-14 km altitude range. So I modified the code that
;  extrapolates the H2O vmr profiles above the highest NCEP
;  level (300 mbar):
;  1) I made H2O decrease faster above 300 mbar level by using a
;      P(k)/P(k-1)**(5.5-P[k]/100) exponent instead of
;      P(k)/P(k-1)**2.5
;     This has the effect of making the extrapolated H2O decrease
;     faster and faster with increasing altitude
;  2) If the H2O partial pressure (vmr * pressure) exceeds the
;     saturation vapor pressure over ice, it is set to the SVP
;     Note that any super-saturation in the low-altitude levels
;     (300-1000 mbar) covered by the NCEP H2O, is left unmodified.
;  3) If the current pressure level is in the stratosphere
;     (i.e. P < Tropopause Pressure) then use the stratospheric vmr
;
;  Also changed the functional of the stratospheric H2O from
;       strat_h2o=9.7e-06-2.5E-06*alog10(10.0+Lev_AT[k])
;     to
;       strat_h2o=7.5E-06*exp(-0.25*zz^2)+0.5E-09*exp(+4.0*zz)
;
;--------------------------------------------------------------
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
;  
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
;11) Added a subroutine (read_global_data_file) to read the NCDF
;  files and to test of the lat/longitude limits of the data.
;  This way the necessary code exists in just one place instead
;  of four(AT, HG, TP, SH).
;
;12) Replaced the linear_interp01 and bilinear_interp01 subroutines
;  by trilinear_interp01 subroputines. This does the Lat/Long/Time
;  interpolation in one step, whereas previously the Lat/Long
;  interpolation was done immediately, reducing the data volume
;  that must be held in memory, and the time interpolation was
;  done later. The new approach is more flexible and simpler.
;  This disadvantage is that all data (AT, GH, TP, SH) must be
;  held in memory together, beacuse it doesn't assume that the
;  site Lat/Long is going to be fixed. The larger memory footprint
;  means that mod_maker9 won't run on a 32-bit machines.
;
;-----------------------------------------------------------
; Note: IDL is an interpretive language, so subroutines must
; be defined *before* they are called by the main program.
; So the subroutines are found at the beginning of the file
; and then the main program is found near the end.
;
;=================================================================
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
;
;------------------------------------------------------------
; SUBROUTINE WRITE_MOD:
pro write_mod, mod_path, Site_Lat, Lev_AT, sat, sgh, stp, ssh
; Creates a GGG-format .mod file when called by mod_maker.pro.
; INPUTS:
;     Site_Lat    ; The latitude of the site
;     Lev_AT      ; The pressure levels on which the data are tabulated
;     sat         ; Site Noon Atmospheric Temperature profile (vector)
;     sgh         ; Site Noon Geometric Height profile in m (vector)
;     ssh         ; Site Noon Specific Humidity profile MMR (vector)
;     stp         ; Site Noon Tropopause Pressure (scalar)
;
;*************************************************************
; Define US Standard Atmosphere (USSA) for use above 10 mbar
  p_ussa=[10.0,  5.0,   2.0,   1.0,   0.1,   0.01,  0.001, 0.0001]
  t_ussa=[227.7, 239.2, 257.9, 270.6, 231.6, 198.0, 189.8, 235.0]
  z_ussa=[31.1,  36.8,  42.4,  47.8,  64.9,  79.3,  92.0,  106.3]

  rmm=28.964/18.02             ; Ratio of Moleculas Masses (Dry_Air/H2O)
  ssh=ssh*rmm/(1-ssh*(1-rmm))  ; Convert Mass Mixing Ratio to VMR
  sgh=sgh/1000                 ; Convert meters to km
    print,Mod_Path,stp/100
; Export the head of .mod
    openw,lunw,Mod_Path,/get_lun
    printf,lunw,'4  5'
    printf,lunw,6378.137,6.000E-05,Site_Lat,9.81, $
      sgh[0],1013.25, stp/10^2, $
      format='(f8.3,1x,e11.4,1x,f7.3,1(1x,f5.3),3(1x,f8.3))'
    printf,lunw, ' mbar        Kelvin         km      g/mole      vmr'
    printf,lunw, 'Pressure  Temperature     Height     MMW        H2O'
  
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
;
;---------------------------------------------------------
; SUBROUTINE TRILINEAR_INTERP0
; Evaluates  fout = fin(xx,yy,tt)
; Result is a scalar
pro trilinear_interp0, fin, fscale_factor, fadd_offset, xx, yy, tt, nt, fout, nxx
  index_xx=fix(xx)
  index_xx_1 = index_xx+1
  fr_xx=xx-index_xx
  index_yy=fix(yy)
  fr_yy=yy-index_yy
  index_tt=fix(tt)
  if index_tt lt 0 then index_tt=0          ; Prevent Jan 1 problem
  if index_tt+1 gt nt-1 then index_tt=nt-2  ; Prevent Dec 31 problem
  fr_tt=tt-index_tt  ;  Should be between 0 and 1 when interpolating in time
  if(fr_tt lt -1 or fr_tt gt 2) then begin
       print, 'Excessive time extrapolation:',fr_tt,' time-steps   =',fr_tt/4,'days'
       print, ' tt= ',tt,'  index_tt=',index_tt,'  fr_tt=',fr_tt
       print, 'An NCEP file doesnt cover the full range of dates'
      stop
  endif
  if(fr_tt lt -0 or fr_tt gt 1) then print, ' Warning: Time extrapolation of ',fr_tt,' time-steps'
  fout= $
   ((fin[index_xx,index_yy,index_tt]*(1.0-fr_xx)+ $
     fin[index_xx_1 mod nxx,index_yy,index_tt]*fr_xx)*(1.0-fr_yy)+ $
    (fin[index_xx,index_yy+1,index_tt]*(1.0-fr_xx)+ $
     fin[index_xx_1 mod nxx,index_yy+1,index_tt]*fr_xx)*fr_yy)*(1.0-fr_tt)+ $
   ((fin[index_xx,index_yy,index_tt+1]*(1.0-fr_xx)+ $
     fin[index_xx_1 mod nxx,index_yy,index_tt+1]*fr_xx)*(1.0-fr_yy)+ $
    (fin[index_xx,index_yy+1,index_tt+1]*(1.0-fr_xx)+ $
     fin[index_xx_1 mod nxx,index_yy+1,index_tt+1]*fr_xx)*fr_yy)*fr_tt 
   fout=reform(fout,/overwrite)*fscale_factor + fadd_offset
end
;
;---------------------------------------------------------------
; SUBROUTINE TRILINEAR_INTERP1
; Evaluates  fout = fin(xx,yy,*,tt) 
; Result is a 1-vector
pro trilinear_interp1, fin, fscale_factor, fadd_offset, xx, yy, tt, nt, fout, nxx
  index_xx=fix(xx)
  index_xx_1 = index_xx+1
  fr_xx=xx-index_xx
  index_yy=fix(yy)
  fr_yy=yy-index_yy
  index_tt=fix(tt)
  if index_tt lt 0 then index_tt=0          ; Prevent Jan 1 problem
  if index_tt+1 gt nt-1 then index_tt=nt-2  ; Prevent Dec 31 problem
  fr_tt=tt-index_tt  ;  Should be between 0 and 1 when interpolating in time
  if(fr_tt lt -1 or fr_tt gt 2) then begin
       print, 'Excessive time extrapolation:',fr_tt,' time-steps   =',fr_tt/4,'days'
       print, ' tt= ',tt,'  index_tt=',index_tt,'  fr_tt=',fr_tt
       print, 'An NCEP file doesnt cover the full range of dates'
      stop
  endif
  if(fr_tt lt -0 or fr_tt gt 1) then print, ' Warning: Time extrapolation of ',fr_tt,' time-steps'
  fout= $
   ((fin[index_xx,index_yy,*,index_tt]*(1.0-fr_xx)+ $
     fin[index_xx_1 mod nxx,index_yy,*,index_tt]*fr_xx)*(1.0-fr_yy)+ $
    (fin[index_xx,index_yy+1,*,index_tt]*(1.0-fr_xx)+ $
     fin[index_xx_1 mod nxx,index_yy+1,*,index_tt]*fr_xx)*fr_yy)*(1.0-fr_tt)+ $ 
   ((fin[index_xx,index_yy,*,index_tt+1]*(1.0-fr_xx)+ $
     fin[index_xx_1 mod nxx,index_yy,*,index_tt+1]*fr_xx)*(1.0-fr_yy)+ $
    (fin[index_xx,index_yy+1,*,index_tt+1]*(1.0-fr_xx)+ $
     fin[index_xx_1 mod nxx,index_yy+1,*,index_tt+1]*fr_xx)*fr_yy)*fr_tt 
   fout=reform(fout,/overwrite)*fscale_factor + fadd_offset
end
;
;-------------------------------------------------------------------
pro read_global_data_file, ncdf_file_XX, qq, Site_Lat, Site_Lon_360, Lev_XX, Lat_XX, Lon_XX, Tim_XX, Global_Data_XX, Global_Data_scale_factor_XX, Global_Data_add_offset_XX
  ncid = NCDF_OPEN(ncdf_file_XX)          ; Open The NetCDF file
  if (qq ne 'pres') then begin
  varid=ncdf_varid(ncid,'level')
  NCDF_VARGET, ncid,  varid, Lev_XX       ; Read in variable 'level'
  endif
  varid=ncdf_varid(ncid,'lat')
  NCDF_VARGET, ncid,  varid, Lat_XX       ; Read in variable 'lat'
  varid=ncdf_varid(ncid,'lon')
  NCDF_VARGET, ncid,  varid, Lon_XX       ; Read in variable 'lon'
  varid=ncdf_varid(ncid,'time')
  NCDF_VARGET, ncid,  varid, Tim_XX       ; Read in variable 'time'
  varid=ncdf_varid(ncid,qq)
  NCDF_VARGET, ncid, varid, Global_Data_XX    ; Read in variable 'Global_Data_XX'
  NCDF_ATTGET, ncid, varid, 'add_offset', Global_Data_add_offset_XX
  NCDF_ATTGET, ncid, varid, 'scale_factor', Global_Data_scale_factor_XX
  NCDF_CLOSE, ncid     ; Close the NetCDF file
;  Global_Data_XX = Global_Data_XX*Float(Global_Data_scale_factor_XX) + Float(Global_Data_add_offset_XX)
; Check that the range of NCEP Lat/Longitudes span the observation locations

 dlon=Lon_XX[N_Elements(Lon_XX)-1]-Lon_XX[0]
 if (Site_Lat-Lat_XX[N_Elements(Lat_XX)-1])*(Site_Lat-Lat_XX[0]) gt 0 or $
    (Site_Lon_360-Lon_XX[N_Elements(Lon_XX)-1])*(Site_Lon_360-Lon_XX[0])*dlon gt 0 then begin
       print,'The data file '+ncdf_XX_file+' does not cover this Site'
       print,'Latitudes: ',Lat_XX
       print,'Longitudes:',Lon_XX
       stop
 endif
print,'read_global_data ',qq,' Nlev, Nlat, Nlong, Ntime =',N_Elements(Lev_XX), N_Elements(Lat_XX), N_Elements(Lon_XX), N_Elements(Tim_XX)
end

;============================================================================

; MAIN PROGRAM:
pro mod_maker9
; Read site abbreviation, lat/long & names of NCDF files from "mod_maker.input"
; Note that it only reads the first 6 lines of this file.
; So you can keep other stuff further down for re-use later.
;
 print,'mod_maker_9.6   2013-05-24   GCT'

 home_path=string(getenv('GGGPATH'))
 occul=' '
 names=' '
 spectrum=' '
 print,format='($,"Enter runlog (CR to use mod_maker.input)")'
 read,occul
 ix=strlen(occul)
 if ix gt 0 then begin     ; New operating mode using user-specified runlog

 print,'Models will be created for all dates in runlog at the appropriate long/latitudes'
 case strmid(occul,ix-3,1) of
   'g': openr,unit,string(home_path+'/runlogs/gnd/'+occul),/get_lun
   'b': openr,unit,string(home_path+'/runlogs/bal/'+occul),/get_lun
   'o': openr,unit,string(home_path+'/runlogs/orb/'+occul),/get_lun
   'l': openr,unit,string(home_path+'/runlogs/lab/'+occul),/get_lun
 endcase
readf,unit,nlhead,ncol
for i=2,nlhead do begin
  readf,unit,names
endfor
;print,strpos(names,'Year')
ll_year=strpos(names,'Year')
rlformat=string('(1x,a57,1x,2i4,f8.4,f8.3,f9.4,2f8.3,f7.4,f8.3,f7.3,f7.2,3f6.4,2i9,f15.11,i9,i3,39x,f8.2)')
yearwas=0
Mod_Was=' '
nspe=0
nmod=0
nyear=0
while(eof(unit) eq 0) do begin
  readf,unit,format=rlformat,$
  spectrum,iyr,idoy,zpdtim,oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bpw,pout
  nspe=nspe+1
;
; Site_Lat should already be in the range -180 to +180, but just in case not...
  if(oblon gt 180.0 ) then Site_Lon_180=oblon-360.0 else Site_Lon_180=oblon
;
; The NCEP files use longitudes ranging from 0 to 360, so define yet another longitude...
  if(oblon lt 0.0) then Site_Lon_360=oblon+360.0 else Site_Lon_360=oblon
;
  julian=julday(1,idoy,iyr,zpdtim)+Site_Lon_180/360   ; GCT 2012-07-28
  caldat,julian,month,day,year
;  print,'caldat:',iyr,idoy,julian,year,month,day
  if oblat lt 0.0  then ns='S'  else ns='N' 
  if Site_Lon_180 lt 0.0  then ew='W'  else ew='E'
  Mod_Name=string('NCEP_',year,month,day,'_',round(abs(oblat)),ns,'_',round(abs(Site_Lon_180)),ew,'.mod',FORMAT='(A5,I4,I2.2,I2.2,A1,i2.2,A1,A1,I3.3,A1,A4)')
  if year ne yearwas then begin  ;  Read a year of NCDF-format NCEP data
     print, ' Reading NCEP data for year ',year
     nyear=nyear+1
     ncdf_AT_file=string(home_path+'/ncdf/air.',year,'.nc',format='(a,i4.4,a3)')
     ncdf_GH_file=string(home_path+'/ncdf/hgt.',year,'.nc',format='(a,i4.4,a3)')
;     ncdf_TP_file=string(home_path+'/ncdf/pres.tropp.',year,'.nc',format='(a,i4.4,a3)')
     ncdf_SH_file=string(home_path+'/ncdf/shum.',year,'.nc',format='(a,i4.4,a3)')
     read_global_data_file,ncdf_AT_file, 'air',oblat,Site_Lon_360,Lev_AT,Lat_AT,LON_AT,TIM_AT,Global_Data_AT, Global_Data_scale_factor_AT, Global_Data_add_offset_AT
     read_global_data_file,ncdf_GH_file, 'hgt',oblat,Site_Lon_360,Lev_GH,Lat_GH,LON_GH,TIM_GH,Global_Data_GH, Global_Data_scale_factor_GH, Global_Data_add_offset_GH
;     read_global_data_file,ncdf_TP_file,'pres',oblat,Site_Lon_360,Lev_TP,Lat_TP,LON_TP,TIM_TP,Global_Data_TP, Global_Data_scale_factor_TP, Global_Data_add_offset_TP
     read_global_data_file,ncdf_SH_file,'shum',oblat,Site_Lon_360,Lev_SH,Lat_SH,LON_SH,TIM_SH,Global_Data_SH, Global_Data_scale_factor_SH, Global_Data_add_offset_SH
;
; Check that the full number of levels have been downloaded
; Geopotential & Air Temperature files
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
  endif
  yearwas=year

  if Mod_Name ne Mod_Was then begin  ; Interpolate in long/latitude/time and Write Model
    nmod=nmod+1
    Site_Noon_Tim=(julday(month,day,year,12)-julday(1,1,1,0))*24.0 - Site_Lon_180/15
    trilinear_interp1,Global_Data_AT, Global_Data_scale_factor_AT, Global_Data_add_offset_AT,(Site_Lon_360-Lon_AT[0])/(Lon_AT[1]-Lon_AT[0]),(oblat-Lat_AT[0])/(Lat_AT[1]-Lat_AT[0]),(Site_Noon_Tim-Tim_AT[0])/(Tim_AT[1]-Tim_AT[0]), N_Elements(Tim_AT), Site_Noon_AT, N_Elements(lon_AT)
    trilinear_interp1,Global_Data_GH, Global_Data_scale_factor_GH, Global_Data_add_offset_GH,(Site_Lon_360-Lon_GH[0])/(Lon_GH[1]-Lon_GH[0]),(oblat-Lat_GH[0])/(Lat_GH[1]-Lat_GH[0]),(Site_Noon_Tim-Tim_GH[0])/(Tim_GH[1]-Tim_GH[0]), N_Elements(Tim_GH), Site_Noon_GH, N_Elements(lon_GH)
;    trilinear_interp0,Global_Data_TP, Global_Data_scale_factor_TP, Global_Data_add_offset_TP,(Site_Lon_360-Lon_TP[0])/(Lon_TP[1]-Lon_TP[0]),(oblat-Lat_TP[0])/(Lat_TP[1]-Lat_TP[0]),(Site_Noon_Tim-Tim_TP[0])/(Tim_TP[1]-Tim_TP[0]), N_Elements(Tim_TP), Site_Noon_TP, N_Elements(lon_TP)
    trilinear_interp1,Global_Data_SH, Global_Data_scale_factor_SH, Global_Data_add_offset_SH,(Site_Lon_360-Lon_SH[0])/(Lon_SH[1]-Lon_SH[0]),(oblat-Lat_SH[0])/(Lat_SH[1]-Lat_SH[0]),(Site_Noon_Tim-Tim_SH[0])/(Tim_SH[1]-Tim_SH[0]), N_Elements(Tim_SH), Site_Noon_SH, N_Elements(lon_SH)

; Determine tropopause altitude from temperature profile
; (instead of using NCEP tropopause pressure product)
; Search upwards for first instance of lapse_rate > -2 K/km above 5000 m.
; This prevents T-inversions in the PBL being mistaken for the
; tropopause, but will occasionally lose a true tropopause below 5000m.
    lr_thresh=-0.002  ;  Lapse Rate Threshold (-2 K/km) = -0.002 K/m)
    k=2
    lapse_was= (Site_Noon_AT(k-1)-Site_Noon_AT(k-2))/(Site_Noon_GH(k-1)-Site_Noon_GH(k-2))  
    lapse_rate=(Site_Noon_AT(k)-Site_Noon_AT(k-1))/(Site_Noon_GH(k)-Site_Noon_GH(k-1))
    zlrwas=0.5*(Site_Noon_GH(k-1)+Site_Noon_GH(k-2)) 
    zlr=0.5*(Site_Noon_GH(k)+Site_Noon_GH(k-1))
    while ( lapse_rate lt lr_thresh or Site_Noon_GH(k) lt 5000.0 ) do begin
       k=k+1
       lapse_was=lapse_rate
       zlrwas=zlr
       lapse_rate=(Site_Noon_AT(k)-Site_Noon_AT(k-1))/(Site_Noon_GH(k)-Site_Noon_GH(k-1))
       zlr=0.5*(Site_Noon_GH(k)+Site_Noon_GH(k-1)) 
    endwhile
    ztrop=zlrwas+(zlr-zlrwas)*(lr_thresh-lapse_was)/(lapse_rate-lapse_was)
    sh=(Site_Noon_GH(k)-Site_Noon_GH(k-2))/alog(Lev_AT(k-2)/Lev_AT(k))
    Site_Noon_TP=100*Lev_AT(k-1)*exp(-(ztrop-Site_Noon_GH(k-1))/sh) ; in Pa
;    print,'ztrop,ptrop=',ztrop,Site_Noon_TP

;  Write the NCEPxxxxxx.mod file into the appropriate sub-directory
    case strmid(occul,ix-3,1) of
      'g': write_mod,home_path+'/models/gnd/'+Mod_Name,oblat,Lev_AT,Site_Noon_AT, $
    Site_Noon_GH,Site_Noon_TP,Site_Noon_SH
      'b': write_mod,home_path+'/models/bal/'+Mod_Name,oblat,Lev_AT,Site_Noon_AT, $
    Site_Noon_GH,Site_Noon_TP,Site_Noon_SH
      'o': write_mod,home_path+'/models/orb/'+Mod_Name,oblat,Lev_AT,Site_Noon_AT, $
    Site_Noon_GH,Site_Noon_TP,Site_Noon_SH
      'l': write_mod,home_path+'/models/lab/'+Mod_Name,oblat,Lev_AT,Site_Noon_AT, $
    Site_Noon_GH,Site_Noon_TP,Site_Noon_SH
    endcase
;    write_mod,home_path+'/models/gnd/'+Mod_Name,oblat,Lev_AT,Site_Noon_AT, $
;    Site_Noon_GH,Site_Noon_TP,Site_Noon_SH
  endif
  Mod_Was=Mod_Name
endwhile
close,unit
free_lun,unit
print,'Number of spectra',nspe
print,'Number of days/models',nmod
print,'Number of years',nyear

endif else begin  ;  Old operating mode using mod_maker.input file

print,'mod_maker.input will be used for Site Lat/Longitudes and dates'
 Site_Abbrev='xx'
 ncdf_AT_file=' '
 ncdf_GH_file=' '
 ncdf_TP_file=' '
 ncdf_SH_file=' '
 openr, unitr,home_path+'/src/idl/mod_maker.input',/get_lun
 readf, unitr,Site_Abbrev
 readf, unitr,Site_Lat, slong
 readf, unitr,ncdf_AT_file
 readf, unitr,ncdf_GH_file
 readf, unitr,ncdf_TP_file
 readf, unitr,ncdf_SH_file
 close, unitr
 free_lun, unitr
; slong should already be in the range -180 to +180, but just in case not...
if(slong gt 180.0 ) then Site_Lon_180=slong-360.0 else Site_Lon_180=slong
;
; The NCEP files use longitudes ranging from 0 to 360, so define yet another longitude...
if(slong lt 0.0) then Site_Lon_360=slong+360.0 else Site_Lon_360=slong

 read_global_data_file,ncdf_AT_file, 'air',Site_Lat,Site_Lon_360,Lev_AT,Lat_AT,LON_AT,TIM_AT,Global_Data_AT,Global_Data_scale_factor_AT, Global_Data_add_offset_AT
 read_global_data_file,ncdf_GH_file, 'hgt',Site_Lat,Site_Lon_360,Lev_GH,Lat_GH,LON_GH,TIM_GH,Global_Data_GH,Global_Data_scale_factor_GH, Global_Data_add_offset_GH
 read_global_data_file,ncdf_TP_file,'pres',Site_Lat,Site_Lon_360,Lev_TP,Lat_TP,LON_TP,TIM_TP,Global_Data_TP,Global_Data_scale_factor_TP, Global_Data_add_offset_TP
 read_global_data_file,ncdf_SH_file,'shum',Site_Lat,Site_Lon_360,Lev_SH,Lat_SH,LON_SH,TIM_SH,Global_Data_SH,Global_Data_scale_factor_SH, Global_Data_add_offset_SH 

; Check that the full number of levels have been downloaded
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

; Define the name of the models
; Naming convention: 14 characters
;   1-2: Site name e.g. 'ca'-Caltech, 'pk'-Park Falls
;   3-10: year+month+day, in the form yyyymmdd
;   11-14: suffix '.mod'
  julian=Tim_AT[0]/24.0+julday(1,1,1,0)
  caldat,julian,month,day,year
  Site_Noon_Tim=(julday(month,day,year,12)-julday(1,1,1,0))*24.0 $
    -Site_Lon_180/15

; Main Loop over model dates
  while ( Site_Noon_Tim le Tim_AT[N_Elements(Tim_AT)-1] $
    and   Site_Noon_Tim le Tim_GH[N_Elements(Tim_GH)-1] $
    and   Site_Noon_Tim le Tim_TP[N_Elements(Tim_TP)-1] $
    and   Site_Noon_Tim le Tim_SH[N_Elements(Tim_SH)-1] ) do begin
  
; Interpolate in longitude (xx), latitude (yy), and time (tt)
    trilinear_interp1,Global_Data_AT,Global_Data_scale_factor_AT, Global_Data_add_offset_AT, (Site_Lon_360-Lon_AT[0])/(Lon_AT[1]-Lon_AT[0]),(Site_Lat-Lat_AT[0])/(Lat_AT[1]-Lat_AT[0]),(Site_Noon_Tim-Tim_AT[0])/(Tim_AT[1]-Tim_AT[0]), N_Elements(Tim_AT), Site_Noon_AT, n_elements(Lon_AT)
    trilinear_interp1,Global_Data_GH,Global_Data_scale_factor_GH, Global_Data_add_offset_GH, (Site_Lon_360-Lon_GH[0])/(Lon_GH[1]-Lon_GH[0]),(Site_Lat-Lat_GH[0])/(Lat_GH[1]-Lat_GH[0]),(Site_Noon_Tim-Tim_GH[0])/(Tim_GH[1]-Tim_GH[0]), N_Elements(Tim_GH), Site_Noon_GH, n_elements(Lon_GH)
    trilinear_interp0,Global_Data_TP,Global_Data_scale_factor_TP, Global_Data_add_offset_TP, (Site_Lon_360-Lon_TP[0])/(Lon_TP[1]-Lon_TP[0]),(Site_Lat-Lat_TP[0])/(Lat_TP[1]-Lat_TP[0]),(Site_Noon_Tim-Tim_TP[0])/(Tim_TP[1]-Tim_TP[0]), N_Elements(Tim_TP), Site_Noon_TP, n_elements(Lon_TP)
    trilinear_interp1,Global_Data_SH,Global_Data_scale_factor_SH, Global_Data_add_offset_SH, (Site_Lon_360-Lon_SH[0])/(Lon_SH[1]-Lon_SH[0]),(Site_Lat-Lat_SH[0])/(Lat_SH[1]-Lat_SH[0]),(Site_Noon_Tim-Tim_SH[0])/(Tim_SH[1]-Tim_SH[0]), N_Elements(Tim_SH), Site_Noon_SH, n_elements(Lon_SH)

;  Write the .mod file 
    caldat,julian,month,day,year
    Mod_Name=string(Site_Abbrev,year,month,day,'.mod',FORMAT='(A2,I4,I2.2,I2.2,A4)')
    write_mod,home_path+'/models/gnd/'+Mod_Name, Site_Lat,Lev_AT,Site_Noon_AT, $
    Site_Noon_GH,Site_Noon_TP,Site_Noon_SH
    julian=julian+1
    Site_Noon_Tim=Site_Noon_Tim+24
  endwhile

endelse
end
