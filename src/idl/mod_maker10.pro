;  ! Number of years of data downloaded mod_maker10.pro plus subroutines
; Overview
;*************************************************************
; Interpolates in time, lat/longitude into the NCEP re-analyses
; to produce a series of model files.  The time/lat/long
; info comes from the user-selected runlog. This provides
; flexibility for runlogs to contain multiple sites or a
; moving instrument (e.g. balloon-borne, air-borne, ship-borne).  
;
; INPUTS:
;  runlog.grl 
;  ~/ggg/ncdf/Site_xxxxxxxxxx_AT.nc    Air Temperature
;  ~/ggg/ncdf/Site_xxxxxxxxxx_GH.nc    Geopotential Height
;  ~/ggg/ncdf/Site_xxxxxxxxxx_SH.nc    Specific Humidity
;
; OUTPUTS:
;  ~/ggg/models/gnd/NCEP_yyyymmdd_xxN_xxE.mod  for all lat/long/dates in runlog
;
;
; Requires that the user already downloaded to $GGGPATH/ncdf/ the necessary .nc files:
; Air Temperature, Geopotential Height, (Tropopause Pressure), Specific Humidity
; from the NOAA website. For example, to download the entire global dataset:
; wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/pressure/air.2013.nc
; wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/pressure/shum.2013.nc
; wget ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/pressure/hgt.2013.nc

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
;  You are now ready to run mod_maker

;  Note: IDL is an interpretive language, so subroutines must
;  be defined *before* they are called by the main program.
;  So the subroutines are found at the beginning of the file
;  and then the main program is found near the end.
;
;----------------------------------------------------------------
; Version History Log
;
; mod_maker_10.8  2019-06-21  GCT
;   Fixed a bug that caused dumplicate entries in the .mod file for
;   the 10 mbar level: one from NCEP and another from the USSA
;----------------------------------------------------------------
; mod_maker_10.7  2018-10-05  GCT
;   Expanded and corrected the USSA data inside the write_mod
;   subroutine.  The altitude of the 5 mbar level was off by 1 km.
;   Added levels at 0.5 and 0.2 mbar.
;----------------------------------------------------------------
; mod_maker_10.6  2017-04-11  GCT
;   Instead of have a fixed value for the mean Molecular weight (MMW)
;   of air, we now used the NCEP H2O to compute a better value using
;   the equation:
;           mmw=28.964*(1-h2o_wmf)+18.02*h2o_wmf
;   at each altitude, where h2o_wmf is the wet mole fraction.
;
;   Also added the mod-maker version number to the header of each .mod file
;
;---------------------------------------------------------------
; mod_maker_10.5  2017-04-04  GCT
;   A third attempt at getting the H2O DMF correct.
;
;   Prior to Feb 2013 the equation was
;      h2o_dmf=sh*rmm 
;   where sh is the NCEP Specific Humidity and
;   rmm=28.964/18.02    ; Ratio of Molecular Masses (Dry_Air/H2O)
;   This neglected to do the dry to wet conversion.
;
;   In Feb 2013 this was changed to:
;      h2o_dmf=sh*rmm/(1-sh*(1-rmm)) ; Convert Mass Mixing Ratio to Mole Fraction
;   which does the dry-to-wet conversion, but took (incorrect)
;   steps to prevent h2o_dmf from ever exceeding 1.0.
;
;   and now the hopefully-correct equation:
;      h2o_dmf=sh*rmm/(1-sh)  ; Convert Wet Mass Mixing Ratio to Dry Mole Fraction
;
;   I was clearly confused in Feb 2013 by the possibility that 
;   the H2O_DMF could exceed 1.0. Now I realize that the H2O 
;   DMF can exceed 1 -- it is the Wet or Real MF that can't.
;
;---------------------------------------------------------------
; mod_maker_10.4  2017-02-18  GCT
;  Modified to read the runlog format from the runlog header,
;  rather than having it hard-wired inside the mod_maker code.
;  If not present in the runlog header, only then uses the
;  hardwired format.
;
;---------------------------------------------------------------
; mod_maker_10.3  2015-08-05  GCT
;  Added if-statement to prevent super-saturation in the 1000-300 mbar
;  range. Previously we relied on NCEP to avoid this problem and only
;  worried about super-saturation at higher altitudes. But there have
;  been reports of occasional very large super-saturations in  NCEP 
;  in the lower troposphere
;
;----------------------------------------------------------------
; mod_maker_10.2  2015-01-16  GCT
;  Removed the code that supported the old way of running mod_maker:
;  using the mod_maker.input file.  This allowed the trilinear_inter0
;  subroutine to be deleted, and various other simplifications.
;
;----------------------------------------------------------------
; mod_maker_10.1  2015-01-15  GCT
;  Removed if-statement at end of read_global_data subroutine
;  that checked that lat/longitude limits of NCDF input file
;  encompassed the site location. This was only testing first
;  spectrum of each year anyway (whenever a new NCDF file was
;  read), and was identifying site longitudes between 357.5 and
;  0.0 as out of range, because it did not recognize the periodic
;  nature of longitue: that long(0)=long(360).
;
;  This if-statement is now replaced by tests on the lat/long
;  indices inside the tri_linear_interp1 subroutines. This has
;  the advantage of being done for every spectrum, not just the
;  first of each year.
;
;  Modified the argument list for the trilinear_interp subroutines.
;  Instead of inputting the lat/ling/time in NCDF indices, they are
;  now input in natural units (deg for lat/long) and the conversion
;  into NCDF indices occurs inside the routines.  This facilitates
;  testing that the requested lat/long/time lies within the bounds
;  of the NCDF file.
;
;----------------------------------------------------------------
; mod_maker_10.0  20-Jun-2014  GCT
; [Yes, this change was performed *before* 9.7]
; Changed the H2O profile in three ways:
; 1) In the 300-1000 mbar range, if the NCEP H2O is such that RH < 30/P_mbar,
;  which corresponds to 3% at 1000 mbar, 6% at 500 mbar, and 10% at 300 mbar,
;  the H2O vmr is replaced with the value that gives RH=30/P_mbar for that T
; 2) For P<300 mbar, the H2O becomes a weighted average of two terms:
;  - trop_h2o: the RH is assumed constant at the 300 mbar value.
;  - strat_h2o: the H2O vmr is a gaussian in log(P) peaking at 7.5 ppm at 1 mbar
; the weights for trop_h2o are (P/300)^3 and for strat_h2o 1-(P/300)^3
; so at the 300 mbar level the strat_h2o has no weight. But by 150 mbar,
; trop_h2o has a weight of only 0.125 and strat_h2o has a weight of 0.875.
; 3) For P<300 mbar, if the H2O RH exceeds 100% then the H2O vmr is set
;  to SVP_H2O(T) so that RH=100%
; Note that in this new formulation the H2O is independent of the
; explicit tropopause altitude, although it does depend on the
; temperature, which implicitly contains tropopause information.
; So the two methods of using mod_maker (mod_maker.input & runlog)
; should now produce identical H2O profiles.
;
;---------------------------------------------------------------
; mod_maker_9.7  24-Oct-2014  DGF (Dietrich Feist)
;  Modified mod_maker9 to handle NCEP netcdf3 (pre Oct 2014) as well as netcdf4 files
;  (post Oct 2014).  Output mod files are identical to previous version of mod_maker9.pro
;  when netcdf3 files are used.  Small rounding differences at the last digit
;  (+/-1 mK for temperature!) may occur when netcdf4 files are used.
;
;---------------------------------------------------------------
; mod_maker_9.6  29-Oct-2013  GCT/DW
; Several minor changes:
; - Got rid of unnecessary fix calls.
; - Modified to handle longitudes between -2.5 and 0.
; - Improved/updated comments at beginning. Deleted some superfluous lines of code.
; - Fixed bug that caused large time extrapolation of NCEP data when day=0 entered
;  into the runlog (because CALDAT-derived year became different from the runlog year (IYR)).
;
;------------------------------------------
; mod_maker_9.5  08-Feb-2013  GCT
; Changed the equation for converting the NCEP Specific Humidity (a MMR)
; into a mole fraction from:
;   rmm=28.964/18.02               ; Ratio of Molecular Masses (Dry_Air/H2O)
;   h2omf=sh*rmm                  ; Convert Mass Mixing Ratio to Mole Fraction
; to
;   rmm=28.964/18.02               ; Ratio of Molecular Masses (Dry_Air/H2O)
;   h2omf=sh*rmm/(1-sh*(1-rmm))  ; Convert Mass Mixing Ratio to Mole Fraction
;
; Differences are only ~1% under warm humid conditions.
; The old equation produced h2omf = 1.607 when MMR = 1.0 which is nonsense.
; The new equation produces h2omf = 1.0 when MMR = 1.0
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
; different sources/users:
; - NCEP requires longitudes that range 0-360 deg
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
;   efflat=Site_Lat+15*cos(2*3.14159265*(idoy-50)/365.25)
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
;
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
;  by trilinear_interp01 subroutines. This does the Lat/Long/Time
;  interpolation in one step, whereas previously the Lat/Long
;  interpolation was done immediately, reducing the data volume
;  that must be held in memory, and the time interpolation was
;  done later. The new approach is more flexible and simpler.
;  The disadvantage is that all data (AT, GH, TP, SH) must be
;  held in memory together, beacuse it doesn't assume that the
;  site Lat/Long is going to be fixed. The larger memory footprint
;  means that mod_maker9 won't run on a 32-bit machines.
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
pro write_mod, mod_path, version, Site_Lat, Lev_AT, sat, sgh, sntp, h2o_dmf
; Creates a GGG-format .mod file when called by mod_maker.pro.
; INPUTS:
;     Site_Lat    ; The latitude of the site
;     Lev_AT      ; The pressure levels on which the data are tabulated
;     sat         ; Site Noon Atmospheric Temperature profile (vector)
;     sgh         ; Site Noon Geopotential Height profile in km (vector)
;     sntp        ; Site Noon Tropopause Pressure (scalar)
;     h2o_dmf     ; Site Noon H2O Dry Mole Fraction (VMR) profile (vector)
;
;  Although the input and outputs to thissubroutine are H2O DMFs,
;  internally WMFs are computed to compare with the SVP to check
;  for humidities that are too high/low. This is because the partial
;  pressure of H2O is the total pressure times the WMF.
;*************************************************************
; Define US Standard Atmosphere (USSA) for use above 10 mbar
; h_ussa is the geopotential altitude, matching those from NCEP
  n_ussa=12  ; number of levels for which USSA is tabulated below.
  p_ussa=[ 10.0,  5.0,  2.0,  1.0,  0.5,  0.2,  0.1, 0.01, .001,.0001,.00001,.000001]
  t_ussa=[227.7,239.3,257.8,270.7,264.3,245.2,231.6,198.0,189.8,235.0, 390.0, 750.0]
  h_ussa=[ 31.0, 35.7, 42.4, 47.8, 53.3, 60.1, 64.9, 79.3, 92.0,106.3, 124.0, 175.0]

; Export the head of .mod
    openw,lunw,Mod_Path,/get_lun
    printf,lunw,'5  6'
    printf,lunw,6378.137,6.000E-05,Site_Lat,9.81, $
      sgh[0],1013.25, sntp/10^2, $
      format='(f8.3,1x,e11.4,1x,f7.3,1(1x,f5.3),3(1x,f8.3))'
    printf,lunw,version
    printf,lunw, ' mbar        Kelvin         km      g/mole      DMF       %'
    printf,lunw, 'Pressure  Temperature    GPHeight    MMW        H2O      RH'
  
; Export the Pressure, Temp and SHum for lower levels (1000 to 300 mbar)
    for k=0,N_elements(h2o_dmf)-1 do begin
      svp_wv_over_ice, sat[k], svp
      h2o_wmf=h2o_dmf[k]/(1+h2o_dmf[k]) ; Wet Mole Faction of H2O
      frh=h2o_wmf*Lev_AT[k]/svp       ; Fractional Relative Humidity
; Relace H2O mole fractions that are too small
      if (frh lt 30./Lev_AT[k] ) then begin
        print, 'Replacing too-small H2O ',mod_path, Lev_AT[k],h2o_wmf,svp*30./Lev_AT[k]/Lev_AT[k],frh,30./Lev_AT[k]
        frh=30./Lev_AT[k]
        h2o_wmf=svp*frh/Lev_AT[k]
        h2o_dmf[k]=h2o_wmf/(1-h2o_wmf)
      endif
; Relace H2O mole fractions that are too large (super-saturated)  GCT 2015-08-05
      if (frh gt 1.0) then begin
        print, 'Replacing too-large H2O ',mod_path,Lev_AT[k],h2o_wmf,svp/Lev_AT[k],frh,1.0
        frh=1.0
        h2o_wmf=svp*frh/Lev_AT[k]
        h2o_dmf[k]=h2o_wmf/(1-h2o_wmf)
      endif
      mmw=28.964*(1-h2o_wmf)+18.02*h2o_wmf
      printf,lunw,Lev_AT[k], sat[k], sgh[k],mmw,h2o_dmf[k],100*frh, $
        format='(e9.3,4x,f7.3,4x,f7.3,4x,f7.4,4x,1e9.3,f6.1)'
    endfor

; Export Pressure and Temp for middle levels (250 to 10 mbar)
; which have no SHum reanalysis.
    ptop=Lev_AT[k-1] ; Top pressure level covered by NCEP H2O
    frh_top=frh  ; remember the FRH at the top (300 mbar) level
;    print,ptop,frh_top 
    for k=N_elements(h2o_dmf),n_elements(Lev_AT)-1 do begin
        zz=alog10(Lev_AT[k])  ; log10[pressure]
        strat_wmf=7.5E-06*exp(-0.16*zz^2)
        svp_wv_over_ice, sat[k], svp

        trop_wmf=frh_top*svp/Lev_AT[k]
        wt=(Lev_AT[k]/ptop)^3
        avg_wmf = trop_wmf*wt + strat_wmf*(1-wt)
        avg_frh = avg_wmf*Lev_AT[k]/svp
        if (avg_frh gt 1.0 ) then begin
          print, 'Replacing super-saturated H2O ',mod_path, Lev_AT[k],avg_wmf,svp*avg_frh/Lev_AT[k],avg_frh,1.0
          avg_frh = 1.0
          avg_wmf= svp*avg_frh/Lev_AT[k]
        endif
        mmw=28.964*(1-avg_wmf)+18.02*avg_wmf
        printf,lunw,Lev_AT[k], sat[k], $
        sgh[k],mmw,avg_wmf/(1-avg_wmf),100*avg_frh, $
        format='(e9.3,4x,f7.3,4x,f7.3,4x,f7.4,4x,1e9.3,f6.1)'
    endfor

; Get the difference between the USSA and given site temperature at 10 mbar,
    Delta_T=sat[16]-t_ussa[0]

; Export the P-T profile above 10mbar, starting at 5 mbar (k=1).
    for k=1,n_ussa-1 do begin
        Delta_T=Delta_T/2
        zz=alog10(p_ussa[k])  ; log10[pressure]
        strat_wmf=7.5E-06*exp(-0.16*zz^2)
        svp_wv_over_ice, sat[k], svp
        mmw=28.964*(1-strat_wmf)+18.02*strat_wmf
        printf,lunw,p_ussa[k],t_ussa[k]+Delta_T,h_ussa[k],mmw,strat_wmf,100*strat_wmf*Lev_AT[k]/svp, $
        format='(e9.3,4x,f7.3,4x,f7.3,4x,f7.4,4x,1e9.3,f6.1)'
    endfor
    free_lun,lunw
end
;
;---------------------------------------------------------------
; SUBROUTINE TRILINEAR_INTERP1
; Evaluates  fout = fin(xx,yy,*,tt) 
; Result is a 1-vector
pro trilinear_interp1, fin, fscale_factor, fadd_offset, Site_Lon_360, Lon_XX, Site_Lat, Lat_XX, Site_Noon_Tim, Tim_XX, fout
;pro trilinear_interp1, fin, fscale_factor, fadd_offset, xx, yy, tt, nt, fout, nxx

dx = Lon_XX[1]-Lon_XX[0]
dy = Lat_XX[1]-Lat_XX[0]
dt = Tim_XX[1]-Tim_XX[0]

xx = (Site_Lon_360-Lon_XX[0])/dx
yy = (Site_Lat-Lat_XX[0])/dy
tt = (Site_Noon_Tim-Tim_XX[0])/dt

nxx =  N_Elements(Lon_XX)
nyy =  N_Elements(Lat_XX)
ntt =  N_Elements(Tim_XX)

  index_xx=fix(xx)
;  if index_xx ge nxx then index_xx eq index_xx-nxx
;  if index_xx lt 0 then index_xx eq index_xx+nxx
  ixpomnxx=(index_xx+1) mod nxx
  fr_xx=xx-index_xx

  index_yy=fix(yy)
  if index_yy gt nyy-2 then index_yy=nyy-2  ;  avoid array-bound violation at SP
  fr_yy=yy-index_yy
;  print,nyy,yy,index_yy,fr_yy

  index_tt=fix(tt)
  if index_tt lt 0 then index_tt=0          ; Prevent Jan 1 problem
  if index_tt+1 gt ntt-1 then index_tt=ntt-2  ; Prevent Dec 31 problem
  fr_tt=tt-index_tt  ;  Should be between 0 and 1 when interpolating in time

  if(fr_tt lt -1 or fr_tt gt 2) then begin
       print, 'Excessive time extrapolation:',fr_tt,' time-steps   =',fr_tt*dt,' days'
       print, ' tt= ',tt,'  index_tt=',index_tt,'  fr_tt=',fr_tt
       print, 'An NCEP file doesnt cover the full range of dates'
      stop
  endif

  if(fr_xx lt 0 or fr_xx gt 1) then begin
       print, 'Excessive longitude extrapolation:',fr_xx,' steps   =',fr_xx*dx,' deg'
       print, ' xx= ',xx,'  index_xx=',index_xx,'  fr_xx=',fr_xx
       print, 'NCEP file doesnt cover the full range of longitudes'
      stop
  endif
  if(fr_yy lt 0 or fr_yy gt 1) then begin
       print, 'Excessive latitude extrapolation:',fr_yy-1,' steps   =',(fr_yy-1)*dy,' deg'
       print, ' yy= ',yy,'  index_yy=',index_yy,'  fr_yy=',fr_yy
       print, 'NCEP file doesnt cover the full range of latitudes'
      stop
  endif
  if(fr_tt lt -0 or fr_tt gt 1) then print, ' Warning: Time extrapolation of ',fr_tt,' time-steps'
  if(fr_xx lt -0 or fr_xx gt 1) then print, ' Warning: Longitude extrapolation of ',fr_xx,' steps'
  if(fr_yy lt -0 or fr_yy gt 1) then print, ' Warning: Latitude extrapolation of ',fr_yy,' steps'
  fout= $
   ((fin[index_xx,index_yy,*,index_tt]*(1.0-fr_xx)+ $
     fin[ixpomnxx,index_yy,*,index_tt]*fr_xx)*(1.0-fr_yy)+ $
    (fin[index_xx,index_yy+1,*,index_tt]*(1.0-fr_xx)+ $
     fin[ixpomnxx,index_yy+1,*,index_tt]*fr_xx)*fr_yy)*(1.0-fr_tt)+ $ 
   ((fin[index_xx,index_yy,*,index_tt+1]*(1.0-fr_xx)+ $
     fin[ixpomnxx,index_yy,*,index_tt+1]*fr_xx)*(1.0-fr_yy)+ $
    (fin[index_xx,index_yy+1,*,index_tt+1]*(1.0-fr_xx)+ $
     fin[ixpomnxx,index_yy+1,*,index_tt+1]*fr_xx)*fr_yy)*fr_tt 
   fout=reform(fout,/overwrite)*fscale_factor + fadd_offset
end
;
;-------------------------------------------------------------------
pro read_global_data_file, ncdf_file_XX, qq, Site_Lat, Site_Lon_360, Lev_XX, Lat_XX, Lon_XX, Tim_XX, Global_Data_XX, Global_Data_scale_factor_XX, Global_Data_add_offset_XX,julday0
  ncid = NCDF_OPEN(ncdf_file_XX)          ; Open The NetCDF file
  if (qq ne 'pres') then begin
  varid=ncdf_varid(ncid,'level')
  NCDF_VARGET, ncid,  varid, Lev_XX       ; Read in variable 'level'
  endif
  varid=ncdf_varid(ncid,'lat')
  NCDF_VARGET, ncid,  varid, Lat_XX       ; Read in variable 'lat'
  varid=ncdf_varid(ncid,'lon')
  NCDF_VARGET, ncid,  varid, Lon_XX       ; Read in variable 'lon'

;
; Modifications by D. Feist to make program work also with NCEP ncdf files after October 2014
;

; Read time
  varid=ncdf_varid(ncid,'time')
  NCDF_VARGET, ncid,  varid, Tim_XX       ; Read in variable 'time'

; Set julday0 => Julian date where time=0
  NCDF_ATTGET, ncid, varid, 'units', time_units  ; string containing definition of time units
  ; Quick and dirty time unit parsing - this could be more intelligent (D.F.)
  IF string(time_units) EQ 'hours since 1-1-1 00:00:0.0' THEN julday0=julday(1,1,1,0)
  IF string(time_units) EQ 'hours since 1800-01-01 00:00:0.0' THEN julday0=julday(1,1,1800,0)

; Read global data
  varid=ncdf_varid(ncid,qq)
  NCDF_VARGET, ncid, varid, Global_Data_XX    ; Read in variable 'Global_Data_XX'

; Initialize add_offset and scale_factor
  Global_Data_add_offset_XX=0
  Global_Data_scale_factor_XX=1.0

; Check variable for offset and scaling factor attributes
  varinq=NCDF_VARINQ(ncid, varid)
  FOR attid=0,varinq.NATTS-1 DO BEGIN
    ; Set Global_Data_add_offset_XX
    IF NCDF_ATTNAME(ncid, varid, attid) EQ 'add_offset' THEN NCDF_ATTGET, ncid, varid, 'add_offset', Global_Data_add_offset_XX
    ; Set Global_Data_scale_factor_XX
    IF NCDF_ATTNAME(ncid, varid, attid) EQ 'scale_factor' THEN NCDF_ATTGET, ncid, varid, 'scale_factor', Global_Data_scale_factor_XX
  ENDFOR
;
; End of modification
;
  NCDF_CLOSE, ncid     ; Close the NetCDF file
end

;============================================================================

; MAIN PROGRAM:
pro mod_maker10, other_args
 version='mod_maker_10.9   2019-07-10   GCT'
 print,version
input_args = command_line_args(count=nargs)
; Edit: Added option to input runlogs in from command line
; e.g. $idl -e "mod_maker9b" -args pa_ggg_benchmark.grl
 

 home_path=string(getenv('GGGPATH'))
 occul=' '
 names=' '
 col1=' '
 spectrum=' '
 apf=' '
 Mod_Was=' '
 rlformat=string('(a1,a57,1x,2i4,f8.4,f8.3,f9.3,2f8.3,1x,f6.4,f8.3,f7.3,f7.2,3(1x,f5.4),2i9,1x,f14.11,i9,i3,1x,f5.3,i5,1x,a2,2(f6.1,f8.2,f5.1),f7.1,f7.4,f6.1,f6.0,f10.3,f7.0,f7.3)')

 if (nargs ne 1) then begin  ; if zero input
  print,format='($,"Enter runlog")'
  read,occul
 endif
 if (nargs eq 1) then occul=input_args

 ix=strlen(occul)
 case strmid(occul,ix-3,1) of
   'g': ext=string('gnd/')
   'b': ext=string('bal/')
   'o': ext=string('orb/')
   'l': ext=string('lab/')
 endcase
keyword=string('       ')
openr,unit,string(home_path+'/runlogs/'+ext+occul),/get_lun
readf,unit,nlhead,ncol
for i=2,nlhead do begin
  readf,unit,names
  if strpos(names,'format=') ge 0 then rlformat=strmid(names,7)
endfor
;print,strpos(names,'Year')
ll_year=strpos(names,'Year')
yearwas=0
nspe=0UL
nmod=0
nyear=0  ; Number of years of data downloaded
 print,'Models will be created for all dates in runlog at the appropriate long/latitudes'
; Loop over spectra
while(eof(unit) eq 0) do begin
   readf,unit,format=rlformat,$
   col1,spectrum,iyr,idoy,zpdtim,Site_Lat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bpw,zoff,snr,apf,tins,pins,hims,tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl
   print, spectrum
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
;   print,'caldat:',iyr,idoy,julian,year,month,day
   if Site_Lat ge 0.0  then ns='N'  else ns='S' 
   if Site_Lon_180 ge 0.0  then ew='E'  else ew='W'
   Mod_Name=string('NCEP_',year,month,day,'_',round(abs(Site_Lat)),ns,'_',round(abs(Site_Lon_180)),ew,'.mod',FORMAT='(A5,I4,I2.2,I2.2,A1,i2.2,A1,A1,I3.3,A1,A4)')
   if year ne yearwas then begin  ;  Read a year of NCDF-format NCEP data
      nyear=nyear+1
      ncdf_AT_file=string(home_path+'/ncdf/air.',year,'.nc',format='(a,i4.4,a3)')
      ncdf_GH_file=string(home_path+'/ncdf/hgt.',year,'.nc',format='(a,i4.4,a3)')
      ncdf_SH_file=string(home_path+'/ncdf/shum.',year,'.nc',format='(a,i4.4,a3)')
;
; Read NCEP files, checking that the full number of levels have been downloaded
      print,' '
      print,'Started reading global data...'
; Air Temperature
      read_global_data_file,ncdf_AT_file, 'air',Site_Lat,Site_Lon_360,Lev_AT,Lat_AT,LON_AT,TIM_AT,Global_Data_AT, Global_Data_scale_factor_AT, Global_Data_add_offset_AT,julday0
      print,'AT: Year, Nlev, Nlat, Nlong, Ntime =',year,N_Elements(Lev_AT),N_Elements(Lat_AT),N_Elements(Lon_AT),N_Elements(Tim_AT),format='(a,5i6)'
      If(N_Elements(Lev_AT) lt 17) then print,' Need 17 levels of AT data: found only ',N_Elements(Lev_AT)

; Geopotential Height
      read_global_data_file,ncdf_GH_file, 'hgt',Site_Lat,Site_Lon_360,Lev_GH,Lat_GH,LON_GH,TIM_GH,Global_Data_GH, Global_Data_scale_factor_GH, Global_Data_add_offset_GH,julday0
      print,'GH: Year, Nlev, Nlat, Nlong, Ntime =',year,N_Elements(Lev_GH),N_Elements(Lat_GH),N_Elements(Lon_GH),N_Elements(Tim_GH),format='(a,5i6)'
      If(N_Elements(Lev_GH) lt 17) then print,' Need 17 levels of GH data: found only ',N_Elements(Lev_GH)

; Specific Humidity
      read_global_data_file,ncdf_SH_file,'shum',Site_Lat,Site_Lon_360,Lev_SH,Lat_SH,LON_SH,TIM_SH,Global_Data_SH, Global_Data_scale_factor_SH, Global_Data_add_offset_SH,julday0
      print,'SH: Year, Nlev, Nlat, Nlong, Ntime =',year,N_Elements(Lev_SH),N_Elements(Lat_SH),N_Elements(Lon_SH),N_Elements(Tim_SH),format='(a,5i6)'
      If(N_Elements(Lev_SH) lt  8) then print,' Need  8 levels of SH data: found only ',N_Elements(Lev_SH)

      yearwas=year
   endif  ; year ne yearwas
;
; Take Julian date of time=0 from netcdf files instead of hardcoded value
; Added by D. Feist to make program work also with NCEP ncdf files after October 2014
; Replaced in following lines: julday(1,1,1,0) => julday0

   if Mod_Name ne Mod_Was then begin  ; Interpolate in long/latitude/time and Write Model
      nmod=nmod+1
      Site_Noon_Tim=(julday(month,day,year,12)-julday0)*24.0 - Site_Lon_180/15
      trilinear_interp1,Global_Data_AT, Global_Data_scale_factor_AT, Global_Data_add_offset_AT,Site_Lon_360,Lon_AT,Site_Lat,Lat_AT,Site_Noon_Tim,Tim_AT, Site_Noon_AT
      trilinear_interp1,Global_Data_GH, Global_Data_scale_factor_GH, Global_Data_add_offset_GH,Site_Lon_360,Lon_GH,Site_Lat,Lat_GH,Site_Noon_Tim,Tim_GH, Site_Noon_GH
      trilinear_interp1,Global_Data_SH, Global_Data_scale_factor_SH, Global_Data_add_offset_SH,Site_Lon_360,Lon_SH,Site_Lat,Lat_SH,Site_Noon_Tim,Tim_SH, Site_Noon_SH

; Convert Specific Humidity, a Wet Mass Mixing Ratio, to Dry Mole Fraction
      Site_Noon_H2ODMF=(28.964/18.02)*Site_Noon_SH/(1-Site_Noon_SH)

      Site_Noon_GH=Site_Noon_GH/1000   ; Convert m to km
      Site_Noon_TP=0.0                 ; Convert m to km

;  Write the NCEPxxxxxx.mod file into the appropriate sub-directory
      print,spectrum,'  ',Mod_Name
      write_mod,home_path+'/models/'+ext+Mod_Name,version,Site_Lat,Lev_AT,Site_Noon_AT,Site_Noon_GH,Site_Noon_TP,Site_Noon_H2ODMF
   endif

   Mod_Was=Mod_Name
endwhile
; End loop over spectra
close,unit
free_lun,unit
print,'Number of spectra = ',nspe
print,'Number of days/models = ',nmod
print,'Number of years = ',nyear
end
