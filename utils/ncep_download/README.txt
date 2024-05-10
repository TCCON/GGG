ncep_download - NCEP download script for Linux/Unix/MacOS etc.

This bash-based shell script downloads all the necessary NCEP data files from
the NOAA site automatically. The ncep_download script should run on all
Unix-like systems that provide bash and wget. All preprocessing steps like
selection for time or geographic region will be done on the NOAA server.
Therefore, the actual download can be very small and the result can directly be
used for running mod_maker.

Instructions

If run without parameters, the script will download NetCDF files for air
temperature (AT), geopotential height (GP), and specific humidity (SH).
The coverage will be global and the time period will be the most
recent day for which this data is available. You can customize ncep_download's
behaviour with several command line options.

Command line options

-S LAT    southern latitude limit (90S to 90N, default: 90S)
-N LAT    northern latitude limit (90S to 90N, default: 90N)
-W LON    western longitude limit (0-360E | 0-360W, default: 0E)
-E LON    eastern longitude limit (0-360E | 0-360W, default: 360E)
-f DATE   from-date (format: YYYYMMDD, default: same as to-date)
-t DATE   to-date (format: YYYYMMDD, default: last available date)
-s STRING site identifier (default: 'NONE')
-n NUM    site latitude [decimal degrees north] for mod_maker input file (default:
          none)
-e NUM    site longitude [decimal degrees east] for mod_maker input file (default:
          none)
-p PATH   optional path for output files (default: none)
-v        enable verbose output
-h        show help screen

Example

ncep_download -s je -n 50.91 -e 11.57 -S 49N -N 53N -W 9E -E 13E -f 20090801 >mod_maker.input

This will download data for the station je (Jena) located at 50.91° N, 11.57° E.
The latitude range is limited from 49° N to 53° N, the longitude range is
limited from 9° E to 13° E. The retrieval starts on August 1, 2009 and ends on
the last day for which data is available. The output of the script is directly
written into the mod_maker input file.

License

This software is free to use and modify for all members of the TCCON community
(everyone who is licensed to use the GGG software). If you make changes or find
bugs, please inform the author so everybody in the TCCON community can profit.

Contact

Dietrich Feist <dfeist@bgc-jena.mpg.de>
