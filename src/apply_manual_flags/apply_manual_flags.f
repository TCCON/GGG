c  Program apply_manual_flags.f
c
c  Purpose: to apply flags defined in a site-specific file
c  which indicate a specific date range is not to be released
c  publicly due to known issues with the instrument.
c
c
c  Input files:
c   runlog.vav.ada.oof
c   xx_manual_flags.dat (in $GGGPATH/tccon)
c
c  Output files:
c   runlog.vav.ada.flg.oof
c
      program apply_manual_flags

      implicit none
      include "../gfit/ggg_int_params.f"
      include "../comn/postproc_params.f"

      integer*4
     & fnbc, lnbc,
     & j, idum,
     & li, ldot, lcom,
     & lunr, lunrcsv, lunw, lunwcsv,
     & nhead, ncol, nheadcsv,
     & nflagged, man_flag, oof_flag, csv_flag,
     & flag_for_spectrum

      logical
     & csv_exists

      parameter (lunr=14,lunrcsv=15,lunw=16,lunwcsv=17)

      character 
     & version*64,
     & dl*1,
     & gggdir*(mpath),
     & input_fmt*40,
     & specname*(nchar), 
     & speccsv*(nchar),
     & inputfile*64,
     & outputfile*64,
     & csvfile*64,
     & outcsvfile*64,
     & hline*1800,
     & hlinecsv*1800,
     & nhcs*16,
     & dline*1800,
     & cdum*100

      version=
     & ' apply_manual_flags   Version 1.00  2021-07-26   JLL'
      write(*,*) version

      call get_ggg_environment(gggdir, dl)


c  Avoid unused parameter warnings
      idum = ncell
      idum = mvmode
      idum = mspeci
      idum = mrow_qc
      idum = mlev
      idum = mgas
      idum = mfilepath
      idum = mcolvsw
      idum = mcolvav
      idum = mauxcol
      idum = maddln
      idum = mcharhead
      cdum = countfmt


c  Get the input file interactively or from the command line

      if (iargc() == 0) then
        write(*,*)
     &'Enter name of input file (e.g. paIn_1.0lm.vav.ada.oof):'
        read(*,'(a)') inputfile
      elseif (iargc() == 1) then
        call getarg(1, inputfile)
      else
        stop 'Usage: $gggpath/bin/apply_manual_flags ooffile'
      endif

      li=lnbc(inputfile)
      if(inputfile(li-3:li).ne.'.oof') write(*,*)
     & 'Warning: input file is not of expected type (.oof)'

c  Create the output name and open the files
      ldot=index(inputfile, '.', .true.)
      outputfile=inputfile(:ldot)//'flg.'//inputfile(ldot+1:li)
c      print *, outputfile
      open(lunr,file=inputfile, status='old')
      open(lunw,file=outputfile,status='unknown')

c  Does the .oof.csv file exist? If so we need to copy that one too
      csv_exists = .false.
      csvfile=inputfile(:lnbc(inputfile))//'.csv'
      outcsvfile=outputfile(:lnbc(outputfile))//'.csv'
      inquire(file=csvfile, exist=csv_exists)

      if (csv_exists) then
         open(lunrcsv, file=csvfile, status='old')
         open(lunwcsv, file=outcsvfile, status='unknown')
      endif


c  Because the .oof file includes the qc and other extra header information,
c  the missing and format lines aren't in the expected place so we're going to
c  have to do this the hard way.
      read(lunr, *) nhead, ncol
      write(lunw, *) nhead+1, ncol
      write(lunw, '(a)') version

      if (csv_exists) then
c  Need to use '(a)' format rather than *, as commas will stop list-directed reads
c  Useful, but it's actually easier to read the line and split on the comma for
c  writing back out.      
         read(lunrcsv, '(a)') hlinecsv
         lcom = index(hlinecsv, ',')
         read(hlinecsv(:lcom-1), *) nheadcsv
         if (nheadcsv .ne. nhead) 
     &     stop '.oof and .oof.csv have different length headers'

         write(nhcs, *) nheadcsv+1
         write(lunwcsv, '(a)') 
     &    nhcs(fnbc(nhcs):lnbc(nhcs))//hlinecsv(lcom:lnbc(hlinecsv))

         write(lunwcsv, '(a)') version
      endif

      input_fmt = ''
      do j=2,nhead
         read(lunr, '(a)') hline
         write(lunw, '(a)') hline(:lnbc(hline))
         if (index(hline, 'format:') .gt. 0) then
            li = index(hline, ':')  
            input_fmt = hline(li+1:)
         endif

         if (csv_exists) then
            read(lunrcsv, '(a)') hlinecsv
            write(lunwcsv, '(a)') hlinecsv(:lnbc(hlinecsv))
         endif
      end do ! j=2,nhead

      if (index(input_fmt ,')') .lt. 1) then
         stop 'Did not find format line in .oof file header'
      endif


c  Check each spectrum and update flags as needed
      nflagged = 0
      do j=1,100000
         read(lunr, '(a)', end=10) dline
c  Flag is not included in format, assume it is three characters long.         
         read(dline(:3), '(i3)') oof_flag
         read(dline(4:), input_fmt) specname

         man_flag = flag_for_spectrum(specname, gggdir, dl, j.eq.1)

c  A flag < 0 indicates that no manual flag was specified. In that case,
c  we still need to add spaces to the beginning of the line to pad it 
c  (since lines with a manual flag will require more space). Otherwise
c  overwrite the existing flag assuming that it is the first three characters.
         if (man_flag .gt. 0) then
            oof_flag = oof_flag + 1000*man_flag
            nflagged = nflagged + 1
         endif
         
         write(lunw, '(i5,1x,a)') oof_flag, dline(5:lnbc(dline))
         
         if (csv_exists) then
            read(lunrcsv, '(a)', end=11) dline
            read(dline, *) csv_flag, speccsv
            if (speccsv(:lnbc(speccsv)) .ne. 
     &          specname(:lnbc(specname))) then
               write(*,*) 'Error on data row', j
               stop 'oof and oof.csv spectra mismatch'
            endif

c  Reuse dummy character var - it was only used before to avoid unused variable
c  warnings
            if (man_flag .gt. 0) then
               write(cdum, *) csv_flag + 1000*man_flag
            else
               write(cdum, *) csv_flag
            endif
            write(lunwcsv, '(a,a)')
     &        cdum(fnbc(cdum):lnbc(cdum)),
     &        dline(index(dline, ','):lnbc(dline)) 
         endif
      end do

11    stop 'oof and oof.csv have different numbers of rows'

10    write(*,*) 'apply_manual_flags done, flagged ',
     & nflagged, ' spectra.' 

      close(lunr)
      close(lunw)
      if (csv_exists) then
         close(lunrcsv)
         close(lunwcsv)
      endif
      end program


      function flag_for_spectrum(spectrum, gggdir, dl, first_call)
c  Get the manual flag for a given spectrum. 
c  Inputs:
c     spectrum: CHARACTER(*) - The spectrum name that we want to check for 
c        manual flags for
c
c  Returns:
c     INTEGER*4 - the flag value. Negative values indicate no flag was found.

         integer*4 flag_for_spectrum, mrow, jday_spec, lunf, 
     &    nhead, ncol, j, jmax, lc
         character site*2, site_was*2, datestr*8, dl*1, flag_file*128,
     &    hline*64, fmt*64    
         character(*) spectrum, gggdir
         logical first_call

         parameter(mrow=500, lunf=20)
         
         integer*4 y1, y2, m1, m2, d1, d2, 
     &    jday_start(mrow), jday_end(mrow), flag(mrow),
     &    fnbc, lnbc

         save site_was, jday_start, jday_end, flag, jmax

c  First figure out the site and date from the spectrum name - assumes CIT
c  naming convention: xxYYYYMMDD
         spectrum = spectrum(fnbc(spectrum):)
         site = spectrum(:2)
         datestr = spectrum(3:10)
         read(datestr, '(i4,i2,i2)') y1, m1, d1
         call julian(y1, m1, d1, jday_spec)

c         print *, first_call, site .ne. site_was
         if (first_call .or. site .ne. site_was) then
            flag_file = gggdir(:lnbc(gggdir))//dl//'tccon'//dl//site//
     &      '_manual_flagging.dat'

c  Open the flag file, if not present, tell the user to create it.
            write(*,*) 'Reading ', flag_file
            site_was = site
            open(lunf, file=flag_file, status='old', err=200)
            read(lunf, *) nhead, ncol

            do j=2,nhead
               read(lunf, '(a)') hline
               if (index(hline, 'format:') .gt. 0) then
                  lc = index(hline, ':')
                  fmt = hline(lc+1:)
               endif
            end do

c  Should now be past the header and ready to read in the date and flag arrays
            jmax = 0
            do j=1,mrow
               read(lunf, fmt, end=100) y1, m1, d1, y2, m2, d2, flag(j)
               call julian(y1, m1, d1, jday_start(j))
               call julian(y2, m2, d2, jday_end(j))

c  If there is a blank line at the end of the file, flag(j) will get read in as
c  zero, which causes this message to be printed confusingly. Checking for y1 > 0
c  ensure that the line wasn't blank.
               if ((flag(j) .lt. 1 .or. flag(j) .gt. 9)
     &              .and. y1 .gt. 0) then
                  write(*,*) 'Manual flags < 1 or > 9 not allowed.'
                  write(*,*) 'Will be replaced with 9 (other). Row=', j
                  flag(j) = 9
               endif

               if (y1 .gt. 0) jmax = jmax + 1
            end do   
            
            stop 'increase mrow'

100         close(lunf)
         endif ! if (first_call .or. site .ne. site_was)

c         print *, 'jday_spec=', jday_spec
c         print *, 'jday_start=', jday_start(:jmax)
c         print *, 'jday_end=',jday_end(:jmax)
c         print *, 'flag=', flag(:jmax)
         do j=1,jmax
            if (jday_spec .ge. jday_start(j) .and.
     &          jday_spec .le. jday_end(j)) then

               flag_for_spectrum = flag(j)
c               print *, spectrum, ' manual flag ', flag(j)
               return
            endif
         end do

c  Did not find any flags for this day
c         print *, spectrum, ' no manual flag'
         flag_for_spectrum = -9999
         return

200      write(*,*) flag_file//
     & ' not found - required for apply_manual_flags.'
         stop
      end function
