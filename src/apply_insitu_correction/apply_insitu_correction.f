c  Program: apply_insitu_correction.f
c
c  Purpose: To apply the (airmass-independent) in situ correction to the data
c  in the selected runlog.vav.ada file
c  
c  yout(i,j) = yin(i,j) / AICF(j)
c
c  where
c   yin(i,j) is the DMF of the j'th gas from the i'th spectrum
c   yout(i,j) is the DMF of the j'th gas from the i'th spectrum
c   AICF(j) is the Airmass-Independent Correction Factor for the j'th gas
c   which is read from the file ~/ggg/tccon/corrections.dat
c
c  Input Files:
c       runlog.vav.ada 
c       corrections.dat  
c
c  Output Files:
c       runlog.vav.ada.aia
c       
      implicit none
      include "../gfit/ggg_int_params.f"
      include "../comn/postproc_params.f"

      integer*4 lunr,luns,lunw,ncoml,ncolvav,kcolvav,icol,j,
     & kgas,lnbc,irow,naux,ngas,nrow,li,k,ncolcorr,naddn,has_err
      parameter (lunr=14,luns=15,lunw=16)
      real*8 yrow(mcolvav),ymiss,
     & aicf(mgas),aicf_err(mgas),cf(mcolvav)

      character header*1800,headarr(mcolvav)*40,gasname(mgas)*20,dl*1,
     & gggdir*(mpath),inputfile*64,outputfile*64, version*62,gaserr*32
      character output_fmt*40, input_fmt*40, specname*(nchar), c1*1,
     & addn_lines(maddln)*(mcharhead)
      
      integer specflag,idum
      specflag=0

      idum=mauxcol   ! Avoid compiler warning (unused parameter)
      idum=mcolvsw   ! Avoid compiler warning (unused parameter)
      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mlev      ! Avoid compiler warning (unused parameter)
      idum=mrow_qc   ! Avoid compiler warning (unused parameter)
      idum=mspeci    ! Avoid compiler warning (unused parameter)
      idum=mvmode    ! Avoid compiler warning (unused parameter)
      idum=ncell     ! Avoid compiler warning (unused parameter)

      version=
     & ' apply_insitu_correction   Version 1.38  2020-03-20   GCT,JLL'
      write(*,*) version

      call get_ggg_environment(gggdir, dl)

c  Open the corrections.dat file and read in the Airmass-Dependent
c  Correction Factors (ADCF) and the Airmass-Independent Correction
c  Factors (AICF) for each gas. Only the former is used by this prog.
      open(luns,file=gggdir(:lnbc(gggdir))//'tccon'//dl
     & //'corrections_insitu_postavg.dat', status='old')
      read(luns,*)ncoml,ncolcorr
      do k=2,ncoml
         read(luns,*)
      end do

      if(ncolcorr .eq. 2) then
        has_err = 0
      elseif(ncolcorr .eq. 3) then
        has_err = 1
      else
         write(*,*)'ncol=',ncolcorr
         stop 'Unrecognized NCOLCORR value'
      end if


      if(has_err .eq. 0) then
         do k=1,mgas
            read(luns,*,end=88) gasname(k),aicf(k)
         end do
      elseif(has_err .eq. 1) then
         do k=1,mgas
            read(luns,*,end=88) gasname(k),
     &     aicf(k),aicf_err(k)
         end do
      else
      endif
      stop 'increase parameter MGAS'
88    ngas=k-1

      if (iargc() == 0) then
         write(*,*)'Enter name of input file (e.g. paIn_1.0lm.vav.ada):'
         read(*,'(a)') inputfile
      elseif (iargc() == 1) then
         call getarg(1, inputfile)
      else
          stop 'Usage: $gggpath/bin/apply_insitu_correction adafile'
      endif
      li=lnbc(inputfile)
      if(inputfile(li-3:li).ne.'.ada') write(*,*)
     & 'Warning: input file is not of expected type (.ada)'
      outputfile=inputfile(:li)//'.aia'
      open(lunr,file=inputfile, status='old')
      open(lunw,file=outputfile,status='unknown')

c  Read the header of the .ada file and figure out the
c  mapping between the gases in the corrections.dat
c  and those in the .vav.ada file header. Prepend this 
c  program's version and append the insitu corrections
c  to the header.
      call read_postproc_header(lunr, ncoml, ncolvav, nrow, naux,
     & ymiss, input_fmt, addn_lines, naddn)
      read(lunr,'(a)') header

      do j=naddn,1,-1
        addn_lines(j+1) = addn_lines(j)
      end do
      addn_lines(1) = version
      naddn = naddn + 1

      write(addn_lines(naddn+1),'(a,i2,1x,i1)')
     & 'Airmass-Independent/In-Situ Correction Factors:', ngas, has_err
      do k=1,ngas
         if(has_err .eq. 0) write(addn_lines(naddn+1+k),'(a,f9.4)') 
     & gasname(k), aicf(k)
         if(has_err .eq. 1) write(addn_lines(naddn+1+k),'(a,2f9.4)') 
     & gasname(k), aicf(k), aicf_err(k)
      end do
      naddn = naddn + ngas + 1

      if(ncolvav.gt.mcolvav) stop 'increase mcolvav'

c  Write the header of the output file.
      output_fmt=input_fmt
      if (index(header,'spectrum') .gt. 0) specflag=1
      if(specflag .eq. 1) output_fmt(6:7) = '1x'

      call write_postproc_header(lunw, ncolvav, nrow, naux,
     & ymiss, output_fmt, addn_lines, naddn, 0)

      write(lunw,'(a)') header(:lnbc(header))
      call substr(header,headarr,mcolvav,kcolvav)
      if(kcolvav.ne.ncolvav ) stop 'ncolvav/kcolvav mismatch'
      do icol=1,ncolvav
         cf(icol)=1.0
         do kgas=1,ngas
c         write(*,*)kgas,icol,gasname(kgas),headarr(icol),aicf(kgas),cf(icol)
           gaserr=gasname(kgas)(:lnbc(gasname(kgas)))//'_error'
           if( headarr(icol) .eq. gasname(kgas) ) cf(icol)=aicf(kgas)
           if( headarr(icol) .eq. gaserr        ) cf(icol)=aicf(kgas)
         end do
         if(icol.gt.naux+specflag) write(*,'(i4,1x,a16,1x,f9.4)')
     &   icol,headarr(icol),cf(icol)
      end do

c      if (specflag .eq. 1) then
c         output_fmt='(a35,f14.8,NNf13.5,200(1pe12.4))'
c         write(output_fmt(12:13),'(i2.2)') naux-1
c      else
c         output_fmt='(f14.8,NNf13.5,200(1pe12.4))'
c         write(output_fmt(8:9),'(i2.2)') naux-2
c      endif

c  Read each day of data into memory and divide XGas values by the
c  appropriate correction factors.
      do irow=1,9999999
         if (specflag .eq. 1) then
            read(lunr,input_fmt,end=99) specname,c1,
     & (yrow(j),j=2,ncolvav)
         else
            read(lunr,input_fmt,end=99) c1,(yrow(j),j=1,ncolvav)
         endif
         do k=naux+1,ncolvav-1,2
              if(yrow(k).lt.ymiss) then  ! Dont correct missing values
                 yrow(k)=yrow(k)/cf(k)
                 yrow(k+1)=yrow(k+1)/cf(k)
              endif  !  (yrow(k).lt.ymiss) then
         end do
         if (specflag .eq. 1) then
            write(lunw,output_fmt) specname,(yrow(j),j=2,ncolvav)
         else
            write(lunw,output_fmt) c1,(yrow(j),j=1,ncolvav)
         endif
      end do         ! do irow=1,9999999
      stop ' irow exceeded 9999999'
99    close (lunr)
      close(lunw)
      stop
      end
