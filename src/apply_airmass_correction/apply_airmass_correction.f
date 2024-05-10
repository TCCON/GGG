C  Program: apply_airmass_dependence.f
c
c  Purpose: To apply the airmass correction to the data in the
c  selected .vav filE
c  
c  yout(i,k) = yin(i,k) / [ 1 + ADCF(k)*SBF(i) ]
c
c  where
c   yin(i,k) is the column measured for the k'th gas from the i'th spectrum
c   yout(i,k) is the DMF of the k'th gas from the i'th spectrum
c   ADCF(j) is the Airmass-Dependent Correction Factor for the k'th gas
c   SBF(i) is the Symmetric Basis Function (symmetric about solar noon)
c   SBF(i) = ((SZA(i)+13)/(90+13))**3 - ((45+13)/(90+13))**3
c   SZA(i) is the Solar Zenith Angle of the i'th spectrum.
c   SBF is zero at SZA=45 deg, -0.122 at SZA=0, and 0.875 at SZA=77 deg.
c   Over the range of SZAs (0-80 deg) encounted at a site, SBF changes by ~1,
c   so ADCF is a good estimate of the size of the airmass correction.
c
c
c  Input Files:
c       runlog.vav 
c       corrections.dat  
c
c  Output Files:
c       runlog.vav.ada
c       
      implicit none
      include "../gfit/ggg_int_params.f"
      include '../comn/postproc_params.f'

      integer*4 lunr,luns,lunw,ncoml,ncolvxx,
     & kcolvxx,jcol,j,nhead,naddn,
     & kgas,ko2,kh2o,kluft,ksza,lnbc,irow,naux,ngas,nrow,li,k,ncolcorr,
     & has_err, has_gp
      parameter (lunr=14,luns=15,lunw=16,nhead=20)
c      parameter (default_g=13.0, default_p=3.0)

      real*8 default_g, default_p
      real*8 yrow(mcolvsw),ymiss,
     & adcf(mgas),adcf_err(mgas),gval_gas(mgas),pval_gas(mgas),
     & cf(mcolvsw),gval_col(mcolvsw),pval_col(mcolvsw),gfu,sbf,vc_air
      character
     & dl*1,
     & gggdir*(mpath),
     & version*62,
     & specname*(nchar),
     & header*5000,
     & headarr(mcolvsw)*(nhead),
     & gasname(mgas)*(nhead),
     & gaserr*(nhead+6),
     & inputfile*64,
     & outputfile*64,
     & output_fmt*40,
     & col1*1,
     & addn_lines(maddln)*(mcharhead),
     & corr_file*40

      integer*4 specflag,idum 
      logical isavg  ! set to .true. if corrected averaged Xgas value, not individual windows
      logical isclose_d  ! function used to test is values are within
                         ! floating point error of fill value
      specflag=0

      idum=mcolvav   ! Avoid compiler warning (unused parameter)
      idum=mauxcol   ! Avoid compiler warning (unused parameter)                                                       
      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mlev      ! Avoid compiler warning (unused parameter) 
      idum=mrow_qc   ! Avoid compiler warning (unused parameter)
      idum=mspeci    ! Avoid compiler warning (unused parameter)
      idum=mvmode    ! Avoid compiler warning (unused parameter)
      idum=ncell     ! Avoid compiler warning (unused parameter)

      version=
     &' apply_airmass_correction   Version 1.38 2020-12-16   GCT,JLL'
      write(*,*) version
      call get_ggg_environment(gggdir, dl)
      ko2=0
      kh2o=0
      kluft=0
      ksza=0
      default_g=13.0
      default_p=3.0


      if (iargc() == 0) then
         write(*,*)'Enter name of input file (e.g. paIn_1.0lm.vav):'
         read(*,'(a)') inputfile
      elseif (iargc() == 1) then
         call getarg(1, inputfile)
      else
         stop 'Usage: $gggpath/bin/apply_airmass_correction vavfile'
      endif

c  Figure out if working with averaged Xgas values or individual windows
      li=lnbc(inputfile)
      if (inputfile(li-3:li) .eq. '.vav') then
        isavg = .true.
        corr_file = 'corrections_airmass_postavg.dat'
      elseif (inputfile(li-3:li) .eq. '.vsw') then
        isavg = .false.
        corr_file = 'corrections_airmass_preavg.dat'
      elseif (inputfile(li-1:li) .eq. 'av') then
        isavg = .true.
        write(*,*) ' Warning: input file is not of expected type'//
     &  '(.vav or .vsw), assuming averaged quantities'
      else
        isavg = .false.
        write(*,*) ' Warning: input file is not of expected type'//
     &  '(.vav or .vsw), assuming unaveraged quantities'
      end if
      outputfile=inputfile(:li)//'.ada'

c  Open the corrections.dat file and read in the Airmass-Dependent
c  Correction Factors (ADCF) for each gas or window, depending if
c  we're working on a .vsw or .vav file
      write(*,*) 'Reading corrections from ', corr_file
      open(luns,file=gggdir(:lnbc(gggdir))//'tccon'//dl
     & //corr_file(:lnbc(corr_file)), status='old')
      read(luns,*)ncoml,ncolcorr
      do k=2,ncoml
         read(luns,*)
      end do

      
      if (ncolcorr .eq. 2) then
        has_err = 0
        has_gp = 0
      elseif (ncolcorr .eq. 3) then
        has_err = 1
        has_gp = 0
      elseif (ncolcorr .eq. 5) then
        has_err = 1
        has_gp = 1
      else
         write(*,*) 'ncolcorr=',ncolcorr
         stop 'Unrecognized ncolcorr value'
      endif

c  JLL: Make sure all parameters are initialized to their
c  default values
      do k=1,mgas
         adcf(k) = 0.0
         adcf_err(k) = 0.0
         gval_gas(k) = default_g
         pval_gas(k) = default_p
      end do

      if(has_err .eq. 1 .and. has_gp .eq. 1) then
         do k=1,mgas
            read(luns,*,end=88) gasname(k),
     &   adcf(k),adcf_err(k),gval_gas(k),pval_gas(k)
         end do
      elseif(has_err .eq. 1) then
         do k=1,mgas
            read(luns,*,end=88) gasname(k),
     &   adcf(k),adcf_err(k)
         end do
      else
         do k=1,mgas
            read(luns,*,end=88) gasname(k),adcf(k)
         end do
      endif
      stop 'increase parameter MGAS'
88    ngas=k-1

      open(lunr,file=inputfile, status='old')
      open(lunw,file=outputfile,status='unknown')

c  Read the header of the .vav or .vsw file
      call read_postproc_header(lunr, ncoml, ncolvxx, nrow, naux,
     & ymiss, output_fmt, addn_lines, naddn)
      read(lunr,'(a)') header

c  Prepend this program's version to the additional lines, then 
c  append all the airmass dependent correction factors
      do j=naddn,1,-1
         addn_lines(j+1) = addn_lines(j)
      end do
      naddn = naddn + 1
      addn_lines(1) = version
      write(addn_lines(naddn+1),'(a,i2,1x,i1)') 
     & ' Airmass-Dependent Correction Factors: ', ngas, has_err+4
      do k=1,ngas
         if(has_err .eq. 1) then
           write(addn_lines(k+naddn+1),'(a,2(1x,f9.4),2(1x,f6.2))') 
     & gasname(k), adcf(k), adcf_err(k), gval_gas(k), pval_gas(k)
         else
           write(addn_lines(k+naddn+1),'(a,1x,f9.4,2(1x,f6.2))') 
     & gasname(k), adcf(k), gval_gas(k), pval_gas(k)
         end if
      end do
      naddn = naddn + ngas + 1

      
c  Check that we have enough columns to read data into
      if(ncolvxx.gt.mcolvsw) stop 'increase mcolvsw'

      call substr(header,headarr,mcolvsw,kcolvxx)
      if (index(header,'spectrum') .gt. 0) specflag=1

c  Go ahead and write the output file's header
      call write_postproc_header(lunw, ncolvxx, nrow, naux,
     & ymiss, output_fmt, addn_lines, naddn, 0)

      do j=naux+1,ncolvxx
         headarr(j)='x'//headarr(j)(:nhead-1)
      end do
      write(lunw,'(400a20)') (headarr(j),j=1,ncolvxx)
      if(kcolvxx.ne.ncolvxx ) stop 'ncolvxx/kcolvxx mismatch'

c  Figure out the mapping between the gases in the corrections.dat
c  and those in the .vav file header.  Avoid assuming that they
c  are all present or in the same order.
      do jcol=1,ncolvxx
         cf(jcol)=0.0
         gval_col(jcol)=default_g
         pval_col(jcol)=default_p
         do kgas=1,ngas
            if(headarr(jcol).eq.gasname(kgas)) then
               cf(jcol)=adcf(kgas)
               gval_col(jcol)=gval_gas(kgas)
               pval_col(jcol)=pval_gas(kgas)
            endif

            gaserr=gasname(kgas)(:lnbc(gasname(kgas)))//'_error'
            if(headarr(jcol).eq.gaserr) then
               cf(jcol)=adcf(kgas)
               gval_col(jcol)=gval_gas(kgas)
               pval_col(jcol)=pval_gas(kgas)
            endif
         end do
         if (isavg) then
            if(headarr(jcol) .eq. 'xo2') ko2=jcol
            if(headarr(jcol) .eq. 'xluft') kluft=jcol
         else
            if(headarr(jcol) .eq. 'xo2_7885') ko2=jcol
            if(headarr(jcol) .eq. 'xluft_6146') kluft=jcol
         end if
c         if(headarr(jcol) .eq. 'xh2o') kh2o=jcol
         if(headarr(jcol) .eq. 'asza' .or. 
     &     headarr(jcol) .eq. 'solzen' ) ksza=jcol
         if(jcol.gt.naux) write(*,'(i3,f8.4,2(1x,f6.2),2x,a)')
     &   jcol,cf(jcol),gval_col(jcol),pval_col(jcol),headarr(jcol)
      end do

      if(ko2.eq.0) stop ' o2 column not found'
      if(kluft.eq.0) stop ' luft column not found'
c      if(kh2o.eq.0) then
c         write(*,*) 'Warning: h2o column not found'
c         write(*,*) 'Assuming no h2o'
c      endif
      if(ksza.eq.0) stop ' asza column not found'

c      if (specflag .eq. 1) then
c         output_fmt='(a38,f14.8,NNf13.5,200(1pe12.4))'
c         write(output_fmt(12:13),'(i2.2)') naux-1-specflag
c         write(output_fmt(20:22),'(i3.3)') ncolvxx-naux
c      else 
c         output_fmt='(f14.8,NNf13.5,MMM(1pe12.4))'
c         write(output_fmt(8:9),'(i2.2)') naux-1-specflag
c         write(output_fmt(16:18),'(i3.3)') ncolvxx-naux
c      endif

c  Read each day of data into memory.
      do irow=1,9999999
         if (specflag .eq. 1) then 
            read(lunr,output_fmt,end=99) specname, col1, 
     & (yrow(j),j=2,ncolvxx)
         else
c  JLL 2020-03-23: note that behavior with specflag==.false. not tested
c  since change of output_fmt(6:7) from a1 to 1x removed
            read(lunr,*,end=99) (yrow(j),j=1,ncolvxx)
         endif
         vc_air=yrow(ko2)/0.2095
c  Compute dry air instead of wet air if kh2o exists
c         if(kh2o.gt.0) yrow(kluft)=yrow(kluft)-yrow(kh2o)*18.02/28.964
         if(yrow(ko2).eq.0.0) then
            write(*,'(a,i6,a)')'Warning: O2_column=0 for spectrum',irow,
     &      '  (lamp run?)  Your output file will contain Inf'
         endif
         do k=naux+1,ncolvxx-1,2
            sbf=((yrow(ksza)+gval_col(k))/(90+gval_col(k)))**pval_col(k)
     &      -((45.+gval_col(k))/(90+gval_col(k)))**pval_col(k)  ! Symmetric Basis Function
            if(.not. isclose_d(yrow(k), ymiss)) then
               if(k.eq.ko2 .or. .not. isavg) then
c  JLL 2020-06-04: error propagation is handled differently 
                  gfu=yrow(k+1)   ! Gas Fractional Uncertainty
c                  write(*,*) 'Not adding O2 error for column ', k ! debugging only
               else
                  gfu=sqrt(yrow(k+1)**2+
     &                    (yrow(k)*yrow(ko2+1)/yrow(ko2))**2)
               endif
               yrow(k)=yrow(k)/vc_air/(1+cf(k)*sbf)   ! apply airmass correction
               yrow(k+1)=gfu/vc_air/(1+cf(k)*sbf)
            endif  ! if(.not. isclose_d(yrow(k), ymiss))
         end do
         if (specflag .eq. 1) then
            write(lunw,output_fmt) specname, col1, (yrow(j),j=2,ncolvxx)
         else
c  JLL 2020-03-23: note that behavior with specflag==.false. not tested
c  since change of output_fmt(6:7) from a1 to 1x removed
            write(lunw,output_fmt) ' ',(yrow(j),j=1,ncolvxx)
         endif
      end do      ! do irow=1,9999999
      stop ' irow exceeded 9999999'
99    close (lunr)
      close (lunw)
      stop
      end
