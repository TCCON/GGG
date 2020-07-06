c  subroutine: prepare_oof_output.f
c
      subroutine prepare_oof_output(inputfile,inputlun,
     & outfmt1, oof_flag, irow,nrow_qc, 
     & ada_ncorr, ada_has_err, ada_gasname, adcf, adcf_err,
     & aia_ncorr, aia_has_err, aia_gasname, aicf, aicf_err,
     & vmin, vmax, pindex, 
     & eflag, kflag, flag,
     & nrow, ncol, mchar, scale,rsc,ofmt, headout,naux)
     

      implicit none
      include "../gfit/ggg_int_params.f"
      include "params.f"
      include "../comn/postproc_params.f"

      integer*4
     & ncoml,ncol,kcol,icol,j,lg,
     & lnbc,nrow,li,k,krow_qc,irow,lof0,lof1,
     & eflag,wrow_flag,wsp_flag,naux,
     & mchar,lh,le,klat,klong,kzobs,ncol_written,
     & jj,ncoml_qc,ncol_qc,nrow_qc,wcol_flag,
     & naddn,ada_ncorr,aia_ncorr,mwin,ada_has_err,aia_has_err
      integer*4 flag(mrow_qc),pindex(mcol),kflag(mrow_qc),
     & oof_flag(mrow)

      parameter (mwin=100)

      character
     & cdum*1,
     & fmtdum*40,
     & dl*1,
     & gggdir*(mpath),
     & headarr(mcol)*20,
     & parname(mrow_qc)*20,
     & inputfile*80,
     & outfmt0*512,
     & outfmt1*512,
     & ofmt(mrow_qc)*4,
     & temp*20,
     & fmt(mrow_qc)*4,
     & unit(mrow_qc)*6,
     & headout*8000,
     & specfmt*3,
     & ada_gasname(mwin)*20,
     & aia_gasname(mwin)*20,
     & addn_lines(maddln)*(mcharhead),
     & ssss*800
      real*4 scale(mrow_qc),rsc(mrow_qc),
     & vmin(mrow_qc),vmax(mrow_qc)
      real*8 adcf(mwin), adcf_err(mwin), aicf(mwin), aicf_err(mwin),
     & fdum
      integer mrow_qc_return, mcol_return
      integer inputlun
      integer i,idum

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)

      idum=lun_qc   ! Avoid compiler warning (unused parameter)
      idum=mcol     ! Avoid compiler warning (unused parameter)
      idum=mrow     ! Avoid compiler warning (unused parameter)
      idum=mluns    ! Avoid compiler warning (unused parameter)
      idum=maddln   ! Avoid compiler warning (unused parameter)
      idum=mcharhead! Avoid compiler warning (unused parameter)
      idum=lun_mul  ! Avoid compiler warning (unused parameter)
      idum=lun_rlg  ! Avoid compiler warning (unused parameter)

      idum=oof_flag(1)

      cdum=gfit_version(1:1)  ! Avoid compiler warning (unused parameter)
      cdum=gsetup_version(1:1) ! Avoid compiler warning (unused parameter)
      cdum=tllsum(1:1)   ! Avoid compiler warning (unused parameter)
      cdum=solarsum(1:1) ! Avoid compiler warning (unused parameter)
      cdum=csformat(1:1) ! Avoid compiler warning (unused parameter)
      cdum=header(1:1)   ! Avoid compiler warning (unused parameter)

      irow=1
      mrow_qc_return = mrow_qc
      mcol_return = mcol
      eflag=0
      
      do i = 1, mrow_qc
         kflag(i) = 0
         vmin(i) = 0
         vmax(i) = 0
      end do
      
      mchar=0
      klat=0
      klong=0
      kzobs=0
      wrow_flag=1
      wsp_flag=1

      call get_ggg_environment(gggdir, dl)
      lg=lnbc(gggdir)

      li=lnbc(inputfile)
c     if(inputfile(li-3:li).ne.'.aia'.or.
c    & inputfile(li-3:li).ne.'.gaa') write(*,*)
c    & 'Warning: input file is not of expected type (.aia or .gaa)'

c  Open the Quality Control (QC) file and read in the information,
c  to find out how many columns there are going to be in the .oof output files
c     write(*,*)gggdir(:lg)//'/tccon/'//inputfile(1:2)//'_qc.dat'
      open(lun_qc,file=
     & gggdir(:lg)//'tccon'//dl//inputfile(1:2)//'_qc.dat',
     & status='old')
c     write(*,*)gggdir(:lg)//'/tccon/'//inputfile(1:2)//'_qc.dat'
      read(lun_qc,*)ncoml_qc,ncol_qc,nrow_qc
      if(nrow_qc.gt.mrow_qc) stop 'nrow_qc > mrow_qc'
      read(lun_qc,*)wrow_flag   ! 0/1   whether to skip/write the out-of-range data records
      read(lun_qc,*)wsp_flag    ! 0/1   whether to skip/write the spectrum names (if present in .aia file)
      do k=4,ncoml_qc
         read(lun_qc,'(a)') header
      end do
      ncol_written=1   ! Always include "flag" in column #1
      wrow_flag=1      ! Always include out-of-range records in the eof
      wsp_flag=1       ! Always include the spectrum names in the eof
      do krow_qc=1,mrow_qc
         read(lun_qc,*,end=77) temp,wcol_flag
         ncol_written=ncol_written+wcol_flag
      end do
77    close(lun_qc)
      if(krow_qc-1.ne.nrow_qc)
     &  stop 'misreading xx_qc.dat file: krow_qc > nrow_qc'

c  Read input file and start writing output files
c       write(*,*) inputfile
      open(inputlun,file=inputfile, status='old')
      call read_postproc_header(inputlun, ncoml, ncol, nrow, naux,
     & fdum, fmtdum, addn_lines, naddn)
      if(ncol.gt.mcol) stop 'increase mcol'
      read(inputlun,'(a)') header ! column headers
      call substr(header,headarr,mcol,kcol)
      if(kcol.ne.ncol) then
        write(*,*) kcol, ncol, inputfile
        stop 'ncol/kcol mismatch'
      end if
      if (index(header,'spectrum') .gt. 0) mchar=1

c  Parse the header looking for the airmass & insitu correction factors
      do j=1,naddn
       if (index(addn_lines(j), 'Airmass-Dependent') .gt. 0) then
        icol = index(addn_lines(j), ':')
        read(addn_lines(j)(icol+1:), *) ada_ncorr, ada_has_err
        do k=1,ada_ncorr
          if(ada_has_err .eq. 1) then
            read(addn_lines(j+k),*) ada_gasname(k), adcf(k), adcf_err(k)
          else
            read(addn_lines(j+k),*) ada_gasname(k), adcf(k)
          end if
        end do
        goto 14
       end if 
      end do
      
      stop 'Could not find ADCF values in .aia file'

 14   continue

      do j=1,naddn
       if (index(addn_lines(j), 'In-Situ') .gt. 0) then
        icol = index(addn_lines(j), ':')
        read(addn_lines(j)(icol+1:), *) aia_ncorr, aia_has_err
        do k=1,aia_ncorr
          if (aia_has_err .eq. 1) then
            read(addn_lines(j+k),*) aia_gasname(k), aicf(k), aicf_err(k)
          else
            read(addn_lines(j+k),*) aia_gasname(k), aicf(k)
          end if
        end do
        goto 19
       end if
      end do

      stop 'Could not find AICF values in .aia file'

 19   continue

      open(lun_qc,file=
     & gggdir(:lg)//'tccon'//dl//inputfile(1:2)//'_qc.dat',
     & status='old')
      read(lun_qc,*)ncoml_qc,ncol_qc,nrow_qc
      do k=2,ncoml_qc
         read(lun_qc,'(a)') header
c        if(k.eq.ncoml_qc) header=' #'//header
      end do

      do krow_qc=1,nrow_qc
         read(lun_qc,'(a)') ssss
         read(ssss,*) parname(krow_qc),flag(krow_qc),scale(krow_qc),
     &   fmt(krow_qc),unit(krow_qc),
     &   vmin(krow_qc),vmax(krow_qc)
c        write(*,*)'unit: ',unit(krow_qc)(2:lnbc(unit(krow_qc))-1)
         if (unit(krow_qc)(1:1).eq.'(') then
             unit(krow_qc)=unit(krow_qc)(1:lnbc(unit(krow_qc))-1)
             unit(krow_qc)(1:1)='_'
         endif
      end do
      close(lun_qc)

      do krow_qc=1,nrow_qc
         rsc(krow_qc)=scale(krow_qc)
      end do
c  Read the header of the .aia file and figure out the
c  mapping between the gases in the corrections.dat
c  and those in the .vav file header
      do icol=1,ncol
         pindex(icol)=0
         temp=headarr(icol)
         if(temp(1:4).eq.'xao2') temp='x'//temp(3:)
         if(temp(1:4).eq.'xbo2') temp='x'//temp(3:)
         do krow_qc=1,nrow_qc
            if( temp .eq. parname(krow_qc) ) pindex(icol)=krow_qc
         end do

c         if(pindex(icol).eq.0) write(*,*)
c     &   ' Parameter missing from QC file: '//headarr(icol)

         if(pindex(icol).eq.0) then
           if(index(headarr(icol),'spectrum').eq.0) then
           write(*,*) ' Parameter missing from QC file: '//headarr(icol)
           endif
         endif

         if(headarr(icol).eq.'lat') klat=icol
         if(headarr(icol).eq.'long') klong=icol
         if(headarr(icol).eq.'zobs') kzobs=icol
      end do

      headout=' flag  spectrum'
      outfmt0='(i3'
      write(specfmt,'(a,i2)')'a',nchar
      outfmt1='(i3,1x,'//specfmt
      lof0=0
      lof1=0
      jj=0
      do icol=1+mchar,ncol
         krow_qc=pindex(icol)
         if(flag(krow_qc).ge.1) then
            lof0=lnbc(outfmt0)
            lof1=lnbc(outfmt1)
            jj=jj+1
            ofmt(jj)=fmt(krow_qc)
            outfmt0=outfmt0(:lof0)//',1x,'//fmt(krow_qc)
            outfmt1=outfmt1(:lof1)//',1x,'//fmt(krow_qc)
            le=index(headarr(icol),'_error')
            lh=lnbc(headarr(icol))
            if(le.eq.0) then
               headout=headout(:lnbc(headout))//'  '//
     &         headarr(icol)(:lh)//unit(krow_qc)
            else
               headout=headout(:lnbc(headout))//'  '//
     &         headarr(icol)(:le-1)//
     &         unit(krow_qc)(:lnbc(unit(krow_qc)))//'_error'
            endif
         endif
      end do
      outfmt0(lof0+9:lof0+9)=')'
      outfmt1(lof1+9:lof1+9)=')'
      end
