c  subroutine read_oneline_oof.f
c

      subroutine read_oneline_oof(inputlun, oof_out,
     & outfmt1, oof_flag, irow, lun_esf, yesf,ncol_esf,nrow_qc,
     & vmin, vmax, pindex,
     & nflag, eflag, kflag, flag,
     & ncol, mchar, scale, rsc, ofmt,gaa_naux,yrow1,specname
     & )

      

      implicit none
      include "../gfit/ggg_int_params.f"
      include "params.f"

      integer*4 ncol, irow, gaa_naux, irow_qc, nrow_qc,idum
      integer*4 nflag, eflag, inputlun, mchar, lnbc,lnblnk
      integer*4 oof_flag(mrow),kcol
      real*4 vmin(mrow_qc), vmax(mrow_qc)
      integer*4 kflag(mrow_qc), flag(mrow_qc)
      integer*4 pindex(mcol),ncol_esf
      character*512 outfmt1
      character*800 oof_out
      character ofmt(mrow_qc)*4, cdum*1
      real*4 scale(mrow_qc),yesf(mcol),rsc(mrow_qc)

C local
      character*(nchar) specname
      integer*4 kmax, krow_qc, nco
      integer*4 j,jj,icol,jcol_gaa,lun_esf
      real*4 yrow(mcol),yrow1(mcol),apesf
      real*8 wlimit
      real*4 dmax, dev

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=nrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)
      idum=lun_qc   ! Avoid compiler warning (unused parameter)
      idum=lun_rlg  ! Avoid compiler warning (unused parameter)
      idum=lun_mul  ! Avoid compiler warning (unused parameter)
      idum=lnblnk(header)


      idum=lun_qc   ! Avoid compiler warning (unused parameter)
      idum=mcol     ! Avoid compiler warning (unused parameter)
      idum=mrow     ! Avoid compiler warning (unused parameter)
      idum=mluns    ! Avoid compiler warning (unused parameter)
      idum=lun_mul  ! Avoid compiler warning (unused parameter)
      idum=lun_rlg  ! Avoid compiler warning (unused parameter)

      cdum=gfit_version(1:1)  ! Avoid compiler warning (unused parameter)
      cdum=gsetup_version(1:1) ! Avoid compiler warning (unused parameter)
      cdum=tllsum(1:1)   ! Avoid compiler warning (unused parameter)
      cdum=solarsum(1:1) ! Avoid compiler warning (unused parameter)
      cdum=csformat(1:1) ! Avoid compiler warning (unused parameter)
      cdum=header(1:1)   ! Avoid compiler warning (unused parameter)
      
      nflag=0
      if (mchar .eq. 1) then
         read(inputlun,*) specname, (yrow(j),j=1+mchar,ncol)
      else
         read(inputlun,*) (yrow(j),j=1,ncol)
      endif

      do j=1,ncol
         jj=j*2-1
         yrow1(j) = yrow(jj+gaa_naux)
      enddo


c  If new day, read in the daily error scale factors and compute
c  new scale factors (RSC) as weighted averages of the a priori
c  ESF factors from the pa_qc.dat file, and the daily values.
c  A priori ESF values are the ratio of the xx_error/xxa scale factors
c  read in from the pa_qc.dat file, with 100% uncertainties assumed.
c  Starts at right-hand column (usually hcl_error) and works to left.
c     write(*,*)'yrow(2)=',yrow(2),' yesf(1)=',yesf(1),'ncol=',ncol
      if(yrow(2).gt.yesf(1)) then
c      write(*,*) 'reading...'
         read(lun_esf,*) (yesf(j),j=1,ncol_esf)
c        write(*,*)'yesf=',yesf(:ncol_esf)
         jcol_gaa = ncol
         do kcol=ncol_esf,4,-2
            irow_qc = pindex(jcol_gaa)
            apesf=scale(irow_qc)/scale(irow_qc-1)  ! xx_error/xx
c           write(*,*)'pindex(jcol_gaa)=',pindex(jcol_gaa),
c    &       ' jcol_gaa=',jcol_gaa,' irow_qc=',irow_qc,' apesf=',apesf
c            write(*,*)kcol,irow_qc,gaa_headarr(irow_qc),apesf
c            write(*,*)'yesf',yesf(kcol-1),yesf(kcol)
            rsc(irow_qc)=scale(irow_qc-1)*      ! error scaling modified
     &      (1/apesf+yesf(kcol-1)/yesf(kcol)**2)/
     &      (1.0/apesf**2+1/yesf(kcol)**2)
c            write(*,*)'scale=',scale(irow_qc),scale(irow_qc-1),
c    &      ' rsc=',rsc(irow_qc)
            jcol_gaa = jcol_gaa-2
         end do
      endif

c  Look within each data record to see if any of the data values are
c  outside their VMIN to VMAX range. If so, set eflag to the index of
c  the variable that was furthest out of range. Then write out the data.
      nco=0
      eflag=0
      kmax=0
      dmax=0.0
      do icol=1+mchar,ncol
         krow_qc=pindex(icol)
         if (yrow(icol).ge.9.000E+29) then
            write(*,*)'Missing value found (>=9E29). Ending program.'
            write(*,*)'You may need to remove missing .col files '//
     & 'from multiggg.sh'
            write(*,*)'and rerun post_processing.sh.'
            stop
         else
            dev=abs((rsc(krow_qc)*yrow(icol)-vmin(krow_qc))/
     &      (vmax(krow_qc)-vmin(krow_qc))-0.5) 
         endif
         if(dev.gt.dmax) then
            dmax=dev
            kmax=krow_qc
         endif
         if(flag(krow_qc).ge.1) then 
c This if statement eliminates columns in the output file that aren't 
c included in the qc.dat file. To include everything by default, 
c remove the if statement
            nco=nco+1
            yrow(nco)=yrow(icol)*rsc(krow_qc)
         endif
      end do  ! do icol=1+mchar,ncol
      if(dmax.gt.0.5) then
         eflag=kmax
         nflag=nflag+1
         kflag(kmax)=kflag(kmax)+1
      endif
      yrow(1)=nint(yrow(1)-yrow(2)/365.25)   ! Year
      yrow(2)=nint(yrow(2)-yrow(3)/24.0)     ! Day

      oof_flag(irow)=eflag
      write(oof_out,outfmt1)eflag,specname,
     &(wlimit(dble(yrow(j)),ofmt(j)),j=1,nco)
      oof_out=oof_out(:lnbc(oof_out))
         

      end

