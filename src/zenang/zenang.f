c  Program computes zenith angle offsets from the CO2 VSF factors.
      include "../gfit/ggg_int_params.f"
      include "../comn/postproc_params.f"

      integer*4 k,ngas,kcol,jgas,luns,lunq,
     $ mobs,iobs,naux,ncol,ndatcol,nrow,lnbc,lr,ncomm,
     $ kco2, kpobs, ksza, idum
      parameter (luns=13)
      parameter (lunq=14)
      parameter (mobs=1500)
      character pabel*800,tavfile*12,clabel(mcolvav)*32,inputfmt*40,
     $ addn_lines(maddln)*(mcharhead),cdum*1
      real*8 yval(mcolvav),vsf,vsf_err,pobs,asza,del,dsbydt,fdum

      idum=mcolvsw   ! Avoid compiler warning (unused parameter)                                                       
      idum=mauxcol   ! Avoid compiler warning (unused parameter)                                                       
      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mlev      ! Avoid compiler warning (unused parameter) 
      idum=mrow_qc   ! Avoid compiler warning (unused parameter)
      idum=mspeci    ! Avoid compiler warning (unused parameter)
      idum=mvmode    ! Avoid compiler warning (unused parameter)
      idum=mgas      ! Avoid compiler warning (unused parameter)
      idum=nchar     ! Avoid compiler warning (unused parameter)
      idum=ncell     ! Avoid compiler warning (unused parameter)

      cdum(1:1) = countfmt(1:1) ! Avoid compiler warning (unused param)
 
      write(6,*)'ZENANG    Version 2.16     2019-04-17    GCT'
      if (iargc() == 0) then
         write(6,101)
 101     format('Enter name of .tav file (e.g.  fts93rat.tav, etc) ',$)
         read(*,85) tavfile 
      elseif (iargc() == 1) then
         call getarg(1, tavfile)
      else
         stop 'Usage: $gggpath/bin/zenang tavfilename'
      endif
 85   format(a)
      lr=lnbc(tavfile)
      if(lr.le.0) stop 'Empty tav file name'
c=====================================================================
c  Read in VSF values and errors from .TAV file
      open(luns,file=tavfile(:lr-3)//'tav',status='old')
      call read_postproc_header(luns, ncomm, ndatcol, nrow, naux,
     & fdum, inputfmt, addn_lines, naddn)
      read(luns,'(a)') pabel
      call substr(pabel,clabel,mcolvav,ncol)
      ngas=(ncol-naux)/2
      kco2 = -1
      kpobs = -1
      ksza = -1
      do kcol=1,ncol
        if(clabel(kcol) .eq. 'co2') kco2 = kcol
        if(clabel(kcol) .eq. 'pout') kpobs = kcol
        if(clabel(kcol) .eq. 'asza') ksza = kcol
        if(clabel(kcol) .eq. 'solzen') ksza = kcol
      end do
      
      if(kco2 .lt. 0) stop 'Did not find CO2'
      if(kpobs .lt. 0) stop 'Did not find pout'
      if(ksza .lt. 0) stop 'Did not find asza'

c     Convert the CO2 index to number of gases past the
c     aux columns, remembering there's a vsf and vsf_error
c     column per gas
      kco2 = (kco2 - naux - 1)/2 + 1
c
      open(lunq,file='zenang.out',status='unknown')
      do iobs=1,mobs
        read(luns,inputfmt,end=14)
     &  (yval(k),k=1,naux),(vsf,vsf_err,jgas=1,kco2)
        pobs=yval(kpobs)
        asza=yval(ksza) 
        del=(vsf-1)/dsbydt(asza)/(1+vsf_err**2)
        write(lunq,'(5f14.8)') yval(1),vsf,pobs,asza,del
      end do  ! iobs=1,mobs
      iobs=iobs+1
      read(luns,*,end=14)
      write(6,*) 'Warning: increase parameter MOBS'
14    close(luns)
      close(lunq)
      stop
      end
