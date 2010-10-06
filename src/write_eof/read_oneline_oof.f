c  Program: create_official_output_file.f
c
c  Purpose: To convert the runlog.vav.ada.aia file to
c  an official format.
c  
c  Input Files:
c       runlog.vav.ada.aia 
c       qc_limits.dat  
c
c  Output Files:
c       runlog.vav.ada.aia.oof
c       

      subroutine read_oneline_oof(inputfile,inputlun, oof_out,
     & outfmt1, oof_flag, irow,
     & vmin, vmax, pindex,
     & nflag, eflag, kflag, flag,
     & nrow, ncol, nchar, scale, ofmt
     & )

      

      implicit none
      include "params.f"
      integer*4 ncol, nrow, irow
      integer*4 nflag, eflag, inputlun, nchar, lnbc
      integer*4 oof_flag(mrow)
      real*4 vmin(mrow_qc), vmax(mrow_qc)
      integer*4 kflag(mrow_qc), flag(mrow_qc)
      integer*4 pindex(mcol)
      character*512 outfmt1
      character*800 oof_out
      character*80 inputfile
      character ofmt(mrow_qc)*4
      real*4 scale(mrow_qc)

      !local
      character*38 specname,cc
      character*800 ssss
      integer*4 kmax, krow_qc, nco
      integer*4 j,k, kcol, icol
      real*4 yrow(mcol)
      character sarr(mcol)*38
      real*8 wlimit
      real*4 dmax, dev
      

         if (nchar .eq. 1) then
             read(inputlun,*) specname, (yrow(j),j=1+nchar,ncol)
         else
             read(inputlun,*) (yrow(j),j=1,ncol)
         endif

c  Look within each data record to see if any of the data values are
c  outside their VMIN to VMAX range. If so, set eflag to the index of
c  the variable that was furthest out of range. Then write out the data.
         nco=0
         eflag=0
         kmax=0
         dmax=0.0
         do icol=1+nchar,ncol
            krow_qc=pindex(icol)
            dev=abs((scale(krow_qc)*yrow(icol)-vmin(krow_qc))/
     &      (vmax(krow_qc)-vmin(krow_qc))-0.5)
            if(dev.gt.dmax) then
               dmax=dev
               kmax=krow_qc
            endif
            if(flag(krow_qc).ge.1) then
               nco=nco+1
               yrow(nco)=yrow(icol)*scale(krow_qc)
            endif
         end do  ! do icol=1,ncol
         if(dmax.gt.0.5) then
            eflag=kmax
            nflag=nflag+1
            kflag(kmax)=kflag(kmax)+1
         endif
         yrow(1)=nint(yrow(1)-yrow(2)/365.25)   ! Year
         yrow(2)=nint(yrow(2)-yrow(3)/24.0)     ! Day

         write(ssss,outfmt1)eflag,specname,
     &   (wlimit(dble(yrow(j)),ofmt(j)),j=1,nco)
         call substr(ssss,sarr,mcol,kcol)
         do k=2,kcol
            cc=sarr(k)
            ssss=ssss(:lnbc(ssss))//','//cc(:lnbc(cc))
         end do
         oof_flag(irow)=eflag
         write(oof_out,outfmt1)eflag,specname,
     &   (wlimit(dble(yrow(j)),ofmt(j)),j=1,nco)
         oof_out=oof_out(:lnbc(oof_out))
         

      end

