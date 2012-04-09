c  subroutine read_oneline_oof.f
c

      subroutine read_oneline_oof(inputlun, oof_out,
     & outfmt1, oof_flag, irow,
     & vmin, vmax, pindex,
     & nflag, eflag, kflag, flag,
     & ncol, mchar, scale, ofmt,gaa_naux,yrow1
     & )

      

      implicit none
      include "../ggg_int_params.f"
      include "params.f"

      integer*4 ncol, irow, gaa_naux
      integer*4 nflag, eflag, inputlun, mchar, lnbc
      integer*4 oof_flag(mrow)
      real*4 vmin(mrow_qc), vmax(mrow_qc)
      integer*4 kflag(mrow_qc), flag(mrow_qc)
      integer*4 pindex(mcol)
      character*512 outfmt1
      character*800 oof_out
      character ofmt(mrow_qc)*4
      real*4 scale(mrow_qc)

      !local
      character*(nchar) specname
      integer*4 kmax, krow_qc, nco
      integer*4 j,jj,icol
      real*4 yrow(mcol),yrow1(mcol)
      real*8 wlimit
      real*4 dmax, dev
      
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
            dev=abs((scale(krow_qc)*yrow(icol)-vmin(krow_qc))/
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
               yrow(nco)=yrow(icol)*scale(krow_qc)
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
     &   (wlimit(dble(yrow(j)),ofmt(j)),j=1,nco)
         oof_out=oof_out(:lnbc(oof_out))
         

      end

