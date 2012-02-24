c  Program to test averaging_with_additive_bias subroutine
      implicit none
      integer*4 lunr,mrow,nrow,irow,mcol,ncol,jcol,nval
      parameter (lunr=15,mcol=10,mrow=12,nval=200)
      real*4 scal,ymiss,yobs(nval),yerr(nval),
     & ybar(mrow),eybar(mrow),bias(mcol),ebias(mcol),
     & rew(mrow),cew(mcol),tew

      ymiss=-999.0
      scal=1.E-35
      open(lunr,file='test.dat',status='unknown')
      read(lunr,*)nrow,ncol
      do irow=1,nrow
         read(lunr,*)(yobs(irow+nrow*(jcol-1)),jcol=1,ncol)
         do jcol=1,ncol
            yobs(irow+nrow*(jcol-1))=scal*yobs(irow+nrow*(jcol-1))
         end do
         write(*,*)(yobs(irow+nrow*(jcol-1)),jcol=1,ncol)
         do jcol=1,ncol
            yerr(irow+nrow*(jcol-1))=0.1*scal+
     &      0.02*yobs(irow+nrow*(jcol-1))
         end do
      end do
      close(lunr)
 
      call average_with_add_bias
     & (ymiss,nrow,ncol,yobs,yerr,ybar,eybar,bias,ebias,rew,cew,tew)
      do irow=1,nrow
       write(*,'(9e12.4)')(ybar(irow)+bias(jcol),jcol=1,ncol),
     &  ybar(irow),eybar(irow),rew(irow)
      end do
      write(*,*)'-----------------------------------------------|'
      write(*,'(9e12.4)')(bias(jcol),jcol=1,ncol)
      write(*,'(9e12.4)')(ebias(jcol),jcol=1,ncol)
      write(*,'(9e12.4)')(cew(jcol),jcol=1,ncol)
      write(*,*)'SQRT(Chi2/N)=',tew
      write(*,*)

      do irow=1,nrow
         write(*,*)(yobs(irow+nrow*(jcol-1)),jcol=1,ncol)
      end do
      write(*,*)
      call average_with_mul_bias
     & (ymiss,nrow,ncol,yobs,yerr,ybar,eybar,bias,ebias,rew,cew,tew)
      do irow=1,nrow
       write(*,'(9e12.4)')(ybar(irow)*bias(jcol),jcol=1,ncol),
     & ybar(irow),eybar(irow),rew(irow)
      end do
      write(*,*)'-----------------------------------------------|'
      write(*,'(9e12.4)')(bias(jcol),jcol=1,ncol)
      write(*,'(9e12.4)')(ebias(jcol),jcol=1,ncol)
      write(*,'(9e12.4)')(cew(jcol),jcol=1,ncol)
      write(*,*)'SQRT(Chi2/N)=',tew

      stop
      end

      include 'average_with_add_bias.f'
      include 'average_with_mul_bias.f'
      include '../comn/vsubs.f'
