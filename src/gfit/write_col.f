      subroutine write_col(lun_col,runlab,lrmax,nit,rms,sgshift,gint,
     & zmin,oloscol,overcol,cx,ex,ntg)
c  Writes one line to the already-opened .col file.

      implicit none
      integer*4 lun_col,nit,ntg,n1,n2,n3,n4,jtg,ktg,lrmax
      character runlab*(*)
      real*4 rms,sgshift,zmin,
     & oloscol(ntg),overcol(ntg),cx(ntg+4),ex(ntg+4)
      real*8 gint,wlimit

      n1=ntg+1
      n2=ntg+2
      n3=ntg+3
      n4=ntg+4
      if(lun_col.eq.6) then
         ktg=min0(1,ntg) ! only first target gas
      else
         ktg=ntg         ! all target gases
      endif

       write(lun_col,76)runlab(:lrmax),nit,
     & wlimit(dble(cx(n1)),'f5.3'),
     & wlimit(100*dble(cx(n2)),'f4.1'),
     & wlimit(1000*gint*dble(cx(n3)),'f4.1'),
     & wlimit(dble(sgshift),'f4.1'),
     & wlimit(dble(cx(n4)),'f5.3'),
     & wlimit(dble(100*abs(rms/cx(n1))),'f6.4'), 
     & wlimit(dble(zmin),'f8.3'),
     &(wlimit(dble(oloscol(jtg)/overcol(jtg)),'f7.3'),
     & overcol(jtg),
     & wlimit(dble(cx(jtg)),'f9.4'),ex(jtg),jtg=1,ktg)

 76    format(1x,a,1x,i2,1x,f5.3,3(1x,f4.1),1x,f5.3,
     & 1x,f6.4,f8.3,15(0pf7.3,1pe11.4,0pf9.4,1pe8.1))
       return
       end
