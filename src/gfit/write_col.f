      subroutine write_col(lun_col,colfile_format,specname,lrmax,nit,
     & rms,sgshift,gint,zmin,oloscol,overcol,cx,ex,ntg,nfp)
c  Writes one line to the already-opened .col file.

      implicit none
      integer*4 lun_col,nit,
     & n1,n2,n3,n4,n5,
     & jtg,ktg,ntg,nfp,lrmax
      character specname*(*),colfile_format*(*)
      real*4 rms,sgshift,zmin,
     & oloscol(ntg),overcol(ntg),cx(nfp),ex(nfp)
      real*8 gint,wlimit

      n1=ntg+1
      n2=ntg+2
      n3=ntg+3
      n4=ntg+4
      n5=ntg+5
      if(lun_col.eq.6) then
         ktg=min0(1,ntg) ! write only first target gas
      else
         ktg=ntg         ! write all target gases
      endif

       write(lun_col,colfile_format) specname(:lrmax),nit,
     & wlimit(dble(cx(n1)),'f5.3'),       ! CL
     & wlimit(100*dble(cx(n2)),'f5.1'),   ! CT
     & wlimit(1000*dble(cx(n3)),'f4.1'),  ! CC
     & wlimit(1000*gint*dble(cx(n4)),'f4.1'),
     & wlimit(dble(sgshift),'f4.1'),
     & wlimit(dble(cx(n5)),'f5.3'),
     & wlimit(dble(100*abs(rms/cx(n1))),'f6.4'), 
     & wlimit(dble(zmin),'f8.3'),
     &(wlimit(dble(oloscol(jtg)/overcol(jtg)),'f7.3'),
     & overcol(jtg),
     & wlimit(dble(cx(jtg)),'f9.4'),ex(jtg),jtg=1,ktg)

c 76    format(1x,a,1x,i2,1x,f5.3,3(1x,f4.1),1x,f5.3,
c     & 1x,f6.4,f8.3,15(0pf7.3,1pe11.4,0pf9.4,1pe8.1))
       return
       end
