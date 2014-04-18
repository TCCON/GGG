      subroutine write_col(lun_col,colfile_format,specname,lrmax,nit,
     & ncbf,xzo,xfs,
     & cont_level,cont_tilt,cont_curv,
     & rmsocl,sgshift,gint,zmin,oloscol,overcol,cx,ex,ntg,nfp)
c  Writes one line to the already-opened .col file.

      implicit none
      integer*4 lun_col,nit,i,
     & ncbf,
     & jtg,ktg,ntg,nfp,lrmax
      character specname*(*),colfile_format*(*)
      real*4 rmsocl,sgshift,zmin,xzo,xfs,
c     & xcx(0:nfp+1),
     & cont_level,cont_tilt,cont_curv,
     & oloscol(ntg),overcol(ntg),cx(nfp),ex(nfp)
      real*8 gint,wlimit

      if(lun_col.eq.6) then
         ktg=min0(1,ntg) ! write only first target gas
      else
         ktg=ntg         ! write all target gases
      endif

       write(lun_col,colfile_format) specname(:lrmax),nit,
     & wlimit(dble(cont_level/1.0),'f5.3'),       ! CL
     & wlimit(-200*1.732*dble(cont_tilt),'f5.1'), ! CT
     & wlimit(2000*dble(cont_curv),'f4.1'),       ! CC
     & wlimit(1000*gint*dble(xfs),'f4.1'),        ! FS
     & wlimit(dble(sgshift),'f5.2'),
     & wlimit(dble(xzo),'f5.3'),
     & wlimit(dble(100*abs(rmsocl)),'f6.4'), 
     & wlimit(dble(zmin),'f8.3'),
     &(wlimit(dble(oloscol(jtg)/overcol(jtg)),'f7.3'),
     & overcol(jtg),
     & wlimit(dble(cx(jtg)),'f9.4'),ex(jtg),jtg=1,ktg)

c 76    format(1x,a,1x,i2,1x,f5.3,3(1x,f4.1),1x,f5.3,
c     & 1x,f6.4,f8.3,15(0pf7.3,1pe11.4,0pf9.4,1pe8.1))
       return
       end
