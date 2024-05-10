      subroutine write_col(lunw_col,colfile_format,ntogo,specname,lrmax,
     & nit,xzo,xfs,
     & cont_level,cont_tilt,cont_curv,
     & rmsocl,xsg,zmin,oloscol,overcol,cx,ex,ntg,nfp)
c  Writes one line to the already-opened .col file.

      implicit none
      integer*4 lunw_col,nit,ntogo,
     & jtg,ktg,ntg,nfp,lrmax
      character specname*(*),colfile_format*(*)
      real*4 rmsocl,xsg,zmin,xzo,xfs,
     & cont_level,cont_tilt,cont_curv,
     & oloscol(ntg),overcol(ntg),cx(nfp),ex(nfp)
      real*8 wlimit

      if(lunw_col.eq.6) then
         ktg=min0(1,ntg) ! write only first target gas to screen
         write(lunw_col,'(i4,'//colfile_format(2:)) ntogo,
     &   specname(:lrmax),nit,
     &   wlimit(dble(cont_level/1.0),'f5.3'),       ! CL
     &   wlimit(-200*1.732*dble(cont_tilt),'f5.1'), ! CT
     &   wlimit(2000*dble(cont_curv),'f4.1'),       ! CC
     &   wlimit(dble(1.E+06*xfs),'f5.2'),           ! FS
     &   wlimit(dble(1.E+06*xsg),'f5.2'),           ! S-G
     &   wlimit(dble(xzo),'f6.4'),
     &   wlimit(dble(100*abs(rmsocl)),'f6.4'), 
     &   wlimit(dble(zmin),'f8.3'),
     &  (wlimit(dble((oloscol(jtg)+.01)/(overcol(jtg)+.01)),'f7.3'),
     &   overcol(jtg),wlimit(dble(cx(jtg)),'f10.5'),ex(jtg),jtg=1,ktg)

      else
         ktg=ntg         ! write all target gases to .col file
         write(lunw_col,colfile_format) 
     &   specname(:lrmax),nit,
     &   wlimit(dble(cont_level/1.0),'f5.3'),       ! CL
     &   wlimit(-200*1.732*dble(cont_tilt),'f5.1'), ! CT
     &   wlimit(2000*dble(cont_curv),'f4.1'),       ! CC
     &   wlimit(dble(1.E+06*xfs),'f5.2'),           ! FS
     &   wlimit(dble(1.E+06*xsg),'f5.2'),           ! S-G
     &   wlimit(dble(xzo),'f6.4'),
     &   wlimit(dble(100*abs(rmsocl)),'f6.4'), 
     &   wlimit(dble(zmin),'f8.3'),
     &  (wlimit(dble((oloscol(jtg)+.01)/(overcol(jtg)+.01)),'f7.3'),
     &   overcol(jtg),wlimit(dble(cx(jtg)),'f10.5'),ex(jtg),jtg=1,ktg)
      endif

      return
      end
