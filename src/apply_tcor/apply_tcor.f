c  apply_t_corr.f
c
c  Applies a T-correction to the contents of a
c  user-prescribed xxxxxx.xsw file created by collate_results.
c
c  Input files:
c    xxxxxx.xsw   ! Raw column amounts
c    tsen.coef3   ! Coefficients describing T-sensitivity of each window
c
c  Output file:
c    xxxxxx.xsw.tcor  ! T-Corrected column amounts
c
c  Program first reads the file containing the tcor coefs
c  for each window, derived by compute_tsen_coef.f. Program
c  then reads the xxxxxxx.vsw file, one line/spectrum at a
c  time, computing the airmass and H2O column.  It then
c  computes Tdel from the h2o windows by regressing their
c  column amounts windows against their T-sensitivities.
c  It then uses Tdel to T-correct the column amounts for
c  *all* windows (not just h2o) and writes them out to
c  xxxxxx.xsw.tcor.
c
c  Since spectra are processed independently, one by one,
c  there is no limit on the number of rows/spectra that
c  can be handled from the .xsw file.
c
c  Assumes that the gas column measured from window i,
c  Y(i), is related to the true column, Yt, by
c      Y(i) = Yt.(1+Tdel.X(i))                           (1)
c  where Tdel is the unknown temperature error, and X(i)
c  is the fractional T-sensitivity (1/Y.dY/dT) of the
c  i'th window.
c     X(i) = a(i)+b(i).airmass+c(i).h2o_bar+d(i)*...
c  The a,b,c,d coefficients are read from tsen.coefN,
c  having previously been determined off-line by the
c  program compute_tsen_coef.f
c
c  Both h2o and th2o windows are used to determine Tdel.
c  The only reason for labeling the high-E" windows as
c  "th2o" is to prevent them being used here for the
c  h2o_bar calculation or subsequently (by average_results)
c  for the xh2o calculation.
c
c  Having determined Tdel from the th2o and h2o windows, the
c  column amounts of all windows/gases can be corrected using
c     Yt(i) = Y(i)/(1+Tdel.X(i))
c
c  These corrected values are then written into xxxxxx.vsw.tcor
c  which is identical in format and content to xxxxxx.vsw
c
c  Integer variables beginning with:
c    k  denotes windows in the tsen_coef file
c    j  denotes columns/windows of the .vsw file.
c    l  denotes rows/spectra of the .vsw file.
c    i  denotes the elements of Yobs, and Yerr.
c          i = (j+1-nauxcol)/2
c          j = 2i-1+nauxcol
c     xsw_window(j) = tsen_window(kwindex(j))

      implicit none
      include "../comn/postproc_params.f"

      integer lspec,i,j,k,mlabel,lr,lrt,
     & lunr_tsen,mcoef,mauxcol,
     & lunr_xsw,lunw_tcor,lnbc,lnblnk,
     & jpasza,jptout,jh2o1,jh2o2,jth2o1,jth2o2,
     & i1,i2,kwin,jf,jr,
     & mwin,nwin_xsw,nwin_tsen,jwin,
     & nauxcol,nlhead,ncol_xsw,ncol_tsen,nfc,nspec_xsw,nss,
     & idum

      parameter (mcoef=5)          ! Max number of coefficients
      parameter (mauxcol=40)       ! Max number of auxiliary columns
      parameter (lunr_tsen=50)     ! Input file (.xsw)
      parameter (lunr_xsw=51)      ! Input file (.xsw)
      parameter (lunw_tcor=54)     ! Output file (.xav)
      parameter (mwin=600)         ! Total number of columns/windows
      parameter (mlabel=18000)     ! Max Number of column lable characters

      integer  spectrum_flag,kwindex(2*mwin+mauxcol)

      character
     & dl*1,                !forward or backward slash
     & gggdir*256,      !ggg directory path (GGGPATH?)
     & version*64,
     & gfit_version*64,
     & gsetup_version*64,
     & collate_version*64,
     & swfile*80,
     & outfile*80,
     & collabel*(mlabel),
     & col1*1,
     & cnfc*1,
     & tsen_window(mwin)*10,
     & ftype*1,
     & coeffile*128,
     & xsw_window(2*mwin+mauxcol)*24,
     & spectrum*57, 
     & textstr*40,
     & data_fmt*80

      real*8 year,dtiny,tnum,tden

      real*4
     & zz(mcoef),
     & tsen(mwin),ts,
     & coef(mcoef,mwin),
     & d2r,
     & tdel,terr,
     & h2o_bar,
     & airmass,
     & yaux(mauxcol),
     & yobs(mwin),yerr(mwin) 

      idum = maddln     ! Avoid compiler warning (unused parameter)
      idum = mcharhead  ! Avoid compiler warning (unused parameter)

      version=
     &' apply_t_corr          Version 1.28    2020-03-12   GCT/JLL'
      write(*,*) version

      spectrum_flag=0
      dtiny=1.D-88
      d2r=3.14159265/180.0

      if (iargc() == 0) then
         write(*,'(a)')
     &   'Enter name of .?sw file whose contents are to be corrected'
         read(*,'(a)') swfile
         write(*,*)'Number of coefficients to describe T-sensitivity'
         read(*,*) nfc
      elseif (iargc() == 2) then
         call getarg(1, swfile)
         call getarg(2, cnfc)
         read(cnfc,*) nfc
      else
         stop 'Usage: $gggpath/bin/apply_t_corr vswfile N'
      endif
      lr=lnblnk(swfile)
      outfile=swfile(:lr)//'.tcor'
      ftype=swfile(lr-2:lr-2)

c     Platform specification:      DG090519
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)     !Length of gggdir

c  Read entire file containing the tcorr coefficients
c  for all windows.
      write(coeffile,'(a,i1)') gggdir(:lrt)//'tsen_coefs/tsen.coef',nfc
      open(lunr_tsen,file=coeffile,status='old')
      read(lunr_tsen,*) nlhead,ncol_tsen
      nfc=ncol_tsen-1
      do k=2,nlhead
         read(lunr_tsen,'(a)') collabel
      end do
      do kwin=1,mwin
         read(lunr_tsen,*,end=77) tsen_window(kwin),
     &   (coef(j,kwin),j=1,nfc)
      end do
      read(lunr_tsen,*,end=77) 
      stop 'nwin_tsen > mwin'
77    nwin_tsen=kwin-1
      close(lunr_tsen)
      write(*,*)'nwin_tsen =',nwin_tsen

c  Read the xxxxxx.xsw file header
      open(lunr_xsw,file=swfile,status='old')
      open (lunw_tcor,file=outfile,status='unknown')
      read(lunr_xsw,*) nlhead,ncol_xsw,nspec_xsw,nauxcol
      write(lunw_tcor,countfmt)
     & nlhead+1,ncol_xsw,nspec_xsw,nauxcol
      write(lunw_tcor,'(a)') version(:lnblnk(version))
      read(lunr_xsw,'(a)') collate_version
      write(lunw_tcor,'(a)') collate_version(:lnblnk(collate_version))
      read(lunr_xsw,'(a)') gfit_version
      write(lunw_tcor,'(a)') gfit_version(:lnblnk(gfit_version))
      read(lunr_xsw,'(a)') gsetup_version
      write(lunw_tcor,'(a)') gsetup_version(:lnblnk(gsetup_version))
      read(lunr_xsw,'(a)') textstr
      write(lunw_tcor,'(a)')textstr(:lnblnk(textstr))
      read(lunr_xsw,'(a)') data_fmt
      write(lunw_tcor,'(a)') data_fmt(:lnblnk(data_fmt))
      read(lunr_xsw,'(a)') collabel
      write(lunw_tcor,'(a)') collabel(:lnblnk(collabel))

      nwin_xsw=(ncol_xsw-nauxcol)/2
      call substr(collabel,xsw_window,2*mwin+mauxcol,nss)
      if(nss.ne.ncol_xsw) then
         write(*,*)' nss, ncol_xsw = ',nss,ncol_xsw
         stop 'nss .ne. ncol_xsw'
      endif

c  Index the windows in tsen.coef3 to those
c  in the header of xxxxx.xsw file such that
c     xsw_window(j) = tsen_window(kwindex(j))
c  This is necessary in case of differences
c  in the order of the contents of the two files.
      call clistindex(nss,xsw_window,nwin_tsen,tsen_window,kwindex)

c  Find starting and ending locations (i.e. column numbers)
c  of the ASZA, h2o, and th2o windows in the XSW file. These
c  will be used to compute the Airmass and the mean H2O column.
c  Ending indices are searched forwardly.
c  Starting indices are searched in reverse.
      jpasza=0  ! Index of ASZA
      jptout=0  ! Index of Tout
      jh2o1=0   ! starting index of H2O windows
      jh2o2=0   ! ending index of H2O windows
      jth2o1=0  ! starting index of tH2O windows
      jth2o2=0  ! ending index of tH2O windows
      do jf=1,nss     ! forward index
         jr=nss-jf+1   ! reverse index
         if( xsw_window(jf)(1:4) .eq. 'tout' ) jptout=jf
         if( xsw_window(jf)(1:4) .eq. 'asza' ) jpasza=jf
         if( xsw_window(jf)(1:5) .eq. 'th2o_') jth2o2=jf
         if( xsw_window(jr)(1:5) .eq. 'th2o_') jth2o1=jr
         if( xsw_window(jf)(1:4) .eq. 'h2o_' ) jh2o2=jf
         if( xsw_window(jr)(1:4) .eq. 'h2o_' ) jh2o1=jr
c         if(kwindex(jf).gt.0) write(*,*)jf,kwindex(jf),xsw_window(jf),
c     &     tsen_window(kwindex(jf))
      end do
      write(*,*)' asza, h2o1 h2o2 =',jpasza,jh2o1,jh2o2,jth2o1,jth2o2
c
c  Check that h2o and th2o windows are consecutive (with h2o first)
      if(jh2o2+1.ne.jth2o1) stop ' jh2o2+1 .ne. jth2o1 '

      if (index(collabel, 'spectrum') .gt. 0) spectrum_flag=1

c  Main loop reading records one-by-one from the XSW file
      do lspec=1,nspec_xsw
         if (spectrum_flag .eq. 1) then
            read(lunr_xsw,data_fmt(8:)) spectrum,col1,year,
     &      (yaux(i),i=3,nauxcol),(yobs(i),yerr(i),i=1,nwin_xsw)
         else
            read(lunr_xsw,data_fmt(8:)) col1,year,
     &      (yaux(i),i=2,nauxcol),(yobs(i),yerr(i),i=1,nwin_xsw)
         endif

c  Compute h2o_bar (averaging over h2o windows, not th2o)
         tnum=0.0d0
         tden=0.0d0
         i1=(jh2o1+1-nauxcol)/2
         i2=(jh2o2-nauxcol)/2
         do i=i1,i2
            tnum=tnum+(yobs(i)/yerr(i))/yerr(i)
            tden=tden+(1.d0/yerr(i))/yerr(i)
         end do
         h2o_bar=tnum/tden

c  Compute "Co-ordinates" of observation.
         airmass=1.01/(0.02+cos(d2r*yaux(jpasza)))  ! Airmass
         zz(1)=1.0                                ! Constant
         zz(2)=sqrt(airmass)                      
         zz(3)=1.0E-22*h2o_bar           ! H2O vertical column (scaled)
         zz(4)=yaux(jptout)              ! Surface temperature (C)
         zz(5)=1.0E-22*h2o_bar*airmass   ! H2O slant column (scaled) 

c  Compute T-sensitivities of th2o & h2o windows and save.
         i2=(jth2o2-nauxcol)/2
         do i=i1,i2
            jwin=2*i-1+nauxcol
            kwin=kwindex(jwin)
            if(kwin.le.0) then
               write(*,*)'Warning: window not found in tsen.coefs3: '//
     &         xsw_window(jwin)
               tsen(i)=0.0
            else
               call vdot(zz,1,coef(1,kwin),1,tsen(i),nfc)
c               tsen(i)=coef(1,kwin)+zz(2)*coef(2,kwin)
c     &         +zz(3)*coef(3,kwin)
            endif
         end do

c  Compute Tdel
         call compute_tdel(i2-i1+1,yobs(i1),yerr(i1),tsen(i1),tdel,terr)
       write(*,'(i4,1x,a,f9.4,1pe9.0,0pf7.0,2f9.4)')lspec,spectrum(:24),
     &   airmass,h2o_bar,273.15+yaux(jptout),tdel,terr

c  Apply Tdel correction to all gases (including H2O).
         do i=1,nwin_xsw
            jwin=2*i-1+nauxcol
            kwin=kwindex(jwin)
            if(kwin.le.0) then
               write(*,*)'Warning: window not found in tsen.coefs?: '//
     &         xsw_window(jwin)
               ts=0.0
            else
               call vdot(zz,1,coef(1,kwin),1,ts,nfc)
c               ts=coef(1,kwin)+zz(2)*coef(2,kwin)
c    &         +zz(3)*coef(3,kwin)
            endif
            yobs(i)=yobs(i)/(1.0+tdel*ts)
            yerr(i)=yerr(i)/(1.0+tdel*ts)
         end do
         
c  Write out corrected values to .xsw.tcor file
         if (spectrum_flag .eq. 1) then   ! Spectrum names are included
            write(lunw_tcor,data_fmt(8:)) spectrum,col1,year,
     &      (yaux(i),i=3,nauxcol),(yobs(i),yerr(i),i=1,nwin_xsw)
         else
            write(lunw_tcor,data_fmt(8:)) col1,year,
     &      (yaux(i),i=2,nauxcol),(yobs(i),yerr(i),i=1,nwin_xsw)
         endif

      end do  !  lspec=1,nspec_xsw
      close(lunr_xsw)
      close(lunw_tcor)
      stop
      end
