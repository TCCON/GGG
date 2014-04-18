c  Program linefinder.f
c
c  Reads a spectrum and uses the 3'rd derivative spectrum
c  to identify the absorption lines.
c
c  Can also read a residual spectrum to identify missing
c  absorption lines.
c
c  Works for any spectrum having a GGG-format runlog,
c  which is used as the input file.
c
c  Byte-reversal handled automatically, if the computer
c  that you are working on has a different endian-ness
c  from the computer that wrote the binary spectra.
c
      implicit none
      include "../ggg_int_params.f"

      integer*4
     & lun_rlg,   ! LUN to read input runlogs from
     & luns,   ! LUN to read binary spectras from
     & lunw,   ! LUN to write ascii spectra to
     & i,j,k,
     & iabpw,  ! absolute values of the bytes per word
     & lnbc,   ! function Last Non-Black Character
     & iend,   ! Endianess of host computer
     & mpts, npts,   ! Number of spectral values
c     & mterm, nterm, ! Number of polynomial terms to fit continuum
     & mhw, nhw, ! Half-width of smoothing operator
     & nlines,
     & lr   ! 

      parameter (lun_rlg=25,luns=15,lunw=16)
      parameter (mpts=1280*1024, mhw=20)
      real*4 bufr4(mpts),deriv1(mpts),deriv3(mpts),sd1(mpts),sd3(mpts),
     & thresh,frx,yyx,smoo_oper(mhw),tm,tc
c     & x(mpts),y(mpts),wt(mpts),jac(mpts,mterm),wk(mterm),c(mterm)

      character 
     & inpath*80,chead*1

      integer*4  nlhead, ncol, np,
     & kgas, kiso,   ! Gas # and Isotope code (first two columns of linelist)
     & kgas1, kiso1,   ! Gas # and Isotope code (first two columns of linelist)
     & kgas2, kiso2,   ! Gas # and Isotope code (first two columns of linelist)
     & kgas3, kiso3,   ! Gas # and Isotope code (first two columns of linelist)
     & istat,        ! status flag (0=success, 1=EOF)
     & iyr,          ! year
     & iset,         ! day of year
     & ifirst,       ! index of first spectral point in disk file
     & ilast,        ! index of last spectral point in disk file
     & m1, m2,       ! starting and ending spectral point indices in output file
     & possp,        ! Length of attached header in bytes
     & bytepw        ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      real*8 rms1d,rms3d,width,tot,tranmax,dnu,ycont,sfac,fwid,
     &   freq,stren,abhw,sbhw,eprime,tdpbhw,pshift,
     &   freq1,stren1,r1,abhw1,sbhw1,eprime1,tdpbhw1,pshift1,
     &   freq2,stren2,r2,abhw2,sbhw2,eprime2,tdpbhw2,pshift2,
     &   freq3,stren3,r3,abhw3,sbhw3,eprime3,tdpbhw3,pshift3,
c     & slope, rms_slope, sbar,
c     & curv, rms_curv,
     & oblat,        ! observation latitude (deg).
     & oblon,        ! observation longitude (deg).
     & obalt,        ! observation altitude (km)
     & asza,         ! astronomical solar zenith angle (unrefracted)
     & azim,         ! azimuth angle
     & opd,          ! Optical path difference (cm) of interferogram
     & graw,         ! spacing of raw spectrum (cm-1) from GETINFO
     & zpdtim,       ! Time of ZPD (UT hours)
     & zenoff,       ! Zenith angle pointing offset (deg)
     & fovi,         ! Internal angular diameter of FOV (radians)
     & fovo,         ! External angular diameter of FOV (radians)
     & amal,         ! angular misalignment of interferometer (radians)
     & zoff,         ! Zero level offset (dimensionless fraction)
     & snr,          ! Signal-to-Noise Ratio (dimensionless)
     & tins,         ! Inside temperature
     & pins,         ! Inside pressure
     & hins,         ! Inside humidity
     & tout,         ! Outside temperature
     & pout,         ! Outside pressure
     & hout,         ! Outside humidity
     & osds,
     & wspd,
     & wdir,
     & sia,          ! Solar Intensity (Average)
     & fvsi,
     & aipl,         ! Airmass-Independent Path Length (km)
     & lasf,         ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr,       ! suntracker frequency (active tracking)
     & nu1, nu2,
     & nus, nue      ! selected frequency range of interest

      character
     & sss*33,sss1*33,sss2*33,sss3*33,nus_str*40,nue_str*40,
     & ttt*33,
     & data_fmt_read_rl*256,
     & col_labels_rl*320,
     & outfile*18,
     & llfmtw*59
      character
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,
     & sllformat*23,
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & rlgfile*120                !name of runlog file

       logical residual_flag

      residual_flag=.false.
c      kgas=65   !  CHF3
c      kiso=1    !  Real linelist (not pseudo)

c      kgas=19    !  OCS
c      kiso=6     !  Real linelist (not pseudo)
c      abhw=0.085
c      sbhw=0.160

c       kgas=56   !  CH3OH
c       kiso=1    !  PLL
c       abhw=0.10
c       sbhw=0.40
c       eprime=400.
c       tdpbhw=0.7 
c       pshift=-0.003

       kgas=70   !  C3H8
       kiso=1    !  Real linelist (not pseudo)
       abhw=0.12
       sbhw=0.17
       eprime= 400.
       tdpbhw=0.65
       pshift=-0.003

c      residual_flag=.true.

      tranmax=1.002  !  Jeremy Harrison's C2H6 spectra
      tranmax=0.997  ! Keeyoon Sung's  C3H8 spectra
      sfac=3.6E-23   ! Strength scale factor
c      tranmax=0.980  ! KP OCS spectra
      tranmax=0.99  ! solar
      tranmax=0.999 ! Jeremy Harrison's CH3OH spectra

      sfac=1.0  ! solar spectrum
      sfac=1.3e-20  ! CH3OH

      write(6,*)
     & ' Linefinder Program    Version 1.43   17-May-2013   GCT'
      call getendian(iend)  ! Find endian-ness of host computer

      llfmtw=
     & '(i2,i1,f12.6,1pe10.3,e10.3,0pf5.4,f5.4,f10.4,f4.2,f8.6,a33)'
      sss1=' generated by linefinder  V 1.41 '
      ttt=' generated by linefinder  extras '

      sllformat='(i3,f13.6,6f9.5,4x,a33)'
c  Interrogate environmental variable GGGPATH to find location
c  of root partition (e.g. "/home/toon/ggg/" ).
      call get_ggg_environment(gggdir, dl)
      lr=lnbc(gggdir)     ! length of root string (e.g. 14)
c
      if (iargc() == 0) then
         Write(*,*)'Enter Starting & Ending frequencies:'
         write(*,*)'Enter 0 99999 to retain original spectral limits'
         read(*,*) nus,nue

         write(*,*)'Enter path to input file/runlog:'
         read(*,'(a)') rlgfile
      elseif (iargc() == 3) then
         call getarg(3, rlgfile)
         call getarg(1, nus_str)
         read(nus_str,*)nus
         call getarg(2, nue_str)
         read(nue_str,*)nue
      else
         write(*,*)'Use: $gggpath/bin/linefinder start_freq end_freq '//
     & 'path/runlog'
         stop
      endif
      open(lun_rlg,file=rlgfile,status='old')
      call read_runlog_header(lun_rlg,data_fmt_read_rl,col_labels_rl)
c      read(lun_rlg,*)  nlhead, ncol
c      do j=2,nlhead
c      read(lun_rlg,*)          ! Skip header line of runlog
c      end do

c  Define Gaussian Smoothing Operator
      width=1.5  ! in spectral points
      nhw=nint(2.0*width)
      if(nhw.gt.mhw) then
        write(*,*) 'mhw, nhw =',mhw,nhw
        stop 'increase parameter MHW'
      endif
      tot=0.0
      do k=1,nhw
         smoo_oper(k)=exp(-(float(k)/width)**2)
         tot=tot+smoo_oper(k)
      end do
      tot=1.0+2*tot  ! Operator is double-sided with center values of 1.0

c  Normalize smoo_oper to unit area
      do k=1,nhw
         smoo_oper(k)=smoo_oper(k)/tot
      end do

c  Read first line of runlog only.
      call read_runlog_data_record(lun_rlg,data_fmt_read_rl,
     & col1,specname,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,azim,osds,
     & opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,
     & tins,pins,hins,tout,pout,hout,
     & sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
      close(lun_rlg)

      write(*,*)'tins,tout=',tins,tout
      graw=graw*(1.0d0+osds*1.0E-06)
      if (nus.lt.ifirst*graw) nus=ifirst*graw
      if (nue.gt.ilast*graw)  nue=ilast*graw
      m1=1+int(nus/graw)-nhw
      m2=int(nue/graw)+nhw

c  Check that buffer will be large enough
      iabpw=iabs(bytepw)
      npts=m2-m1+1
      write(*,*)'ifirst, ilast, graw =',ifirst, ilast, graw
      write(*,*)'m1, m2, npts=',m1, m2, npts
      if(npts.gt.mpts) stop 'Increase parameter MPTS'
      if(npts.lt.5) stop 'NPTS < 5'

      if(graw*(m2-nhw).lt.999.0) then
         write(outfile,'(a3,i3.3,a1,i3.3,a4)')
     &   'lf_',nint(graw*(m1+nhw)),'_',nint(graw*(m2-nhw)),'.101'
      elseif(graw*(m2-nhw).lt.9999.0) then
         write(outfile,'(a3,i4.4,a1,i4.4,a4)')
     &   'lf_',nint(graw*(m1+nhw)),'_',nint(graw*(m2-nhw)),'.101'
      elseif(graw*(m2-nhw).lt.99999.0) then
         write(outfile,'(a3,i5.5,a1,i5.5,a4)')
     &   'lf_',nint(graw*(m1+nhw)),'_',nint(graw*(m2-nhw)),'.101'
      else
         write(outfile,'(a3,i6.6,a1,i6.6,a4)')
     &   'lf_',nint(graw*(m1+nhw)),'_',nint(graw*(m2-nhw)),'.101'
      endif

c  Search for spectrum "runlab"
      call gindfile(gggdir(:lr)//'config'//dl//'data_part.lst',
     &     specname,inpath)
      if(lnbc(inpath).eq.0) then
          write(*,*) ':'//specname
          stop ' Cannot find input spectrum'
      endif

c  Copy spectral fit file from $GGGPATH/spt/zxxxxxxx to local directory
      if(residual_flag) then
          inpath='z'//specname
          iabpw=7
      endif
      write(*,*) inpath
 
      if(iabpw.eq.4) then
c  Open binary spectrum with recl = total length be read
         open(luns,file=inpath,access='direct',status='old',
     &    form='unformatted',recl=possp+iabpw*(m2-ifirst+1))

c  Read spectral header and data values all at once.
         read(luns,rec=1) (chead,j=1,possp+iabpw*(m1-ifirst)),
     &   (bufr4(j),j=1,npts)
         close(luns)

c  If necessary, byte-reverse data
         if(iend*bytepw.lt.0) call rbyte(bufr4,iabpw,npts)
      else     ! ascii spt file 
         open(luns,file=inpath,status='old')
         read(luns,*) nlhead, ncol
         read(luns,*) nu1,nu2,np
         dnu=(nu2-nu1)/(np-1)
         ifirst=nint(nu1/dnu)
         ilast=nint(nu2/dnu)
         do j=3,nlhead+m1-ifirst
            read(luns,*)
         end do
         do j=1,npts
            read(luns,*)freq,tm,tc
            bufr4(j)=1+tm-tc
         end do
         close(luns)
      endif
c
cc  Write out spectrum (ascii format)
c      write(56,*) ' 2  2'
c      write(56,*) ' f  t'
c      do j=1,npts
c          write (56,*) graw*(j-1+m1),bufr4(j)
c      end do

c  Smooth spectrum/residuals to kill ringing
      do j=2,npts-1
         bufr4(j)=0.5*bufr4(j)+0.25*(bufr4(j-1)+bufr4(j+1))
      end do
  
c  Compute 1'st and 3'rd derivative spectra
      deriv1(1)=0.0
      deriv1(2)=(bufr4(3)-bufr4(1))/2
      deriv3(1)=0.0
      deriv3(2)=0.0
      do i=3,npts-2
         deriv1(i)=(bufr4(i+1)-bufr4(i-1))/2
         deriv3(i)=bufr4(i-1)-bufr4(i+1)-(bufr4(i-2)-bufr4(i+2))/2
      end do
      deriv3(npts-1)=0.0
      deriv3(npts)=0.0
      deriv1(npts-1)=(bufr4(npts)-bufr4(npts-2))/2
      deriv1(npts)=0.0
      
c Smooth Derivative spectra and find RMS
      rms1d=0.0
      rms3d=0.0
      do i=1+nhw,npts-nhw
         sd1(i)=deriv1(i)/tot
         sd3(i)=deriv3(i)/tot
         do k=1,nhw
            sd1(i)=sd1(i)+smoo_oper(k)*(deriv1(i-k)+deriv1(i+k))
            sd3(i)=sd3(i)+smoo_oper(k)*(deriv3(i-k)+deriv3(i+k))
         end do
         rms1d=rms1d+sd1(i)**2
         rms3d=rms3d+sd3(i)**2
      end do
      rms1d=sqrt(rms1d/(npts-nhw))
      rms3d=sqrt(rms3d/(npts-nhw))
      write(*,*) 'RMS 1D, 3D =',rms1d,rms3d

c  Find the +ve zero crossings of 3'rd derivative
      open(33,file=outfile,status='unknown')
      open(35,file=outfile(:lnbc(outfile)-4)//'.108',status='unknown')
      open(34,file='lf.out',status='unknown')
      write(34,*) 2,10
      write(34,*) 'f frx s y sd1i sd1p rms1d sd3i sd3p rms3d'
      thresh=0.25*rms3d
      nlines=0
      ycont=bufr4(nhw)
c      ycont=1.0
      do i=1+nhw,npts-nhw
         if( bufr4(i).gt.ycont) then
            ycont=bufr4(i)
         else
            ycont=ycont*0.9999
         endif
         
         if(sd3(i)*sd3(i+1).lt.0.000*rms3d) then   ! zero crossing of 3'rd derivative
            write(67,*)'ZC1',sd3(i),sd3(i+1),graw*(i+0.5+m1-1)
            if(sd3(i+1)-sd3(i).lt.-0.00*rms3d) then  ! sufficiently negative-going zero crossing
            write(67,*)'ZC2',sd3(i),sd3(i+1),graw*(i+0.5+m1-1)
            if(sd1(i+1)-sd1(i).gt.-1.15*rms1d) then  ! increasing slope (+ve curvature)
            write(67,*)'ZC3',sd3(i),sd3(i+1),graw*(i+0.5+m1-1)
               frx=sd3(i)/(sd3(i)-sd3(i+1))
               yyx=(1-frx)*bufr4(i)+frx*bufr4(i+1)
               stren=-log(yyx/tranmax/ycont)
c               stren=1-yyx/tranmax/ycont
c               write(*,*) graw*(i+frx+m1-1),stren
               if(stren.gt.0.001) then
                   nlines=nlines+1
                   freq=graw*(i+frx+m1-1)
                   stren=stren*exp(1.4388*eprime*(1.0/tout-1.0/296))
                   write(34,'(f11.5,2f8.4,7(1pe11.3))') freq,frx,stren,
     &             yyx,sd1(i),sd1(i+1),rms1d,sd3(i),sd3(i+1),rms3d
                   write(33,llfmtw) kgas,kiso,1.00*freq,sfac*stren,
     &             0.0, abhw, sbhw, eprime, tdpbhw, pshift, sss1 
         write(35,sllformat) 56,freq,
     &   stren,0.97*stren,0.07*sqrt(0.03+stren),0.07*sqrt(0.03+stren),
     &   8.5E-06*freq,8.5E-06*freq,sss
               endif
            endif
            endif
         endif
c  Insert extra lines into broad absorption features which are devoid of peaks.
         if( (i+m1-1)*graw .gt. freq+5*width*graw) then  ! fill gaps with extra lines
           yyx=bufr4(i)
           stren=-log(yyx/tranmax/ycont)
c          stren=1-yyx/tranmax/ycont
           if(stren.gt.2.0) then  ! high values reduce number of extras
               nlines=nlines+1
               freq=graw*(i+m1-1)
               write(34,'(f11.5,3f8.4,6(1pe11.3))') freq,frx,stren,
     &         yyx,sd1(i),sd1(i+1),rms1d,sd3(i),sd3(i+1),rms3d
               write(33,llfmtw) kgas,kiso,1.00*freq,sfac*stren,
     &         0.0, abhw, sbhw, eprime, tdpbhw, pshift, ttt
         write(35,sllformat) 56,freq,
     &   stren,0.97*stren,0.07*sqrt(0.03+stren),0.07*sqrt(0.03+stren),
     &   8.5E-06*freq,8.5E-06*freq,ttt

           endif
         endif
      enddo
      close(33)
      close(34)
      write(*,*) nlines

c  Write ASCI derivative spectrum (for diagnostic purposes)
      write(6,*)inpath(:lnbc(inpath))
      open(lunw,file='./der_'//specname,status='unknown')
      write(lunw,*)2,4
      write(lunw,*)' f  S  SD1  SD3 '
      do i=1,npts
         write(lunw,'(f12.6,3(1pe12.4))') graw*(i+m1-1),
     &   bufr4(i)-1,10*sd1(i),250*sd3(i)
      end do
      close(lunw)

cc  Correct line strengths for overlap
c      fwid=4.0*width*graw
c      open(33,file=outfile,status='old')
c      open(34,file='sc'//outfile(3:),status='unknown')
c         read(33,llfmtw) kgas1,kiso1,
c     &   freq1,stren1,r1,abhw1,sbhw1,eprime1,tdpbhw1,pshift1,sss1
c         read(33,llfmtw) kgas2,kiso2,
c     &   freq2,stren2,r2,abhw2,sbhw2,eprime2,tdpbhw2,pshift2,sss2
c         stren1=stren1-0.5*stren2*exp(-((freq2-freq1)/fwid)**2)
c         write(34,llfmtw)kgas1,kiso1,
c     &   freq1,stren1,r1,abhw1,sbhw1,eprime1,tdpbhw1,pshift1,sss1
c         do i=3,nlines
c         read(33,llfmtw) kgas3,kiso3,
c     &   freq3,stren3,r3,abhw3,sbhw3,eprime3,tdpbhw3,pshift3,sss3
c         stren2=stren2-0.5*stren1*exp(-((freq2-freq1)/fwid)**2)
c     &   -0.5*stren3*exp(-((freq3-freq2)/fwid)**2)
c         write(34,llfmtw) kgas2,kiso2,
c     &   freq2,stren2,r2,abhw2,sbhw2,eprime2,tdpbhw2,pshift2,sss2
c         stren1=stren2
c         freq1=freq2
c         eprime1=eprime2
c         stren2=stren3
c         freq2=freq3
c         eprime2=eprime3
c         end do
c         stren2=stren2-0.5*stren1*exp(-((freq2-freq1)/fwid)**2)
c         write(34,llfmtw) kgas2,kiso2,
c     &   freq2,stren2,r2,abhw2,sbhw2,eprime2,tdpbhw2,pshift2,sss2
c      close(33)
c      close(34)

      stop
      end
