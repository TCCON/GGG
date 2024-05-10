c  Program to compare two similar spectra and quantify the differences.
c  Useful in I2S development to see whether changes to a spectrum are
c  significant or simply rounding error, since the "diff" command simply
c  says whether files are identical or not.
c
c  Currently hard-wired for OPUS format binary spectra.
c  Currently hard-wired for MkIV format binary spectra.
c
c  Computes the L2 norm of spectrum 1, the rms difference between
c  spectrum 1 and 2, and their ratio.
c
      implicit none

      integer*4
     & ispec,npt,mpt,ls,lloc,lunr_bs,
     & iy,im,id,lst,j,lnblnk,
     & ncount,      !  Number of spectral values that were different
     & ifirst(2),   !  Index of the first spectral point on disk
     & ilast(2),    !  Index of the last spectral point on disk
     & possp,       !  Bytes before IFIRST'th point (i.e. header length)
     & bpw          !  Byte per data word (2 or 4)

      parameter (lunr_bs=14,mpt=1024*2048)

      real*4 r4buf(mpt,2),diff
      integer*2 i2buf(mpt,2)
      equivalence (i2buf,r4buf)

      real*8 delwav,
     & tb,td,
     & lse,lsu,lsf,dip,mvd,
     $ gmt,fovi,lasf,
     & opd,
     $ snr,wavtkr

      character
     & cdum*1,
     & spfmt*2,
     & apf*2,
     & dl*1,
     & specpath(2)*128,
     & version*48                 !current program version


      version=' spec_diff    Version 1.05     2019-04-01    GCT'
      write(*,'(a)') version
      dl='/'

      spfmt='op'  ! OPUS
c      spfmt='m4'  ! MkIV

      if (iargc() == 0) then
         write(6,'(a)') 'Enter full paths to spectra:'
         read(*,'(a)') specpath(1),specpath(2)
      elseif (iargc() == 2) then
         call getarg(1, specpath(1))
         call getarg(2, specpath(2))
      else
         stop 'Usage: $gggpath/bin/specdiff specpath1 specpath2'
      endif

      do ispec=1,2
         ls=lloc(specpath(ispec),dl)
c  Read spectrum header.
         call rdsphead(spfmt,specpath(ispec)(ls+1:),
     &   specpath(ispec)(:ls),ifirst(ispec),ilast(ispec),possp,bpw,apf,
     &   delwav,opd,fovi,snr,iy,im,id,gmt,lasf,wavtkr,lst,lse,lsu,lsf,
     &   dip,mvd)
         npt=ilast(1)-ifirst(1)+1
         if(npt.gt.mpt) stop 'npt > MPT'

c  Open binary spectrum with recl = total length be read
         open(lunr_bs,file=specpath(ispec),access='direct',status='old',
     &   form='unformatted',recl=possp+npt*iabs(bpw))

c  Read spectral header and data values all at once, then close.
         if(iabs(bpw).eq.4) then
            read(lunr_bs,rec=1)(cdum,j=1,possp),(r4buf(j,ispec),j=1,npt)
         elseif(iabs(bpw).eq.2) then
            read(lunr_bs,rec=1)(cdum,j=1,possp),(i2buf(j,ispec),j=1,npt)
            do j=npt,1,-1
               r4buf(j,ispec)=float(i2buf(j,ispec))
            end do
         else
            stop 'unknown BPW value'
         end if
         close(lunr_bs)
      end do   !  do ispec=1,2

c  Check that the spectral regions match
      if(ifirst(1).ne.ifirst(2)) then
         write(*,*) 'IFIRST = ',ifirst(1),ifirst(2)
         stop 'IFIRSTs differ'
      endif
      if(ilast(1).ne.ilast(2)) then
         write(*,*) 'ILAST = ',ilast(1),ilast(2)
         stop 'ILASTs differ'
      endif

c  Compute spectral differences
      td=0.0d0
      tb=0.0d0
      ncount=0
      do j=1,npt
         diff=abs(r4buf(j,2)-r4buf(j,1))
         if(diff.gt.0.0) ncount=ncount+1
         td=td+diff
         tb=tb+abs(r4buf(j,1))
      end do
      write(*,*)
     &' Spectrum     NDIFFS     NPTS     RMS Value       RMS Diff    
     &  Ratio: Diff/Val'
      write(*,'(a30,2(1x,i7),1x,f11.7,2(1x,1pe13.7))') 
     & specpath(1)(ls+1:lnblnk(specpath(1))),ncount,npt,
     & sngl(tb/npt),sngl(td/npt),sngl(td/tb)
      stop
      end
