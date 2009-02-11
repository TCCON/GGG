c  Program: derive_airmass_dependence.f
c
c  Purpose: To apply the (airmass-independent) in situ correction to the data
c  in the selected runlog.vav.ada file
c  
c  yout(i,j) = AICF(j) * yin(i,j)
c
c  where
c   yin(i,j) is the DMF of the j'th gas from the i'th spectrum
c   yout(i,j) is the DMF of the j'th gas from the i'th spectrum
c   AICF(j) is the Airmass-Independent Correction Factor for the j'th gas
c   which is read from the file ~/ggg/tccon/corrections.dat
c
c  Input Files:
c       runlog.vav.ada 
c       corrections.dat  
c
c  Output Files:
c       runlog.vav.ada.aia
c       
      implicit none
      integer*4 lunr,luns,lunw,ncoml,ncol,mcol,kcol,icol,j,
     & kgas,lnbc,irow,naux,mgas,ngas,nrow,li,k
      parameter (lunr=14,luns=15,lunw=16,mcol=150,mgas=9)
      character header*800,headarr(mcol)*20,gasname(mgas)*20,
     & gggdir*80,inputfile*40,outputfile*40, version*60,gaserr*32
      real*8 yrow(mcol),adcf(mgas),aicf(mgas),cf(mcol)

      version=
     & ' apply_insitu_correction    Version 1.0.4   2009-02-04   GCT'

      call getenv('GGGPATH',gggdir)

c  Open the corrections.dat file and read in the Airmass-Dependent
c  Correction Factors (ADCF) and the Airmass-Independent Correction
c  Factors (AICF) for each gas. Only the former is used by this prog.
      open(luns,file=gggdir(:lnbc(gggdir))//'/tccon/corrections.dat',
     & status='old')
      read(luns,*)ncoml,ncol
      do k=2,ncoml
         read(luns,*)
      end do
      do k=1,mgas
         read(luns,*,end=88) gasname(k),adcf(k),aicf(k)
      end do
      stop 'increase parameter MGAS'
88    ngas=k-1

      write(*,*)'Enter name of input file (e.g. paIn_1.0lm.vav.ada):'
      read(*,'(a)') inputfile
      li=lnbc(inputfile)
      if(inputfile(li-3:li).ne.'.ada') write(*,*)
     & 'Warning: input file is not of expected type (.ada)'
      outputfile=inputfile(:li)//'.aia'
      open(lunr,file=inputfile, status='old')
      open(lunw,file=outputfile,status='unknown')

c  Read the header of the .ada file and figure out the
c  mapping between the gases in the corrections.dat
c  and those in the .vav file header
      read(lunr,'(i2,i4,i7,i4)') ncoml,ncol,nrow,naux
      write(lunw,'(i2,i4,i7,i4)') ncoml+1,ncol,nrow,naux
      write(lunw,'(a)') version
      if(ncol.gt.mcol) stop 'increase mcol'
      do j=2,ncoml
         read(lunr,'(a)') header
         write(lunw,'(a)') header(:lnbc(header))
      end do
      call substr(header,headarr,mcol,kcol)
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'
      do icol=1,ncol
         cf(icol)=1.0
         do kgas=1,ngas
c         write(*,*)kgas,icol,gasname(kgas),headarr(icol),aicf(kgas),cf(icol)
           gaserr=gasname(kgas)(:lnbc(gasname(kgas)))//'_error'
           if( headarr(icol) .eq. gasname(kgas) ) cf(icol)=aicf(kgas)
           if( headarr(icol) .eq. gaserr        ) cf(icol)=aicf(kgas)
         end do
         if(icol.gt.naux) write(*,'(i4,a16,f8.3)')
     &   icol,headarr(icol),cf(icol)
      end do

c  Read each day of data into memory and divide XGas values by the
c  appropriate correction factors.
      do irow=1,9999999
         read(lunr,*,end=99) (yrow(j),j=1,ncol)
         do k=naux+1,ncol-1,2
              yrow(k)=yrow(k)/cf(k)
              yrow(k+1)=yrow(k+1)/cf(k)
         end do
         write(lunw,65) (yrow(j),j=1,ncol)
65       format (f14.8,18f13.5,200(1pe12.4))
      end do         ! do irow=1,9999999
      stop ' irow exceeded 9999999'
99    close (lunr)
      close(lunw)
      stop
      end
