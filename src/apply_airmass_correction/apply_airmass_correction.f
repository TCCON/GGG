c  Program: derive_airmass_dependence.f
c
c  Purpose: To apply the airmass correction to the data in the
c  selected .vav file
c  
c  yout(i,k) = yin(i,k) / [ 1 + ADCF(k)*SBF(i) ]
c
c  where
c   yin(i,k) is the column measured for the k'th gas from the i'th spectrum
c   yout(i,k) is the DMF of the k'th gas from the i'th spectrum
c   ADCF(j) is the Airmass-Dependent Correction Factor for the k'th gas
c   SBF(i) = ((SZA(i)+13)/(90+13))**3 - ((45+13)/(90+13))**3
c   SZA(i) is the Solar Zenith Angle of the i'th spectrum.
c
c  SBF is zero at SZA=45 deg, is -0.122 at SZA=0, is 0.875 at SZA=77
c
c
c  Input Files:
c       runlog.vav 
c       corrections.dat  
c
c  Output Files:
c       runlog.vav.ada
c       
      implicit none
      integer*4 lunr,luns,lunw,ncoml,ncol,mcol,kcol,jcol,j,
     & kgas,ko2,ksza,lnbc,irow,naux,mgas,ngas,nrow,li,k
      parameter (lunr=14,luns=15,lunw=16,mcol=150,mgas=9)
      character header*800,headarr(mcol)*20,gasname(mgas)*20,gggdir*80,
     & inputfile*40,outputfile*40, version*62,gaserr*32,output_fmt*32
      real*8 yrow(mcol),adcf(mgas),aicf(mgas),cf(mcol),fu,sbf,vc_air
      character specname*38

      integer nchar 
      nchar=0

      version=
     &' apply_airmass_correction     Version 1.1.2   2009-11-07   GCT'
      write(*,*) version
      call getenv('GGGPATH',gggdir)
      ko2=0
      ksza=0

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

      write(*,*)'Enter name of input file (e.g. paIn_1.0lm.vav):'
      read(*,'(a)') inputfile
      li=lnbc(inputfile)
      if(inputfile(li-3:li).ne.'.vav') write(*,*)
     &  ' Warning: input file is not of expected type (.vav) '
      outputfile=inputfile(:li)//'.ada'
      open(lunr,file=inputfile, status='old')
      open(lunw,file=outputfile,status='unknown')

c  Read the header of the .vav file and figure out the
c  mapping between the gases in the corrections.dat
c  and those in the .vav file header
      read(lunr,'(i2,i4,i7,i4)') ncoml,ncol,nrow,naux
      write(lunw,'(i2,i4,i7,i4)') ncoml+1+ngas+1,ncol,nrow,naux
      write(lunw,'(a)') version
      if(ncol.gt.mcol) stop 'increase mcol'
      do j=2,ncoml-1
         read(lunr,'(a)') header
         write(lunw,'(a)') header(:lnbc(header))
      end do
      write(lunw,'(a)') ' Airmass-Dependent Correction Factors: '
      do k=1,ngas
         write(lunw,'(a,f9.4)') gasname(k),adcf(k)
      end do
      read(lunr,'(a)') header
      call substr(header,headarr,mcol,kcol)
      if (index(header,'Spectrum') .gt. 0) nchar=1
      write(*,*) index(header,'Spectrum') 
      do j=naux+nchar+1,ncol
        headarr(j)='x'//headarr(j)
      end do
      write(lunw,'(100a12)') (headarr(j),j=1,ncol)
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'
      do jcol=1,ncol
         cf(jcol)=0.0
         do kgas=1,ngas
            if(headarr(jcol).eq.gasname(kgas)) cf(jcol)=adcf(kgas)
            gaserr=gasname(kgas)(:lnbc(gasname(kgas)))//'_error'
            if(headarr(jcol).eq.gaserr) cf(jcol)=adcf(kgas)
         end do
         if(headarr(jcol) .eq. 'xo2') ko2=jcol
         if(headarr(jcol) .eq. 'asza') ksza=jcol
         if(jcol.gt.naux+nchar) write(*,'(i3,f8.4,2x,a)')
     &   jcol,cf(jcol),headarr(jcol)
      end do

      if(ko2.eq.0) stop ' o2 column not found'
      if(ksza.eq.0) stop ' asza column not found'

      if (nchar .eq. 1) then
         output_fmt='(a35,f14.8,NNf13.5,200(1pe12.4))'
         write(output_fmt(12:13),'(i2.2)') naux-1
      else 
         output_fmt='(f14.8,NNf13.5,200(1pe12.4))'
         write(output_fmt(8:9),'(i2.2)') naux-1
      endif
c  Read each day of data into memory.
      do irow=1,9999999
         if (nchar .eq. 1) then 
             read(lunr,*,end=99) specname, (yrow(j),j=1+nchar,ncol)
         else
             read(lunr,*,end=99) (yrow(j),j=1,ncol)
         endif
         vc_air=yrow(ko2)/0.2095
         if(yrow(ko2).eq.0.0) then
           write(*,'(a,i6,a)')'Warning: O2_column=0 for spectrum',irow,
     &     '  (lamp run?)  Your output file will contain Inf'
         endif
         sbf=((yrow(ksza)+13)/(90+13))**3-((45.+13)/(90+13))**3  ! Symmetric Basis Function
         do k=naux+nchar+1,ncol-1,2
            if(k.eq.ko2) then
              fu=yrow(k+1)/yrow(k)                 ! fractional uncertainty
            else
              fu=sqrt((yrow(k+1)/yrow(k))**2+(yrow(ko2+1)/yrow(ko2))**2)
            endif
            yrow(k)=yrow(k)/vc_air/(1+cf(k)*sbf)   ! apply airmass correction
            yrow(k+1)=yrow(k)*fu
         end do
         if (nchar .eq. 1) then
            write(lunw,output_fmt) specname, (yrow(j),j=1+nchar,ncol)
         else
            write(lunw,output_fmt) (yrow(j),j=1,ncol)
         endif
      end do      ! do irow=1,9999999
      stop ' irow exceeded 9999999'
99    close (lunr)
      close (lunw)
      stop
      end
