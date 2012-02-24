C  Program: apply_airmass_dependence.f
c
c  Purpose: To apply the airmass correction to the data in the
c  selected .vav filE
c  
c  yout(i,k) = yin(i,k) / [ 1 + ADCF(k)*SBF(i) ]
c
c  where
c   yin(i,k) is the column measured for the k'th gas from the i'th spectrum
c   yout(i,k) is the DMF of the k'th gas from the i'th spectrum
c   ADCF(j) is the Airmass-Dependent Correction Factor for the k'th gas
c   SBF(i) is the Symmetric Basis Function (symmetric about solar noon)
c   SBF(i) = ((SZA(i)+13)/(90+13))**3 - ((45+13)/(90+13))**3
c   SZA(i) is the Solar Zenith Angle of the i'th spectrum.
c   SBF is zero at SZA=45 deg, -0.122 at SZA=0, and 0.875 at SZA=77 deg.
c   Over the range of SZAs (0-80 deg) encounted at a site, SBF changes by ~1,
c   so ADCF is a good estimate of the size of the airmass correction.
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
      include "../ggg_int_params.f"

      integer*4 lunr,luns,lunw,ncoml,ncolvav,
     & kcolvav,jcol,j,
     & kgas,ko2,ksza,lnbc,irow,naux,mgas,ngas,nrow,li,k,ncolcorr
      parameter (lunr=14,luns=15,lunw=16,mgas=9)
      real*8 yrow(mcolvav),ymiss,
     & adcf(mgas),adcf_err(mgas),
     & aicf(mgas),aicf_err(mgas),
     & cf(mcolvav),fu,sbf,vc_air
      character
     & dl*1,
     & gggdir*(mpath),
     & version*62,
     & specname*57,
     & header*800,
     & headarr(mcolvav)*20,
     & gasname(mgas)*20,
     & inputfile*40,
     & outputfile*40,
     & gaserr*32,
     & output_fmt*40

      integer*4 specflag 
      specflag=0

      version=
     &' apply_airmass_correction     Version 1.1.7   2011-11-05   GCT'
      write(*,*) version
      call get_ggg_environment(gggdir, dl)
      ko2=0
      ksza=0

c  Open the corrections.dat file and read in the Airmass-Dependent
c  Correction Factors (ADCF) and the Airmass-Independent Correction
c  Factors (AICF) for each gas. Only the former is used by this prog.
      open(luns,file=gggdir(:lnbc(gggdir))//'tccon'//dl
     & //'corrections.dat', status='old')
      read(luns,*)ncoml,ncolcorr
      do k=2,ncoml
         read(luns,*)
      end do

      if(ncolcorr.eq.3) then
         do k=1,mgas
            read(luns,*,end=88) gasname(k),adcf(k),aicf(k)
         end do
      elseif(ncolcorr.eq.5) then
         do k=1,mgas
            read(luns,*,end=88) gasname(k),
     &   adcf(k),adcf_err(k),aicf(k),aicf_err(k)
         end do
      else
         write(*,*) 'ncolcorr=',ncolcorr
         stop 'Unrecognized ncolcorr value'
      endif
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

c  Read the header of the .vav file and copy it to the output
c  file, adding the correction coefficients
      read(lunr,'(i2,i4,i7,i4)') ncoml,ncolvav,nrow,naux
      write(lunw,'(i2,i4,i7,i4)') ncoml+1+ngas+1,ncolvav,nrow,naux
      write(lunw,'(a)') version
      if(ncolvav.gt.mcolvav) stop 'increase mcolvav'
      do j=2,ncoml-3
         read(lunr,'(a)') header
         write(lunw,'(a)') header(:lnbc(header))
      end do
      write(lunw,'(a)') ' Airmass-Dependent Correction Factors: '
      do k=1,ngas
         if(ncolcorr.eq.3) write(lunw,'(a,f9.4)') gasname(k),adcf(k)
         if(ncolcorr.eq.5) write(lunw,'(a,2f9.4)') gasname(k),adcf(k),
     &   adcf_err(k)
      end do
      read(lunr,*) ymiss
      ymiss=ymiss-0.0001E+29
      write(lunw,*)ymiss
      read(lunr,'(a)') output_fmt
      read(lunr,'(a)') header
      call substr(header,headarr,mcolvav,kcolvav)
      if (index(header,'spectrum') .gt. 0) specflag=1
      if (specflag.eq.1) output_fmt(6:7)='1x'
      write(lunw,'(a)') output_fmt
c      write(*,*) index(header,'spectrum') 
      do j=naux+1,ncolvav
        headarr(j)='x'//headarr(j)
      end do
      write(lunw,'(100a12)') (headarr(j),j=1,ncolvav)
      if(kcolvav.ne.ncolvav ) stop 'ncolvav/kcolvav mismatch'

c  Figure out the mapping between the gases in the corrections.dat
c  and those in the .vav file header
      do jcol=1,ncolvav
         cf(jcol)=0.0
         do kgas=1,ngas
            if(headarr(jcol).eq.gasname(kgas)) cf(jcol)=adcf(kgas)
            gaserr=gasname(kgas)(:lnbc(gasname(kgas)))//'_error'
            if(headarr(jcol).eq.gaserr) cf(jcol)=adcf(kgas)
         end do
         if(headarr(jcol) .eq. 'xo2') ko2=jcol
         if(headarr(jcol) .eq. 'asza') ksza=jcol
         if(jcol.gt.naux) write(*,'(i3,f8.4,2x,a)')
     &   jcol,cf(jcol),headarr(jcol)
      end do

      if(ko2.eq.0) stop ' o2 column not found'
      if(ksza.eq.0) stop ' asza column not found'

c      if (specflag .eq. 1) then
c         output_fmt='(a38,f14.8,NNf13.5,200(1pe12.4))'
c         write(output_fmt(12:13),'(i2.2)') naux-1-specflag
c         write(output_fmt(20:22),'(i3.3)') ncolvav-naux
c      else 
c         output_fmt='(f14.8,NNf13.5,MMM(1pe12.4))'
c         write(output_fmt(8:9),'(i2.2)') naux-1-specflag
c         write(output_fmt(16:18),'(i3.3)') ncolvav-naux
c      endif

c  Read each day of data into memory.
      do irow=1,9999999
         if (specflag .eq. 1) then 
             read(lunr,*,end=99) specname, (yrow(j),j=2,ncolvav)
         else
             read(lunr,*,end=99) (yrow(j),j=1,ncolvav)
         endif
         vc_air=yrow(ko2)/0.2095
         if(yrow(ko2).eq.0.0) then
           write(*,'(a,i6,a)')'Warning: O2_column=0 for spectrum',irow,
     &     '  (lamp run?)  Your output file will contain Inf'
         endif
         sbf=((yrow(ksza)+13)/(90+13))**3-((45.+13)/(90+13))**3  ! Symmetric Basis Function
         do k=naux+1,ncolvav-1,2
            if(yrow(k).lt.ymiss) then
            if(k.eq.ko2) then
              fu=yrow(k+1)/yrow(k)                 ! fractional uncertainty
            else
              fu=sqrt((yrow(k+1)/yrow(k))**2+(yrow(ko2+1)/yrow(ko2))**2)
            endif
            yrow(k)=yrow(k)/vc_air/(1+cf(k)*sbf)   ! apply airmass correction
            yrow(k+1)=yrow(k)*fu
            endif  ! if(yrow(k).lt.ymiss
         end do
         if (specflag .eq. 1) then
            write(lunw,output_fmt) specname, (yrow(j),j=2,ncolvav)
         else
            write(lunw,output_fmt) ' ',(yrow(j),j=1,ncolvav)
         endif
      end do      ! do irow=1,9999999
      stop ' irow exceeded 9999999'
99    close (lunr)
      close (lunw)
      stop
      end
