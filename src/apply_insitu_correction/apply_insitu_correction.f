c  Program: apply_insitu_correction.f
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
      include "../ggg_int_params.f"

      integer*4 lunr,luns,lunw,ncoml,ncolvav,kcolvav,icol,j,
     & kgas,lnbc,irow,naux,mgas,ngas,nrow,li,k,ncolcorr
      parameter (lunr=14,luns=15,lunw=16,mgas=9)
      real*8 yrow(mcolvav),adcf(mgas),adcf_err(mgas),ymiss,
     & aicf(mgas),aicf_err(mgas),cf(mcolvav)

      character header*800,headarr(mcolvav)*40,gasname(mgas)*20,dl*1,
     & gggdir*(mpath),inputfile*40,outputfile*40, version*62,gaserr*32
      character output_fmt*40, input_fmt*40, specname*(nchar), c1*1
      
      integer specflag
      specflag=0

      version=
     & ' apply_insitu_correction      Version 1.3.4   2011-11-05   GCT'

      call get_ggg_environment(gggdir, dl)

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
     &      adcf(k),adcf_err(k),aicf(k),aicf_err(k)
         end do
      else
         write(*,*)'ncol=',ncolcorr
         stop 'Unrecognized NCOLCORR value'
      endif
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
c  and those in the .vav.ada file header
      read(lunr,'(i2,i4,i7,i4)') ncoml,ncolvav,nrow,naux
      write(lunw,'(i2,i4,i7,i4)') ncoml+1+ngas+1,ncolvav,nrow,naux
      write(lunw,'(a)') version
      if(ncolvav.gt.mcolvav) stop 'increase mcolvav'
      do j=2,ncoml-3
         read(lunr,'(a)') header
         write(lunw,'(a)') header(:lnbc(header))
      end do
      
      write(lunw,'(a)')'Airmass-Independent/In-Situ Correction Factors:'
      do k=1,ngas
         if(ncolcorr.eq.3) write(lunw,'(a,f9.4)') gasname(k),aicf(k)
         if(ncolcorr.eq.5) write(lunw,'(a,2f9.4)') gasname(k),aicf(k),
     &   aicf_err(k)
      end do

      read(lunr,*) ymiss
      write(lunw,*) ymiss
      read(lunr,'(a)') input_fmt
      write(lunw,'(a)') input_fmt
      read(lunr,'(a)') header
      write(lunw,'(a)') header(:lnbc(header))
      if (index(header,'spectrum') .gt. 0) specflag=1
      call substr(header,headarr,mcolvav,kcolvav)
      if(kcolvav.ne.ncolvav ) stop 'ncolvav/kcolvav mismatch'
      do icol=1,ncolvav
         cf(icol)=1.0
         do kgas=1,ngas
c         write(*,*)kgas,icol,gasname(kgas),headarr(icol),aicf(kgas),cf(icol)
           gaserr=gasname(kgas)(:lnbc(gasname(kgas)))//'_error'
           if( headarr(icol) .eq. gasname(kgas) ) cf(icol)=aicf(kgas)
           if( headarr(icol) .eq. gaserr        ) cf(icol)=aicf(kgas)
         end do
         if(icol.gt.naux+specflag) write(*,'(i4,a16,f8.3)')
     &   icol,headarr(icol),cf(icol)
      end do

c      if (specflag .eq. 1) then
c         output_fmt='(a35,f14.8,NNf13.5,200(1pe12.4))'
c         write(output_fmt(12:13),'(i2.2)') naux-1
c      else
c         output_fmt='(f14.8,NNf13.5,200(1pe12.4))'
c         write(output_fmt(8:9),'(i2.2)') naux-2
c      endif
      output_fmt=input_fmt

c  Read each day of data into memory and divide XGas values by the
c  appropriate correction factors.
      do irow=1,9999999
         if (specflag .eq. 1) then
            read(lunr,input_fmt,end=99) specname, (yrow(j),j=2,ncolvav)
         else
            read(lunr,input_fmt,end=99) c1,(yrow(j),j=1,ncolvav)
         endif
         do k=naux+1,ncolvav-1,2
              if(yrow(k).lt.ymiss) then  ! Dont correct missing values
                 yrow(k)=yrow(k)/cf(k)
                 yrow(k+1)=yrow(k+1)/cf(k)
              endif  !  (yrow(k).lt.ymiss) then
         end do
         if (specflag .eq. 1) then
            write(lunw,output_fmt) specname,(yrow(j),j=2,ncolvav)
         else
            write(lunw,output_fmt) c1,(yrow(j),j=1,ncolvav)
         endif
      end do         ! do irow=1,9999999
      stop ' irow exceeded 9999999'
99    close (lunr)
      close(lunw)
      stop
      end
