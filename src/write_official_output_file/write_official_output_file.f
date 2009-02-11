c  Program: create_official_output_file.f
c
c  Purpose: To conver the runlog.vav.ada.aia file to
c  an official format.
c  
c  Input Files:
c       runlog.vav.ada.aia 
c       qc_limits.dat  
c
c  Output Files:
c       runlog.vav.ada.aia.oof
c       
      implicit none
      integer*4 lunr,luns,lunw,lunc,ncoml,ncol,mcol,kcol,icol,j,lg,
     & lnbc,nrow,li,ncolout,k,kmax,mpar,kpar,npar,irow,lo,eflag,ii
      parameter (lunr=14,luns=15,lunw=16,lunc=17,mcol=150,mpar=150)
      integer*4 flag(mpar),pindex(mcol)
      character header*800,headarr(mcol)*20,parname(mpar)*20,gggdir*80,
     & inputfile*80,outputfile*80,csvfile*80,version*62,
     & outputfmt*150,temp*20, fmt(mpar)*4,unit(mpar)*6,headout*800,
     & ssss*800,sarr(mcol)*20,cc*20
      real*4 yrow(mcol),dev,dmax,scale(mpar),
     & vmin(mpar),vmax(mpar)

      version=
     & ' write_official_output_file   Version 1.0.4   2009-02-06   GCT'

      call getenv('GGGPATH',gggdir)
      lg=lnbc(gggdir)

      write(*,*)' Name of input file (e.g. paIn_1.0lm.vav.ada.aia):'
      read(*,'(a)') inputfile
      li=lnbc(inputfile)
      if(inputfile(li-3:li).ne.'.aia') write(*,*)
     & 'Warning: input file is not of expected type (.aia)'
      outputfile=inputfile(:li)//'.oof'
      csvfile=inputfile(:li)//'.oof.csv'

c  Open the Quality Control (QC) file and read in the information
      open(lunr,file=inputfile, status='old')
      open(lunw,file=outputfile,status='unknown')
      open(lunc,file=csvfile,status='unknown')
      open(luns,file=gggdir(:lg)//'/tccon/'//inputfile(1:2)//'_qc.dat',
     & status='old')
      read(luns,*)ncoml,ncol
      do k=2,ncoml
         read(luns,*)
      end do
      do kpar=1,mpar
         read(luns,'(a)',end=88) ssss
         write(lunw,'(a)') ssss(:lnbc(ssss))
         write(lunc,'(a)') ssss(:lnbc(ssss))
         read(ssss,*) ii,parname(kpar),flag(kpar),scale(kpar),
     &   fmt(kpar),unit(kpar),vmin(kpar),vmax(kpar)
         if(ii.ne.kpar) write(*,*)'Warning: ii.ne.kgas ',ii,kpar
      end do
      stop 'increase parameter MGAS'
88    npar=kpar-1

c  Read the header of the .aia file and figure out the
c  mapping between the gases in the corrections.dat
c  and those in the .vav file header
      read(lunr,'(i2,i4,i7)') ncoml,ncol,nrow
      if(ncol.gt.mcol) stop 'increase mcol'
      do j=2,ncoml
         read(lunr,'(a)') header
      end do
      call substr(header,headarr,mcol,kcol)
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'
      do icol=1,ncol
         pindex(icol)=0
         do kpar=1,npar
            temp=headarr(icol)
            if(temp(1:4).eq.'xao2') temp='x'//temp(3:)
            if(temp(1:4).eq.'xbo2') temp='x'//temp(3:)
            if( temp .eq. parname(kpar) ) pindex(icol)=kpar
         end do
         if(pindex(icol).eq.0) write(*,*)
     &   ' Parameter missing from QC file: '//headarr(icol)
      end do

      headout=' flag'
      outputfmt='(i3'
      lo=0
      do icol=1,ncol
         kpar=pindex(icol)
         if(flag(kpar).ge.1) then
            lo=lnbc(outputfmt)
            outputfmt=outputfmt(:lo)//','//fmt(kpar)
            headout=headout(:lnbc(headout))//'  '//
     &      headarr(icol)(:lnbc(headarr(icol)))//unit(kpar)
         endif
      end do
      outputfmt(lo+6:lo+6)=')'

      write(lunw,'(a)') headout(:lnbc(headout))
      write(lunc,'(a)') headout(:lnbc(headout))

c  Read each day of data into memory and multiply XGas values by the
c  appropriate correction factors.
      do irow=1,9999999
         read(lunr,*,end=99) (yrow(j),j=1,ncol)
         ncolout=0
         eflag=0
         kmax=0
         dmax=0.0
         do icol=1,ncol
            kpar=pindex(icol)
            dev=abs((scale(kpar)*yrow(icol)-vmin(kpar))/
     &      (vmax(kpar)-vmin(kpar))-0.5)
            if(dev.gt.dmax) then
               dmax=dev
               kmax=kpar
            endif
            if(flag(kpar).ge.1) then
               ncolout=ncolout+1
               yrow(ncolout)=yrow(icol)*scale(kpar)
            endif
         end do  ! do icol=1,ncol
         if(dmax.gt.0.5) eflag=kmax
         yrow(1)=int(yrow(1))   ! year
         yrow(2)=int(yrow(2))   ! Day
         write(lunw,outputfmt) eflag,(yrow(j),j=1,ncolout)
         write(ssss,outputfmt) eflag,(yrow(j),j=1,ncolout)
         call substr(ssss,sarr,mcol,kcol)
         if(kcol.ne.ncolout+1) stop 'kcol.ne.ncolout+1'
         ssss=sarr(1)
         do k=2,kcol
            cc=sarr(k)
            ssss=ssss(:lnbc(ssss))//','//cc(:lnbc(cc))
         end do
         write(lunc,'(a)')ssss(:lnbc(ssss))
      end do         ! do irow=1,9999999
      stop ' irow exceeded 9999999'
99    close(lunr)
      close(lunw)
      close(lunc)
      stop
      end

