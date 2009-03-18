c  Program: create_official_output_file.f
c
c  Purpose: To convert the runlog.vav.ada.aia file to
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
      integer*4 lunr,luns,lunw,lunc,lunh,
     & ncoml,ncol,mcol,kcol,icol,j,lg,ncoml_head,
     & lnbc,nrow,li,nco,k,kmax,mpar,kpar,npar,irow,lo,eflag,ii,
     & nchar,
     & jj,nflag,ncoml_qc,ncol_qc,nrow_qc
      parameter (lunr=14,luns=15,lunw=16,lunc=17,lunh=18,mcol=150,
     & mpar=150)
      integer*4 flag(mpar),pindex(mcol),kflag(mpar)
      character header*800,headarr(mcol)*20,parname(mpar)*20,gggdir*80,
     & inputfile*80,outputfile*80,csvfile*80,version*62,outfmt*240,
     & ofmt(mpar)*4,temp*20, fmt(mpar)*4,unit(mpar)*6,headout*800,
     & specname,
     & ssss*800,sarr(mcol)*20,cc*20
      real*4 yrow(mcol),dev,dmax,scale(mpar),
     & vmin(mpar),vmax(mpar)
      real*8 wlimit

      data kflag/mpar*0/
      version=
     & ' write_official_output_file   Version 1.1.1   2009-03-04   GCT'
      write(*,*) version

      nchar=0

      call getenv('GGGPATH',gggdir)
      lg=lnbc(gggdir)

      write(*,*)' Name of input file (e.g. paIn_1.0lm.vav.ada.aia):'
      read(*,'(a)') inputfile
      li=lnbc(inputfile)
      if(inputfile(li-3:li).ne.'.aia') write(*,*)
     & 'Warning: input file is not of expected type (.aia)'
      outputfile=inputfile(:li)//'.oof'
      csvfile=inputfile(:li)//'.oof.csv'

c  Find length of site_oof_header.dat file
      open(lunh,
     & file=gggdir(:lg)//'/tccon/'//inputfile(1:2)//'_oof_header.dat',
     & status='old')
      do j=1,999
         read(lunh,'(a)',end=66),ssss
      end do
66    ncoml_head=j-1
      close(lunh)

c  Open the Quality Control (QC) file and read in the information
      open(lunr,file=inputfile, status='old')
      read(lunr,'(i2,i4,i7)') ncoml,ncol,nrow
      if(ncol.gt.mcol) stop 'increase mcol'
      open(lunw,file=outputfile,status='unknown')
      open(lunc,file=csvfile,status='unknown')
      open(luns,file=gggdir(:lg)//'/tccon/'//inputfile(1:2)//'_qc.dat',
     & status='old')
      read(luns,*)ncoml_qc,ncol_qc,nrow_qc
      write(lunw,*) ncoml_head+ncoml_qc+ncoml+nrow_qc,ncol
      write(lunc,*) ncoml_head+ncoml_qc+ncoml+nrow_qc,ncol
      write(lunw,'(a)') version
      write(lunc,'(a)') version
      do j=2,ncoml-2
         read(lunr,'(a)') header
         write(lunw,'(a)') header(:lnbc(header))
         write(lunc,'(a)') header(:lnbc(header))
      end do
      read(lunr,'(a)') header ! missing values
      read(lunr,'(a)') header ! column headers
      call substr(header,headarr,mcol,kcol)
      if(kcol.ne.ncol) stop 'ncol/kcol mismatch'
      if (index(header,'Spectrum') .gt. 0) nchar=1

      open(lunh,
     & file=gggdir(:lg)//'/tccon/'//inputfile(1:2)//'_oof_header.dat',
     & status='old')
      do j=1,ncoml_head
         read(lunh,'(a)'),ssss
         write(lunw,'(a)') ssss(:lnbc(ssss))
         write(lunc,'(a)') ssss(:lnbc(ssss))
      end do
      close(lunh)

      do k=2,ncoml_qc
         read(luns,'(a)') header
         write(lunw,'(a)') header(:lnbc(header))
         write(lunc,'(a)') header(:lnbc(header))
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
      outfmt='(i3'
      lo=0
      jj=0
      do icol=1+nchar,ncol
         kpar=pindex(icol)
         if(flag(kpar).ge.1) then
            lo=lnbc(outfmt)
            jj=jj+1
            ofmt(jj)=fmt(kpar)
            outfmt=outfmt(:lo)//',1x,'//fmt(kpar)
            headout=headout(:lnbc(headout))//'  '//
     &      headarr(icol)(:lnbc(headarr(icol)))//unit(kpar)
         endif
      end do
      outfmt(lo+9:lo+9)=')'

      write(lunw,'(a)') headout(:lnbc(headout))
      write(lunc,'(a)') headout(:lnbc(headout))

      write(lunw,'(a)') '-----------------------------------------------
     &-----------------------------------------------------------------'
      write(lunc,'(a)') '-----------------------------------------------
     &-----------------------------------------------------------------'
c  Read each day of data into memory and multiply XGas values by the
c  appropriate correction factors.
      nflag=0
      do irow=1,9999999
         if (nchar .eq. 1) then
             read(lunr,*,end=99) specname, (yrow(j),j=1+nchar,ncol)
         else
             read(lunr,*,end=99) (yrow(j),j=1,ncol)
         endif
         nco=0
         eflag=0
         kmax=0
         dmax=0.0
         do icol=1+nchar,ncol
            kpar=pindex(icol)
            dev=abs((scale(kpar)*yrow(icol)-vmin(kpar))/
     &      (vmax(kpar)-vmin(kpar))-0.5)
            if(dev.gt.dmax) then
               dmax=dev
               kmax=kpar
            endif
            if(flag(kpar).ge.1) then
               nco=nco+1
               yrow(nco)=yrow(icol)*scale(kpar)
            endif
         end do  ! do icol=1,ncol
         if(dmax.gt.0.5) then
            eflag=kmax
            nflag=nflag+1
            kflag(kmax)=kflag(kmax)+1
         endif
c         yrow(1)=int(yrow(1))   ! year
c         yrow(2)=int(yrow(2))   ! Day
         yrow(1)=nint(yrow(1)-yrow(2)/365.25)   ! Year
         yrow(2)=nint(yrow(2)-yrow(3)/24.0)     ! Day

         write(lunw,outfmt)eflag,(wlimit(dble(yrow(j)),ofmt(j)),j=1,nco)
         write(ssss,outfmt)eflag,(wlimit(dble(yrow(j)),ofmt(j)),j=1,nco)
         call substr(ssss,sarr,mcol,kcol)
         if(kcol.ne.nco+1) stop 'kcol.ne.nco+1'
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
      if(irow-1.ne.nrow) stop 'nrow mismatch'
      write(*,*)nrow,' data records, of which',nflag,' flagged as bad'
      write(*,*)
      write(*,*)' Listed below are the fields with error flags exceeding
     & the allowed range and the number of such occurrences'
      write(*,*)
      write(*,*)' #   Parameter            N_flag     %'
      do k=1,npar
          if(kflag(k).gt.0) write(*,'(i3,3x,a,i6,f8.1)') k,
     &  parname(k), kflag(k), 100*float(kflag(k))/nrow
      end do
      stop
      end

