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
      include "../ggg_int_params.f"

      integer*4 lunr,lunr_qc,lunw,lunw_csv,lunr_oof,
     & ncoml,ncol,mcol,kcol,icol,j,lg,ncoml_head,
     & lnbc,nrow,li,nco,k,kmax,krow_qc,irow,lof0,lof1,
     & eflag,wrow_flag,wsp_flag,nhead_gaa,nrow_gaa,ngas_gaa,lunr_gaa,
     & spectrum_flag,lh,le,klat,klong,kzobs,ncol_written,
     & jj,nflag,ncoml_qc,ncol_qc,nrow_qc,wcol_flag
      parameter (lunr=14,lunr_qc=15,lunw=16,lunw_csv=17,lunr_oof=18,
     & lunr_gaa=20,mcol=150)
      integer*4 flag(mrow_qc),pindex(mcol),kflag(mrow_qc)
      character
     & dl*1,
     & gggdir*(mpath),
     & version*62,
     & specname*(nchar),
     & csvfile*80,
     & header*800,
     & headarr(mcol)*20,
     & parname(mrow_qc)*20,
     & inputfile*80,
     & outputfile*80,
     & outfmt0*400,
     & outfmt1*400,
     & ofmt(mrow_qc)*4,
     & temp*20,
     & fmt(mrow_qc)*4,
     & unit(mrow_qc)*6,
     & headout*800,
     & headoutcsv*800,
     & sitefile*200,
     & ssss*800,
     & sarr(mcol)*57,
     & cc*20
      real*4 yrow(mcol),dev,dmax,scale(mrow_qc),
     & vmin(mrow_qc),vmax(mrow_qc),ymiss
      real*8 wlimit

      data kflag/mrow_qc*0/
      
      version=
     & ' write_official_output_file   Version 1.3.3   2010-12-04   GCT'
      write(*,*) version

      spectrum_flag=0    ! initialise to avoid compiler warnings
      klat=0
      klong=0
      kzobs=0
      wrow_flag=1
      wsp_flag=1

      call get_ggg_environment(gggdir, dl)
      lg=lnbc(gggdir)

      write(*,*)' Name of input file (e.g. paIn_1.0lm.vav.ada.aia.gaa):'
      read(*,'(a)') inputfile
      li=lnbc(inputfile)
      if(inputfile(li-3:li).ne.'.aia'.and.
     & inputfile(li-3:li).ne.'.gaa') write(*,*)
     & 'Warning: input file is not of expected type (.aia or .gaa)'
      outputfile=inputfile(:li)//'.oof'
      csvfile=inputfile(:li)//'.oof.csv'

c  Find length of site_oof_header.dat file
      sitefile=gggdir(:lg)//'tccon'//dl//inputfile(1:2)
     & //'_oof_header.dat'
      open(lunr_oof, file=sitefile, status='old')
      do j=1,9999
         read(lunr_oof,'(a)',end=66),ssss
      end do
66    ncoml_head=j-1
      close(lunr_oof)
      write(*,'(a,a,i4,a)') sitefile(:lnbc(sitefile)),
     & ' contains',ncoml_head,' lines'

c  Open the Quality Control (QC) file and read in the information,
c  to find out how many columns there are going to be in the .oof output files
      open(lunr_qc,file=
     & gggdir(:lg)//'tccon'//dl//inputfile(1:2)//'_qc.dat',
     & status='old')
      read(lunr_qc,*)ncoml_qc,ncol_qc,nrow_qc
      if(nrow_qc.gt.mrow_qc) stop 'nrow_qc > mrow_qc'
      read(lunr_qc,*)wrow_flag   ! 0/1   whether to skip/write the out-of-range data records
      read(lunr_qc,*)wsp_flag    ! 0/1   whether to skip/write the spectrum names (if present in .aia file)
      do k=4,ncoml_qc
         read(lunr_qc,'(a)') header
      end do
      ncol_written=1+wsp_flag   ! Always include "flag" in column #1
      do krow_qc=1,mrow_qc
         read(lunr_qc,*,end=77) temp,wcol_flag
         ncol_written=ncol_written+wcol_flag
      end do
77    close(lunr_qc)
      if(krow_qc-1.ne.nrow_qc)
     &  stop 'misreading xx_qc.dat file: krow_qc > nrow_qc'

c  Open the Ghost Correction file and figure out how many lines are
c  contained
      

c  Read input file and start writing output files
      open(lunr_gaa,file=
     & gggdir(:lg)//'tccon'//dl//inputfile(1:2)//'_ghost_corr.dat',
     & status='old')
      read(lunr_gaa,*)nhead_gaa,ngas_gaa,nrow_gaa
      close(lunr_gaa)

c     write(*,*)'inputfile=',inputfile
      open(lunr,file=inputfile, status='old')
      read(lunr,'(i2,i4,i7)') ncoml,ncol,nrow
      if(ncol.gt.mcol) stop 'increase mcol'
      open(lunw,file=outputfile,status='unknown')
      open(lunw_csv,file=csvfile,status='unknown')
      open(lunr_qc,file=
     & gggdir(:lg)//'tccon'//dl//inputfile(1:2)//'_qc.dat',
     & status='old')
      read(lunr_qc,*)ncoml_qc,ncol_qc,nrow_qc
      write(lunw,*) ncoml_head+ncoml_qc-0+ncoml+nrow_qc+4,ncol_written
      write(lunw_csv,*) ncoml_head+ncoml_qc-0+ncoml+nrow_qc+4,
     &       ncol_written
      write(lunw,'(a)') version
      write(lunw_csv,'(a)') version
      do j=2,ncoml-1
         read(lunr,'(a)') header
         write(lunw,'(a)') header(:lnbc(header))
         write(lunw_csv,'(a)') header(:lnbc(header))
         if (j.eq.ncoml-(nhead_gaa+nrow_gaa-2)) then
           read(header(:lnbc(header)),*)ymiss
         endif
c        if (j.eq.ncoml-4) write(*,*)'ymiss=',ymiss
      end do
      read(lunr,'(a)') header ! column headers
      call substr(header,headarr,mcol,kcol)
      if(kcol.ne.ncol) stop 'ncol/kcol mismatch'
      if (index(header,'spectrum') .gt. 0) spectrum_flag=1

      open(lunr_oof,
     & file=gggdir(:lg)//'tccon'//dl//inputfile(1:2)//'_oof_header.dat',
     & status='old')
      do j=1,ncoml_head
         read(lunr_oof,'(a)'),ssss
         write(lunw,'(a)') ssss(:lnbc(ssss))
         write(lunw_csv,'(a)') ssss(:lnbc(ssss))
      end do
      close(lunr_oof)

      do k=2,ncoml_qc
         read(lunr_qc,'(a)') header
         if(k.eq.ncoml_qc) header=' #'//header
         write(lunw,'(a)') header(:lnbc(header))
         write(lunw_csv,'(a)') header(:lnbc(header))
      end do

      do krow_qc=1,nrow_qc
         read(lunr_qc,'(a)') ssss
         write(lunw,'(i2,a)') krow_qc,ssss(:lnbc(ssss))
         write(lunw_csv,'(i2,a)') krow_qc,ssss(:lnbc(ssss))
         read(ssss,*) parname(krow_qc),flag(krow_qc),scale(krow_qc),
     &   fmt(krow_qc),unit(krow_qc),vmin(krow_qc),vmax(krow_qc)
      end do
      close(lunr_qc)

c  Read the header of the .aia file and figure out the
c  mapping between the gases in the corrections.dat
c  and those in the .vav file header
      do icol=1,ncol
         pindex(icol)=0
         temp=headarr(icol)
         if(temp(1:4).eq.'xao2') temp='x'//temp(3:)
         if(temp(1:4).eq.'xbo2') temp='x'//temp(3:)
         do krow_qc=1,nrow_qc
            if( temp .eq. parname(krow_qc) ) pindex(icol)=krow_qc
         end do
         if(pindex(icol).eq.0) then
           if(index(headarr(icol),'spectrum').eq.0) then
           write(*,*) ' Parameter missing from QC file: '//headarr(icol)
         endif
         endif
         if(headarr(icol).eq.'lat') klat=icol
         if(headarr(icol).eq.'long') klong=icol
         if(headarr(icol).eq.'zobs') kzobs=icol
      end do

      if(spectrum_flag.eq.0  .or. wsp_flag.eq.0) then
         headout=' flag'
         headoutcsv=' flag'
      else
         headout=' flag  spectrum'
         headoutcsv=' flag,spectrum'
      endif
      outfmt0='(i3'
      outfmt1='(i3,1x,a57'
      lof0=0
      lof1=0
      jj=0
      do icol=1+spectrum_flag,ncol
         krow_qc=pindex(icol)
         if(flag(krow_qc).ge.1) then
            lof0=lnbc(outfmt0)
            lof1=lnbc(outfmt1)
            jj=jj+1
            ofmt(jj)=fmt(krow_qc)
            outfmt0=outfmt0(:lof0)//',1x,'//fmt(krow_qc)
            outfmt1=outfmt1(:lof1)//',1x,'//fmt(krow_qc)
            le=index(headarr(icol),'_error')
            lh=lnbc(headarr(icol))
            if(le.eq.0) then
               headout=headout(:lnbc(headout))//'  '//
     &         headarr(icol)(:lh)//unit(krow_qc)
               headoutcsv=headoutcsv(:lnbc(headoutcsv))//','//
     &         headarr(icol)(:lh)//unit(krow_qc)
            else
               headout=headout(:lnbc(headout))//'  '//
     &         headarr(icol)(:le-1)//
     &         unit(krow_qc)(:lnbc(unit(krow_qc)))//'_error'
               headoutcsv=headoutcsv(:lnbc(headoutcsv))//','//
     &         headarr(icol)(:le-1)//
     &         unit(krow_qc)(:lnbc(unit(krow_qc)))//'_error'
c     &      headarr(icol)(:lnbc(headarr(icol)))//unit(krow_qc)
            endif
         endif
      end do
      outfmt0(lof0+9:lof0+9)=')'
      outfmt1(lof1+9:lof1+9)=')'

      write(lunw,'(120a1)') ('-',j=1,120)
      write(lunw_csv,'(120a1)') ('-',j=1,120)
c  Read each day of data into memory and multiply XGas values by the
c  appropriate correction factors.
      nflag=0
      do irow=1,9999999
         if (spectrum_flag .eq. 1) then
             read(lunr,*,end=99) specname, (yrow(j),j=2,ncol)
         else
c             write(*,*)'irow,ncol=',irow,ncol
             read(lunr,*,end=99) (yrow(j),j=1,ncol)
         endif
         if (irow .eq. 1) then
             write(lunw,'(a)') 'Latitude  Longitude  Altitude  SiteID'
             write(lunw_csv,'(a)')
     &             'Latitude  Longitude  Altitude  SiteID'
             write(lunw,'(3f9.3,4x,(a))')
     &       yrow(klat),yrow(klong),yrow(kzobs),inputfile(1:2)
             write(lunw_csv,'(3f9.3,4x,(a))')
     &       yrow(klat),yrow(klong),yrow(kzobs),inputfile(1:2)
             write(lunw,'(120a1)') ('-',j=1,120)
             write(lunw_csv,'(120a1)') ('-',j=1,120)
             write(lunw,'(a)') headout(:lnbc(headout))
             write(lunw_csv,'(a)') headoutcsv(:lnbc(headoutcsv))
         endif

c  Look within each data record to see if any of the data values are
c  outside their VMIN to VMAX range. If so, set eflag to the index of
c  the variable that was furthest out of range. Then write out the data.
         nco=0
         eflag=0
         kmax=0
         dmax=0.0
         do icol=1+spectrum_flag,ncol
            krow_qc=pindex(icol)
c           write(*,*)'scale,yrow,vmax,vmin,ymiss=',
c    & scale(krow_qc),yrow(icol),vmax(krow_qc),vmin(krow_qc),ymiss
            if (yrow(icol).ge.0.9*ymiss) then
c              write(*,*)'Missing field found.'
               dev=0.0 ! Don't include missing fields in the flags
            else
               dev=abs((scale(krow_qc)*yrow(icol)-vmin(krow_qc))/
     &         (vmax(krow_qc)-vmin(krow_qc))-0.5)
            endif
            if(dev.gt.dmax) then
               dmax=dev
               kmax=krow_qc
            endif
            if(flag(krow_qc).ge.1) then
               nco=nco+1
               if (yrow(icol).ge.0.9*ymiss) then
                  yrow(nco)=yrow(icol)
               else
                  yrow(nco)=yrow(icol)*scale(krow_qc)
               endif
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

         if(eflag.eq.0.or.wrow_flag.gt.0) then
            if(spectrum_flag.eq.0 .or. wsp_flag.eq.0) then
                write(lunw,outfmt0)eflag,
     &          (wlimit(dble(yrow(j)),ofmt(j)),j=1,nco)
            else
                write(lunw,outfmt1)eflag,specname,
     &          (wlimit(dble(yrow(j)),ofmt(j)),j=1,nco)
            endif
         endif
         if(spectrum_flag.eq.0 .or. wsp_flag.eq.0) then
            write(ssss,outfmt0)eflag,
     &      (wlimit(dble(yrow(j)),ofmt(j)),j=1,nco)
         else
            write(ssss,outfmt1)eflag,specname,
     &      (wlimit(dble(yrow(j)),ofmt(j)),j=1,nco)
         endif
         call substr(ssss,sarr,mcol,kcol)
c         if(kcol.ne.nco+wsp_flag+spectrum_flag) then
c            write(*,*)ssss
c            write(*,*)kcol,nco,wsp_flag,spectrum_flag
c            stop 'kcol.ne.nco+wsp_flag+spectrum_flag'
c         endif
         ssss=sarr(1)
         do k=2,kcol
            cc=sarr(k)(:20)
            ssss=ssss(:lnbc(ssss))//','//cc(:lnbc(cc))
         end do
         if(eflag.eq.0.or.wrow_flag.gt.0)
     &   write(lunw_csv,'(a)')ssss(:lnbc(ssss))
      end do         ! do irow=1,9999999
      stop ' irow exceeded 9999999'
99    close(lunr)
      close(lunw)
      close(lunw_csv)
      if(irow-1.ne.nrow) then
         write(*,*)'irow-1,nrow=',irow-1,nrow
         stop 'irow/nrow mismatch'
      endif
      write(*,*)nrow,' data records, of which',nflag,' flagged as bad'
      write(*,*)
      write(*,*)' Listed below are the fields with error flags exceeding
     & the allowed range and the number of such occurrences'
      write(*,*)
      write(*,*)' #   Parameter            N_flag     %'
      do k=1,nrow_qc
          if(kflag(k).gt.0) write(*,'(i3,3x,a,i6,f8.1)') k,
     &  parname(k), kflag(k), 100*float(kflag(k))/nrow
      end do
      stop
      end

