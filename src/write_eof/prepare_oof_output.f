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

c subroutine read_oneline_oof(inputfile,inputlun, oof_out,
c & outfmt1, oof_flag, irow,
c & vmin, vmax, pindex,
c & nflag, eflag, kflag, flag,
c & nrow, ncol, mchar, scale, ofmt
c & )

      subroutine prepare_oof_output(inputfile,inputlun,
     & outfmt1, oof_flag, irow, 
     & vmin, vmax, pindex, 
     & nflag, eflag, kflag, flag,
     & nrow, ncol, mchar, scale, ofmt, headout,naux)
     

      implicit none
      include "../ggg_int_params.f"
      include "params.f"

      integer*4
     & ncoml,ncol,kcol,icol,j,lg,
     & lnbc,nrow,li,k,krow_qc,irow,lof0,lof1,
     & eflag,wrow_flag,wsp_flag,naux,
     & mchar,lh,le,klat,klong,kzobs,ncol_written,
     & jj,nflag,ncoml_qc,ncol_qc,nrow_qc,wcol_flag
      integer*4 flag(mrow_qc),pindex(mcol),kflag(mrow_qc),
     & oof_flag(mrow)
      character
     & dl*1,
     & gggdir*(mpath),
c    & version*62,
     & headarr(mcol)*20,
     & parname(mrow_qc)*20,
     & inputfile*80,
     & outfmt0*512,
     & outfmt1*512,
     & ofmt(mrow_qc)*4,
     & temp*20,
     & fmt(mrow_qc)*4,
     & unit(mrow_qc)*6,
     & headout*8000,
     & ssss*800
      real*4 scale(mrow_qc),
     & vmin(mrow_qc),vmax(mrow_qc)
      integer mrow_qc_return, mcol_return
      integer inputlun
      integer i

      irow=1
      mrow_qc_return = mrow_qc
      mcol_return = mcol
      eflag=0
      
      do i = 1, mrow_qc
          kflag(i) = 0
          vmin(i) = 0
          vmax(i) = 0
      end do
      
c      version=
c     & ' write_official_output_file   Version 1.2.2   2009-11-07   GCT'
c      write(*,*) version

      mchar=0
      klat=0
      klong=0
      kzobs=0
      wrow_flag=1
      wsp_flag=1

      call get_ggg_environment(gggdir, dl)
      lg=lnbc(gggdir)

      li=lnbc(inputfile)
c     if(inputfile(li-3:li).ne.'.aia'.or.
c    & inputfile(li-3:li).ne.'.gaa') write(*,*)
c    & 'Warning: input file is not of expected type (.aia or .gaa)'

c  Open the Quality Control (QC) file and read in the information,
c  to find out how many columns there are going to be in the .oof output files
c     write(*,*)gggdir(:lg)//'/tccon/'//inputfile(1:2)//'_qc.dat'
      open(lun_qc,file=
     & gggdir(:lg)//'tccon'//dl//inputfile(1:2)//'_qc.dat',
     & status='old')
c     write(*,*)gggdir(:lg)//'/tccon/'//inputfile(1:2)//'_qc.dat'
      read(lun_qc,*)ncoml_qc,ncol_qc,nrow_qc
      if(nrow_qc.gt.mrow_qc) stop 'nrow_qc > mrow_qc'
      read(lun_qc,*)wrow_flag   ! 0/1   whether to skip/write the out-of-range data records
      read(lun_qc,*)wsp_flag    ! 0/1   whether to skip/write the spectrum names (if present in .aia file)
      do k=4,ncoml_qc
         read(lun_qc,'(a)') header
      end do
      ncol_written=1   ! Always include "flag" in column #1
      wrow_flag=1      ! Always include out-of-range records in the eof
      wsp_flag=1       ! Always include the spectrum names in the eof
      do krow_qc=1,mrow_qc
         read(lun_qc,*,end=77) temp,wcol_flag
         ncol_written=ncol_written+wcol_flag
      end do
77    close(lun_qc)
      if(krow_qc-1.ne.nrow_qc)
     &  stop 'misreading xx_qc.dat file: krow_qc > nrow_qc'

c  Read input file and start writing output files
c       write(*,*) inputfile
      open(inputlun,file=inputfile, status='old')
      read(inputlun,'(i2,i4,i7,i4)') ncoml,ncol,nrow,naux
c     write(*,*)'naux_prepare',naux
      if(ncol.gt.mcol) stop 'increase mcol'
      open(lun_qc,file=
     & gggdir(:lg)//'tccon'//dl//inputfile(1:2)//'_qc.dat',
     & status='old')
      read(lun_qc,*)ncoml_qc,ncol_qc,nrow_qc
      do j=2,ncoml-1
         read(inputlun,'(a)') header
      end do
      read(inputlun,'(a)') header ! column headers
      call substr(header,headarr,mcol,kcol)
      if(kcol.ne.ncol) stop 'ncol/kcol mismatch'
      if (index(header,'spectrum') .gt. 0) mchar=1

      do k=2,ncoml_qc
         read(lun_qc,'(a)') header
         if(k.eq.ncoml_qc) header=' #'//header
      end do

      do krow_qc=1,nrow_qc
         read(lun_qc,'(a)') ssss
         read(ssss,*) parname(krow_qc),flag(krow_qc),scale(krow_qc),
     &   fmt(krow_qc),unit(krow_qc),
     &   vmin(krow_qc),vmax(krow_qc)
c        write(*,*)'unit: ',unit(krow_qc)(2:lnbc(unit(krow_qc))-1)
         if (unit(krow_qc)(1:1).eq.'(') then
             unit(krow_qc)=unit(krow_qc)(1:lnbc(unit(krow_qc))-1)
             unit(krow_qc)(1:1)='_'
         endif
      end do
      close(lun_qc)

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

c         if(pindex(icol).eq.0) write(*,*)
c     &   ' Parameter missing from QC file: '//headarr(icol)

         if(pindex(icol).eq.0) then
           if(index(headarr(icol),'spectrum').eq.0) then
           write(*,*) ' Parameter missing from QC file: '//headarr(icol)
           endif
         endif

         if(headarr(icol).eq.'lat') klat=icol
         if(headarr(icol).eq.'long') klong=icol
         if(headarr(icol).eq.'zobs') kzobs=icol
      end do

      headout=' flag  spectrum'
      outfmt0='(i3'
      outfmt1='(i3,1x,a57'
      lof0=0
      lof1=0
      jj=0
      do icol=1+mchar,ncol
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
            else
               headout=headout(:lnbc(headout))//'  '//
     &         headarr(icol)(:le-1)//
     &         unit(krow_qc)(:lnbc(unit(krow_qc)))//'_error'
            endif
         endif
      end do
      outfmt0(lof0+9:lof0+9)=')'
      outfmt1(lof1+9:lof1+9)=')'

c  Read each day of data into memory and multiply XGas values by the
c  appropriate correction factors.
      nflag=0

CCC LOOP GOES HERE.


CCC

      
c      if(irow-1.ne.nrow) stop 'nrow mismatch'
c      write(*,*)nrow,' data records, of which',nflag,' flagged as bad'
c      write(*,*)
c      write(*,*)' Listed below are the fields with error flags exceeding
c     & the allowed range and the number of such occurrences'
c      write(*,*)
c      write(*,*)' #   Parameter            N_flag     %'
c      do k=1,nrow_qc
c          if(kflag(k).gt.0) write(*,'(i3,3x,a,i6,f8.1)') k,
c     &  parname(k), kflag(k), 100*float(kflag(k))/nrow
c      end do
cc     stop
      end

