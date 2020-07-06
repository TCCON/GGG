c  Writes out a version of the oof file but with the kitchen sink. 

      implicit none
      include "../gfit/ggg_int_params.f"
      include "params.f"
      include "../comn/postproc_params.f"
      integer*4 lun_tav,lun_vav,lun_ada,lun_out,lun_outc,lun_aia,
     & li,lnbc,ncoml,ncolt,ncolv,ncola,nrow,naux,j,lnit,lun_esf,
     & lspace,jj,mhead,ncomlc,ncolc,lun_cor,flag,
     & k,maux,mwin,nwin,NN2,NN3,NN4,NN5,NN0,NN6,NN7,gaa_naux,
     & col_ncol,kcol,gh,nlhead_esf,ncol_esf,nrow_qc,ls,nrow1_qc,
     & nheaders,ncolumns,lunr_lse,ndum,mnum,nnum,lse_count,
     & ncoml_qc,krow_qc,ncol1_qc,flag_qc(mrow_qc),
     & ada_ncorr, ada_has_err, aia_ncorr, aia_has_err
      integer i
      parameter (lun_tav=12)      ! input file (.tav)
      parameter (lun_vav=13)      ! input file (.vav)
      parameter (lun_cor=14)      ! input file (corrections.dat)
      parameter (lun_ada=15)      ! input file (.ada)
c      parameter (lun_gaa=31)      ! input file (.gaa)
      parameter (lunr_lse=17)     ! input file (.lse)
      parameter (lun_aia=32)      ! input file (.aia)
      parameter (lun_out=16)      ! output file (.eof)
      parameter (lun_outc=18)     ! output file (.eof.csv)
      parameter (lun_esf=29)     ! output file (.eof.csv)
      parameter (mnum=5)
      parameter (maux=25)
      parameter (mwin=100)
      parameter (mhead=110)
      character we_version*62,tavfile*80,vavfile*80,adafile*90,
     & aiafile*90,outfile*80,specname*(nchar),tavheader*800,
     & vavheader*800,outfilec*80,lsefile*(mfilepath),lseformat*100,
     & adaheader*800,
c     aiaheader*800,gaafile*90,ghost_head(mhead)*(24),yobs_ghost(maux),
     & esffile*90,
     & inputfmt*1024,adafmt*1024,specname_lse*(nchar),
     & outfmt*1024,hdr(mhead)*(24),headaux(mhead)*(24),
     & ssss*16584,sarr(mrow)*(nchar),csvarr(mrow)*(nchar),
     & headtav(mhead)*(24),headvav(mhead)*(24),headada(mhead)*(24),
     & headlse(mhead)*(24),lseheader*256,specname_lse_was*(nchar),
     & ada_gasname(mwin)*20,aia_gasname(mwin)*20,gggdir*(mpath),
     & ada_corr_head*900,aia_corr_head*900,corr_head*900,cc*(nchar),
     & dl*1,ooffmt*512,oof_header*8000,oof_out*800,aswfile*80,
     & col_header*36000,col_out*36000,headercsv*35000,qc_header*1800,
     & specfmt*3,esf_header*8000,specname_oof*(nchar),ssssc*800,
     & parname(mrow_qc)*20

      real*4 yauxt(maux),yobst(mwin),yerrt(mwin),yauxv(maux),
     & yobsv(mwin),yerrv(mwin),
c     & yauxai(maux),yobsai(mwin),yerrai(mwin),
     & yauxi(maux),yobsi(mwin),yerri(mwin),yobsga(mwin)

      real*8 year,adcf(mwin),adcf_err(mwin),aicf(mwin),aicf_err(mwin),
     & dum4,lse,lsu,lsf,dip,mvd
      integer*4 oof_flag(mrow),ktg(mcol),lst
      integer*4 ncol, mchar, eflag, irow, nflag,dum2,dum3
      integer*4 kflag(mrow_qc), oflag(mrow_qc), 
     & pindex(mrow_qc)
      real*4 vmin(mrow_qc),vmax(mrow_qc),scale_qc(mrow_qc)
      real*4 vmin_qc(mrow_qc),vmax_qc(mrow_qc),rsc_qc(mrow_qc)
      real*4 scale(mrow_qc), yesf(mcol),rsc(mrow_qc)
      character ofmt(mrow_qc)*4,fmt_qc(mrow_qc)*4,unit_qc(mrow_qc)*6
      integer*4 idum
      integer*4 percent_complete, percent_complete_target
     
      integer*4 luns(mluns)
      save yesf

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)
      idum=maddln   ! Avoid compiler warning (unused parameter)
      idum=mcharhead! Avoid compiler warning (unused parameter)
      idum=lun_qc   ! Avoid compiler warning (unused parameter)
      idum=lun_rlg  ! Avoid compiler warning (unused parameter)
      idum=lun_mul  ! Avoid compiler warning (unused parameter)

c      call getlun(inputlun)
      we_version=
     & ' WRITE_EOF        Version 0.49a   2020-03-12   GCT/DW/JLL '

      write(*,*) we_version

      if (iargc() == 0) then
         write(*,'(a)') 'Enter name of .tav file'
         read(*,'(a)')tavfile
      elseif (iargc() == 1) then
         call getarg(1, tavfile)
      else
         stop 'Usage: $gggpath/bin/write_eof tavfilename'
      endif
      li=lnbc(tavfile)
      if(tavfile(li-3:li).ne.'.tav') write(*,*)
     &  ' Warning: input file is not of expected type (.tav) '

c Set the names of the input files to be read, and the output files to be written
      vavfile=tavfile(:li-3)//'vav'
      aiafile=vavfile(:li)//'.ada.aia'
c      gaafile=vavfile(:li)//'.ada.aia.gaa'
      esffile=vavfile(:li)//'.ada.aia.daily_error.out'

      adafile=vavfile(:li)//'.ada'
      aswfile=vavfile(:li-3)//'asw'

      outfile=vavfile(:li)//'.ada.aia.eof'
      outfilec=vavfile(:li)//'.ada.aia.eof.csv'
     
      specname_lse = 'first_spectrum.000'
c Open the other files
      open(lun_vav,file=vavfile,status='old')
      open(lun_tav,file=tavfile,status='old')
c     open(lun_aia,file=aiafile,status='old')
c     open(lun_gaa,file=gaafile,status='old')
      open(lun_ada,file=adafile,status='old')
      open(lun_out,file=outfile,status='unknown')
      open(lun_outc,file=outfilec,status='unknown')

      open(lun_esf,file=esffile,status='old')
      read(lun_esf,*) nlhead_esf,ncol_esf
      do i=2,nlhead_esf
         read(lun_esf,'(a)') esf_header
      end do
      do j=1,mcol
         yesf(j)=0.0
      end do

C prepare oof and col files for reading

      call prepare_oof_output(aiafile,lun_aia,
     & ooffmt, oof_flag, irow,nrow_qc,
     & ada_ncorr, ada_has_err, ada_gasname, adcf, adcf_err, 
     & aia_ncorr, aia_has_err, aia_gasname, aicf, aicf_err, 
     & vmin, vmax, pindex, 
     & eflag, kflag, oflag,
     & nrow, ncol, mchar, scale,rsc, ofmt, oof_header,gaa_naux)
c
c verify that the adcf and aicf both have or do not have error -
c not implemented if one does and the other doesn't
c
      if(ada_has_err .ne. aia_has_err)
     & stop 'Only one of adcf/aicf has an error, not implemented'
CC
      call prepare_collate_all(col_header, luns, col_ncol,ktg, 
     & gfit_version,gsetup_version, tllsum, solarsum,lnit,csformat)

c Reorganize the corrections read from the aia file
      call get_ggg_environment(gggdir, dl)
      
      if(ada_has_err .eq. 1) then
         write(ada_corr_head,*)
     & (ada_gasname(k)(:lnbc(ada_gasname(k)))//'_adcf ', 
     &  ada_gasname(k)(:lnbc(ada_gasname(k)))//'_adcf_err ',
     &  k=1, ada_ncorr)
      else
         write(ada_corr_head,*)
     & (ada_gasname(k)(:lnbc(ada_gasname(k)))//'_adcf ', 
     &  k=1, ada_ncorr)
      endif

      if(aia_has_err .eq. 1) then
        write(aia_corr_head,*)
     & (aia_gasname(k)(:lnbc(aia_gasname(k)))//'_aicf ',
     &  aia_gasname(k)(:lnbc(aia_gasname(k)))//'_aicf_err ',
     &  k=1, aia_ncorr)
      else
        write(aia_corr_head,*)
     & (aia_gasname(k)(:lnbc(aia_gasname(k)))//'_aicf ',
     &  k=1, aia_ncorr)
      endif

      write(corr_head,'(a,1x,a)')
     & ada_corr_head(:lnbc(ada_corr_head)),
     & aia_corr_head(:lnbc(aia_corr_head))

c Parse the tav file header
      read(lun_tav,countfmt) ncoml,ncolt,nrow,naux
      do j=2,ncoml-1
         read(lun_tav,'(a)') tavheader
      end do
      read(lun_tav,'(a)') tavheader
c     write(*,*)'tavheader='//tavheader
      if(index(tavheader,'spectrum').le.0) then
         write(*,*) ' Error: write_eof requires that spectrum names be
     &   present in input files'
         stop ': spectrum names absent from .tav file'
      endif
      k=1
      flag=0
      jj=1
      do j=1,99999
         if(tavheader(k:k).eq.' ') then
            k=k+1
         else
            lspace=index(tavheader(k:),' ')+k
            if(index(tavheader(k:),' ')+k.ge.lnbc(tavheader)) then
               hdr(jj)=tavheader(k:)
               flag=1
            else
               hdr(jj)=tavheader(k:lspace)
            endif
            jj=jj+1
            k=lspace+1
         endif
         if(flag.eq.1) then
            goto 77
         endif
      enddo ! j=1,99999
77    continue

c Write the headers for the output files
      do k=1,naux
         write(headaux(k),'(a)')hdr(k)(:lnbc(hdr(k)))
      enddo
      gh=0
      do k=naux+1,jj-1
         j=k-(naux+1)+1
         write(headtav(j),'(a)')'vsf_'//hdr(k)(:lnbc(hdr(k)))
         write(headvav(j),'(a)')'column_'//hdr(k)(:lnbc(hdr(k)))
         write(headada(j),'(a)')'ada_x'//hdr(k)(:lnbc(hdr(k)))
         if(index(hdr(k),'_error').gt.0) then
         else  ! we only have corrections for the values, not the errors
            gh = gh+1
c            write(ghost_head(gh),'(a)')'x'//hdr(k)(:lnbc(hdr(k)))//
c     &      '_ghost_corr'
         endif
      enddo
      
c Read the .vav, .vav.ada files for inclusion
      read(lun_vav,countfmt) ncoml,ncolv,nrow,naux
      do j=2,ncoml-1
         read(lun_vav,'(a)') vavheader
      end do
      read(lun_vav,'(a)') vavheader
      if(index(vavheader,'spectrum').le.0) then
         write(*,*) ' Error: write_eof requires that spectrum names be
     &   present in input files'
         stop ': spectrum names absent from .tav file'
      endif

      read(lun_ada,countfmt) ncoml,ncola,nrow,naux
      do j=2,ncoml-1
         read(lun_ada,'(a)') adaheader
      end do
      read(lun_ada,'(a)') adaheader
      
      lsefile=gggdir(:lnbc(gggdir))//'lse'//dl
     & //'gnd'//dl//tavfile(:li-3)//'lse'
      open(lunr_lse,file=lsefile,status='old')
      read(lunr_lse,'(i2,i4)') nheaders,ncolumns
      do j=2,nheaders-1
         read(lunr_lse,'(a)') lseformat
c         write(*,'(a)') ' lse format:'//lseformat
      enddo
      read(lunr_lse,'(a)') lseheader
c      write(*,*) ' lse header:'//lseheader
      call substr(lseheader,headlse,ncolumns,ndum)
c     write(*,*)'headlse=',(headlse(k)(:lnbc(headlse(k))),k=1,ncolumns)

c     read(lun_aia,'(i2,i4,i7,i4)') ncoml,ncola,nrow,naux
c     do j=2,ncoml-1
c        read(lun_aia,'(a)') aiaheader
c     end do
c     read(lun_aia,'(a)') aiaheader

c     write(*,*)'lun_gaa=',lun_gaa
c     read(lun_gaa,'(a)') gaaheader
c     write(*,*)'gaaheader',gaaheader
c     read(lun_gaa,'(i2,i4,i7,i4)') ncoml,ncola,nrow,naux
c     do j=2,ncoml-1
c        read(lun_gaa,'(a)') gaaheader
c     end do
c     read(lun_gaa,'(a)') gaaheader

c Create the output header from the various file headers
      write(header,*)oof_header(:lnbc(oof_header)),' ',
     & (headtav(k)(:lnbc(headtav(k))+1),k=1,jj-naux-1),
     & (headvav(k)(:lnbc(headvav(k))+1),k=1,jj-naux-1),
     & (headada(k)(:lnbc(headada(k))+1),k=1,jj-naux-1),
     & corr_head(:lnbc(corr_head)),' ',
     & (headlse(k)(:lnbc(headlse(k))+1),k=5,10),
c    & (ghost_head(k)(:lnbc(ghost_head(k))+1),k=1,gh),
     & col_header(:lnbc(col_header))
c Create the output csv header
      headercsv=header
      
      call substr(headercsv(:lnbc(headercsv)),csvarr,mrow,kcol)
      headercsv=csvarr(1)
      do k=2,kcol
         cc=csvarr(k)
         headercsv=headercsv(:lnbc(headercsv))//','//cc(:lnbc(cc))
      enddo

c Get the number of header lines from the qc.dat file

c      write(*,*)gggdir(:lnbc(gggdir))//
c     & '/tccon/'//tavfile(1:2)//'_qc.dat'
      open(lun_qc,file=
     & gggdir(:lnbc(gggdir))//'tccon'//dl//
     & tavfile(1:2)//'_qc.dat',
     & status='old')
      read(lun_qc,*)ncoml_qc,ncol1_qc,nrow1_qc

c Write the top part of the output files
      write(lun_out,'(i2,x,i5,x,i7)') 4+nrow1_qc,kcol,nrow
      write(lun_out,'(a)') we_version

      write(ssss,*) 4+nrow1_qc,kcol,nrow
      call substr(ssss,sarr,mnum,nnum)
      ssss=sarr(1)
      do k=2,nnum
         cc=sarr(k)(:20)
         ssss=ssss(:lnbc(ssss))//','//cc(:lnbc(cc))
      end do
      write(lun_outc,'(a)') ssss(:lnbc(ssss))
c      write(lun_outc,'(i2,x,i5,x,i7)') 3,kcol,nrow
      write(lun_outc,'(a)') we_version

c add the qc.dat file output
      do k=2,ncoml_qc
         read(lun_qc,'(a)') qc_header
         if(k.eq.ncoml_qc) then
             qc_header=' #'//qc_header(:1800-2)
             write(lun_out,'(a)') qc_header(:lnbc(qc_header))
             write(lun_outc,'(a)') qc_header(:lnbc(qc_header))
         end if
c         write(lun_out,'(a)') qc_header(:lnbc(qc_header))
c         write(lun_outc,'(a)') qc_header(:lnbc(qc_header))
      end do

      do krow_qc=1,nrow1_qc
         read(lun_qc,'(a)') ssssc
         write(lun_out,'(i3,a)') krow_qc,ssssc(:lnbc(ssssc))
         write(lun_outc,'(i3,a)') krow_qc,ssssc(:lnbc(ssssc))
         read(ssssc,*) parname(krow_qc),flag_qc(krow_qc),
     &   scale_qc(krow_qc),fmt_qc(krow_qc),unit_qc(krow_qc),
     &   vmin_qc(krow_qc),vmax_qc(krow_qc)
         rsc_qc(krow_qc)=scale_qc(krow_qc)
      end do
      close(lun_qc)

      write(lun_out,'(a)')header(:lnbc(header))
      write(lun_outc,'(a)')headercsv(:lnbc(headercsv))

c Set the formats of the files
      write(specfmt,'(a,i2)')'a',nchar
      inputfmt='('//specfmt//',1x,f13.8,NNf13.5,800(1pe12.4))'
      write(inputfmt(15:16),'(i2.2)')naux-1
      adafmt='('//specfmt//',1x,f13.8,NNf13.5,200(1pe12.4))'
      write(adafmt(15:16),'(i2.2)') naux-1

      nwin=(ncolt-naux)/2

C end prepare files

c Read in each line of the files, and write the string to the output file
      percent_complete_target=0.0
      do j=1,nrow
         percent_complete = j*100/nrow
         if(percent_complete .gt. percent_complete_target) then
            write(*,'(a,x,i4,x,a)') "Completed", 
     &      percent_complete_target,"%"
            percent_complete_target=percent_complete_target+1+100/nrow
         endif
         oof_out(:)=""

c read one line of the oof file
         call read_oneline_oof(lun_aia, oof_out,
     & ooffmt, oof_flag, irow, lun_esf, yesf,ncol_esf,nrow_qc,
     & vmin, vmax, pindex,
     & nflag, eflag, kflag, oflag,
     & ncol, mchar, scale,rsc, ofmt, gaa_naux,yobsga,specname_oof
     & )

c read one line of the col file(s)
         call read_oneline_col(col_out, luns, col_ncol,csformat,
     & lnit,ktg,specname_oof, 
     & gfit_version,
     & gsetup_version,tllsum,solarsum)

c set the output file format
         outfmt='(aNNX,'
     &    //'NN(1pe12.4),'
     &    //'NN(1pe12.4),'
     &    //'NN(1pe12.4),0p,NN(f8.4),1x,'
     &    //lseformat(21:65)//',1x,'
     &    //'aNNXXX)'!,1x,a32,1x,a32,1x,a32,1x,a32,1x,a32)'
         NN0 = index(outfmt,'NN')
         write(outfmt(NN0:NN0+2),'(i3.3)')lnbc(oof_out)
         NN2 = index(outfmt(NN0+2:),'NN')+NN0+1
         write(outfmt(NN2:NN2+1),'(i2.2)')nwin*2
         NN3 = index(outfmt(NN2+2:),'NN')+NN2+1
         write(outfmt(NN3:NN3+1),'(i2.2)')nwin*2
         NN4 = index(outfmt(NN3+2:),'NN')+NN3+1
         write(outfmt(NN4:NN4+1),'(i2.2)')nwin*2
         NN5 = index(outfmt(NN4+2:),'NN')+NN4+1
         write(outfmt(NN5:NN5+1),'(i2.2)')
     & (1+ada_has_err)*ada_ncorr + (1+aia_has_err)*aia_ncorr

c         NN6 = index(outfmt(NN5+2:),'NN')+NN5+1
c         write(outfmt(NN6:NN6+1),'(i2.2)')nwin
         NN6 = NN5
         NN7 = index(outfmt(NN6+2:),'NN')+NN6+1
         write(outfmt(NN7:NN7+4),'(i5.5)')lnbc(col_out)
c        write(*,*)'outfmt=',outfmt

c The above assumes that all the oof and col file lines are the same length!

         read(lun_tav,inputfmt) 
     &   specname,year,
     &   (yauxt(k),k=1,naux-2),
     &   (yobst(k),yerrt(k),k=1,nwin)
         read(lun_vav,inputfmt) 
     &   specname,year,
     &   (yauxv(k),k=1,naux-2),
     &   (yobsv(k),yerrv(k),k=1,nwin)

         read(lun_ada,adafmt) 
     &   specname,year,
     &   (yauxi(k),k=1,naux-2),
     &   (yobsi(k),yerri(k),k=1,nwin)

         specname_lse_was = specname_lse
         read(lunr_lse,lseformat) 
     &   specname_lse,dum2,dum3,dum4,
     &   lst,lse,lsu,lsf,dip,mvd
         lse_count = 0
         if(specname_lse.ne.specname) then ! specname doesn't match
            do while (specname_lse.ne.specname)
               lse_count = lse_count + 1
               read(lunr_lse,lseformat) 
     &         specname_lse,dum2,dum3,dum4,lst,lse,lsu,lsf,dip,mvd
               ls=lnbc(specname)
               if (lse_count.gt.4) then
                  write(*,*)'ERROR: Cannot find ',specname(:ls),
     &            ' within 4 lines of ',specname_lse_was(:ls),' in '//
     &            'your .lse file.'
                  write(*,*) 'There may be something wrong with '//
     &            'your .lse file.'
                  write(*,*) 'Please check that the .lse file is '//
     &            'consistent with your runlog and re-run write_eof.'
                  stop
               endif
            enddo
         endif
c        read(lun_aia,adafmt) 
c    &   specname,year,
c    &   (yauxai(k),k=1,naux-2),
c    &   (yobsai(k),yerrai(k),k=1,nwin)

c        read(lun_gaa,adafmt) 
c    &   specname,year,
c    &   (yauxga(k),k=1,naux-2),
c    &   (yobsga(k),yerrga(k),k=1,nwin)

C differences between the aia file and the gaa file are the additional ghost corrections
c      do k=1,nwin
c         yobs_ghost(k) = yobsga(k) - yobsai(k)
c      enddo
c      write(*,*)'yobs_ghost=',yobs_ghost(:nwin)

c Write the output to a string for parsing into csv later

         if(ada_has_err .eq. 1 .and. aia_has_err .eq. 1) then
            write(ssss,outfmt)oof_out(:lnbc(oof_out)),
     &      (yobst(k),yerrt(k),k=1,nwin),
     &      (yobsv(k),yerrv(k),k=1,nwin),
     &      (yobsi(k),yerri(k),k=1,nwin),
     &      (adcf(k),adcf_err(k),k=1,ada_ncorr),
     &      (aicf(k),aicf_err(k),k=1,aia_ncorr),
     &      lst,lse,lsu,lsf,dip,mvd,
     &      col_out(:lnbc(col_out))
         elseif(ada_has_err .eq. 1 .and. aia_has_err .eq. 1) then
            write(ssss,outfmt)oof_out(:lnbc(oof_out)),
     &      (yobst(k),yerrt(k),k=1,nwin),
     &      (yobsv(k),yerrv(k),k=1,nwin),
     &      (yobsi(k),yerri(k),k=1,nwin),
     &      (adcf(k),k=1,ada_ncorr),
     &      (aicf(k),k=1,aia_ncorr),
     &      lst,lse,lsu,lsf,dip,mvd,
     &      col_out(:lnbc(col_out))
         else
            stop 'Not implemented: aicf/adcf_has_err mismatch'
         endif
         write(lun_out,'(a)')ssss(:lnbc(ssss))

c Parse the string to create the comma-separated string
         call substr(ssss(:lnbc(ssss)),sarr,mrow,kcol)
         
         ssss=sarr(1)
         i=lnbc(ssss)
         i=i+1
         jj=0
         do k=2,kcol
            cc=sarr(k)
            jj=lnbc(cc)
            ssss(i:i)=','
            ssss(i+1:i+1+jj)=cc(:jj)
            i=i+jj+1
c            ssss=ssss(:lnbc(ssss))//','//cc(:lnbc(cc))
         enddo
         write(lun_outc,'(a)')ssss(:lnbc(ssss))
         ssss=''
      enddo ! do j=1,nrow
      write(*,'(a,x,i4,x,a)') "Completed  100 %"
c Close files
c      close(lun_oof)
c      close(lun_col)

c      call freelun(inputlun)
      do j=1, mluns
         if (luns(j) .gt. 0) call freelun(luns(j))
      end do
      
      close(lun_tav)
      close(lun_vav)
      close(lun_ada)
      close(lun_aia)
c      close(lun_gaa)
      close(lunr_lse)
      close(lun_out)
      close(lun_outc)
      end
