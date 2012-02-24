c  Writes out a version of the oof file but with the kitchen sink. 

      implicit none
      include "../ggg_int_params.f"
      include "params.f"
      integer*4 lun_tav,lun_vav,lun_ada,lun_out,lun_outc,lun_aia,
     & li,lnbc,ncoml,ncolt,ncolv,ncola,nrow,naux,j,
     & lspace,jj,mhead,ncomlc,ncolc,ngas,lun_cor,flag,lun_asw,
     & k,maux,mwin,nwin,NN2,NN3,NN4,NN5,NN0,NN6,NN7,gaa_naux,
     & col_ncol, mcsv,kcol,gh
      integer i
      parameter (lun_tav=12)      ! input file (.tav)
      parameter (lun_vav=13)      ! input file (.vav)
      parameter (lun_cor=14)      ! input file (corrections.dat)
      parameter (lun_ada=15)      ! input file (.ada)
c      parameter (lun_gaa=31)      ! input file (.gaa)
      parameter (lun_aia=32)      ! input file (.aia)
      parameter (lun_asw=17)      ! input file (.asw)
      parameter (lun_out=16)      ! output file (.eof)
      parameter (lun_outc=18)     ! output file (.eof.csv)
      parameter (maux=25)
      parameter (mwin=20)
      parameter (mhead=110)
      parameter (mcsv=99999)
      character we_version*62,tavfile*80,vavfile*80,adafile*90,
     & aiafile*90,outfile*80,specname*(nchar),tavheader*550,
     & vavheader*550,outfilec*80,
     & adaheader*550,gaafile*90,aiaheader*550,
     & inputfmt*1024,adafmt*1024,
     & outfmt*1024,hdr(mhead)*(24),headaux(mhead)*(24),
     & corfmt*1024,ssss*16584,sarr(mrow)*57,csvarr(mrow)*57,
     & headtav(mhead)*(24),headvav(mhead)*(24),headada(mhead)*(24),
     & gasname(mwin)*20,gggdir*(mpath),corr_head*300,cc*57,dl*1,
     & ooffmt*512,oof_header*8000,oof_out*800,aswfile*80,
     & col_header*20000,col_out*10000,headercsv*20000,
     & ghost_head(mhead)*(24)
      real*4 yauxt(maux),yobst(mwin),yerrt(mwin),yauxv(maux),
     & yobsv(mwin),yerrv(mwin),
     & yauxi(maux),yobsi(mwin),yerri(mwin),yauxai(maux),yobsai(mwin),
     & yerrai(mwin),yobs_ghost(maux),yobsga(mwin)

      real*8 year,adcf(mwin),adcf_err(mwin),aicf(mwin),aicf_err(mwin)
      integer*4 oof_flag(mrow)
      integer*4 ncol, mchar, eflag, irow, nflag
      integer*4 kflag(mrow_qc), oflag(mrow_qc), 
     & pindex(mrow_qc)
      real*4 vmin(mrow_qc), vmax(mrow_qc)
      real*4 scale(mrow_qc)
      character ofmt(mrow_qc)*4
      integer*4 inputlun
      integer*4 percent_complete, percent_complete_target
     
      integer*4 luns(mluns), noc

      integer grl_array_size
      parameter (grl_array_size=10)
      character specname_grl_array(grl_array_size)*57
      integer iyr_array(grl_array_size), doy_array(grl_array_size)
      integer delta_t_array(grl_array_size)
      integer zpdtim_array(grl_array_size)
      integer grl_array_counter


      call getlun(inputlun)
      we_version=
     &' write_eof  Version 0.3.2   2012-01-17   DW/CL'
      write(*,*) we_version

      write(*,'(a)')
     & 'Enter name of .tav file'
      read(*,'(a)')tavfile
      li=lnbc(tavfile)
      if(tavfile(li-3:li).ne.'.tav') write(*,*)
     &  ' Warning: input file is not of expected type (.tav) '

c Set the names of the input files to be read, and the output files to be written
      vavfile=tavfile(:li-3)//'vav'
      aiafile=vavfile(:li)//'.ada.aia'
      gaafile=vavfile(:li)//'.ada.aia.gaa'

      adafile=vavfile(:li)//'.ada'
      aswfile=vavfile(:li-3)//'asw'

      outfile=vavfile(:li)//'.ada.aia.gaa.eof'
      outfilec=vavfile(:li)//'.ada.aia.gaa.eof.csv'
     
c Open the other files
      open(lun_vav,file=vavfile,status='old')
      open(lun_tav,file=tavfile,status='old')
      open(lun_aia,file=aiafile,status='old')
c     open(lun_gaa,file=gaafile,status='old')
      open(lun_ada,file=adafile,status='old')
      open(lun_out,file=outfile,status='unknown')
      open(lun_outc,file=outfilec,status='unknown')

CC prepare oof

       call prepare_oof_output(gaafile,inputlun,
     & ooffmt, oof_flag, irow, 
     & vmin, vmax, pindex, 
     & nflag, eflag, kflag, oflag,
     & nrow, ncol, mchar, scale, ofmt, oof_header,gaa_naux)
CC
       call prepare_collate_all(col_header, luns, col_ncol, noc,
     & grl_array_size, specname_grl_array, iyr_array, doy_array,
     & delta_t_array, zpdtim_array, grl_array_counter,gfit_version,
     & gsetup_version, atmsum, gctsum, fciasum, sciasum, solarsum)

c Read and parse the corrections.dat file
      call get_ggg_environment(gggdir, dl)
      open(lun_cor,file=gggdir(:lnbc(gggdir))//'tccon'//dl
     & //'corrections.dat', status='old')
      read(lun_cor,*)ncomlc,ncolc
      do k=2,ncomlc
         read(lun_cor,*)
      end do
      if(ncolc.eq.3) then
         do k=1,mwin
            read(lun_cor,*,end=88) gasname(k),adcf(k),aicf(k)
         end do
      elseif(ncolc.eq.5) then
         do k=1,mwin
            read(lun_cor,*,end=88) gasname(k),
     &      adcf(k),adcf_err(k),aicf(k),aicf_err(k)
         end do
      else
         write(*,*)'ncolc=',ncolc
         stop 'Unrecognized NCOLC value'
      endif
88    ngas=k-1
      if(ncolc.eq.3) then
         write(corr_head,*)(gasname(k)(:lnbc(gasname(k)))//'_adcf ',
     & gasname(k)(:lnbc(gasname(k)))//'_aicf ',k=1,ngas)
      elseif(ncolc.eq.5) then
         write(corr_head,*)(gasname(k)(:lnbc(gasname(k)))//'_adcf ',
     & gasname(k)(:lnbc(gasname(k)))//'_adcf_err ',
     & gasname(k)(:lnbc(gasname(k)))//'_aicf ',
     & gasname(k)(:lnbc(gasname(k)))//'_aicf_err ',
     & k=1,ngas)
      endif

c Parse the tav file header
      read(lun_tav,'(i2,i4,i7,i4)') ncoml,ncolt,nrow,naux
      do j=2,ncoml-1
         read(lun_tav,'(a)') tavheader
      end do
      read(lun_tav,'(a)') tavheader
c     write(*,*)'tavheader='//tavheader
      if(index(tavheader,'spectrum').le.0) then
         write(*,*) ' Error: write_eof requires that spectrum names be
     &  present in input files'
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
          write(ghost_head(gh),'(a)')'x'//hdr(k)(:lnbc(hdr(k)))//
     &    '_ghost_corr'
          endif
      enddo
      
c Read the .vav, .vav.ada files for inclusion
      read(lun_vav,'(i2,i4,i7,i4)') ncoml,ncolv,nrow,naux
      do j=2,ncoml-1
         read(lun_vav,'(a)') vavheader
      end do
      read(lun_vav,'(a)') vavheader
      if(index(vavheader,'spectrum').le.0) then
         write(*,*) ' Error: write_eof requires that spectrum names be
     &  present in input files'
         stop ': spectrum names absent from .tav file'
      endif

      read(lun_ada,'(i2,i4,i7,i4)') ncoml,ncola,nrow,naux
      do j=2,ncoml-1
         read(lun_ada,'(a)') adaheader
      end do
      read(lun_ada,'(a)') adaheader
      
      read(lun_aia,'(i2,i4,i7,i4)') ncoml,ncola,nrow,naux
      do j=2,ncoml-1
         read(lun_aia,'(a)') aiaheader
      end do
      read(lun_aia,'(a)') aiaheader

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
     & (ghost_head(k)(:lnbc(ghost_head(k))+1),k=1,gh),
     & col_header(:lnbc(col_header))
c Create the output csv header
      headercsv=header
      
      call substr(headercsv(:lnbc(headercsv)),csvarr,mrow,kcol)
      headercsv=csvarr(1)
      do k=2,kcol
         cc=csvarr(k)
         headercsv=headercsv(:lnbc(headercsv))//','//cc(:lnbc(cc))
      enddo

c Write the top part of the output files
      write(lun_out,'(i2,x,i5,x,i7)') 3,kcol,nrow
      write(lun_out,'(a)') we_version
      write(lun_out,'(a)')header(:lnbc(header))

      write(lun_outc,'(i2,x,i5,x,i7)') 3,kcol,nrow
      write(lun_outc,'(a)') we_version
      write(lun_outc,'(a)')headercsv(:lnbc(headercsv))

c Set the formats of the files
      inputfmt='(a57,1x,f13.8,NNf13.5,800(1pe12.4))'
      write(inputfmt(15:16),'(i2.2)')naux-1

      adafmt='(a57,1x,f13.8,NNf13.5,200(1pe12.4))'
      write(adafmt(15:16),'(i2.2)') naux-1

      corfmt='(NNf8.4)'
      write(corfmt(2:3),'(i2.2)') ngas*2

      nwin=(ncolt-naux)/2

c Read in each line of the files, and write the string to the output file

C prepare oof and col files for reading

  

   
C end prepare files

      percent_complete_target=0.0
      do j=1,nrow
          percent_complete = j*100/nrow
          if(percent_complete .gt. percent_complete_target) then
              write(*,'(a,x,i4,x,a)') "Completed", 
     &                          percent_complete_target,"%"
              percent_complete_target=percent_complete_target+1+100/nrow
           endif
          oof_out(:)=""
         call read_oneline_oof(gaafile,inputlun, oof_out,
     & ooffmt, oof_flag, irow,
     & vmin, vmax, pindex,
     & nflag, eflag, kflag, oflag,
     & nrow, ncol, mchar, scale, ofmt, gaa_naux,yobsga
     & )

         call read_oneline_col(col_out, luns, col_ncol, noc,
     & grl_array_size, specname_grl_array, iyr_array, doy_array,
     & delta_t_array, zpdtim_array, grl_array_counter, gfit_version,
     & gsetup_version,atmsum,gctsum,fciasum,sciasum,solarsum)

         outfmt='(aNNX,'
     &    //'NN(1pe12.4),'
     &    //'NN(1pe12.4),'
     &    //'NN(1pe12.4),0p,NN(f8.4),1x,NN(e10.3)'
     &    //'aNNXX)'!,1x,a32,1x,a32,1x,a32,1x,a32,1x,a32)'
         NN0 = index(outfmt,'NN')
         write(outfmt(NN0:NN0+2),'(i3.3)')lnbc(oof_out)
         NN2 = index(outfmt(NN0+2:),'NN')+NN0+1
         write(outfmt(NN2:NN2+1),'(i2.2)')nwin*2
         NN3 = index(outfmt(NN2+2:),'NN')+NN2+1
         write(outfmt(NN3:NN3+1),'(i2.2)')nwin*2
         NN4 = index(outfmt(NN3+2:),'NN')+NN3+1
         write(outfmt(NN4:NN4+1),'(i2.2)')nwin*2
         NN5 = index(outfmt(NN4+2:),'NN')+NN4+1

         if(ncolc.eq.3) then
             write(outfmt(NN5:NN5+1),'(i2.2)')ngas*2
         elseif(ncolc.eq.5) then
             write(outfmt(NN5:NN5+1),'(i2.2)')ngas*4
         endif

         NN6 = index(outfmt(NN5+2:),'NN')+NN5+1
         write(outfmt(NN6:NN6+1),'(i2.2)')nwin
         NN7 = index(outfmt(NN6+2:),'NN')+NN6+1
         write(outfmt(NN7:NN7+3),'(i4.4)')lnbc(col_out) 
c        write(*,*)'outfmt=',outfmt

c The above assumes that all the oof and col file lines are the same length!

c read one line of the oof file

c read one line of the col file(s)

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

         read(lun_aia,adafmt) 
     &   specname,year,
     &   (yauxai(k),k=1,naux-2),
     &   (yobsai(k),yerrai(k),k=1,nwin)

c        read(lun_gaa,adafmt) 
c    &   specname,year,
c    &   (yauxga(k),k=1,naux-2),
c    &   (yobsga(k),yerrga(k),k=1,nwin)

C differences between the aia file and the gaa file are the additional ghost corrections
      do k=1,nwin
         yobs_ghost(k) = yobsga(k) - yobsai(k)
      enddo
c      write(*,*)'yobs_ghost=',yobs_ghost(:nwin)

c Write the output to a string for parsing into csv later

      if(ncolc.eq.3) then
         write(ssss,outfmt)oof_out(:lnbc(oof_out)),
     &   (yobst(k),yerrt(k),k=1,nwin),
     &   (yobsv(k),yerrv(k),k=1,nwin),
     &   (yobsi(k),yerri(k),k=1,nwin),
     &   (adcf(k),aicf(k),k=1,ngas),
     &   (yobs_ghost(k),k=1,nwin),
     &   col_out(:lnbc(col_out))
      elseif(ncolc.eq.5) then
         write(ssss,outfmt)oof_out(:lnbc(oof_out)),
     &   (yobst(k),yerrt(k),k=1,nwin),
     &   (yobsv(k),yerrv(k),k=1,nwin),
     &   (yobsi(k),yerri(k),k=1,nwin),
     &   (adcf(k),adcf_err(k),aicf(k),aicf_err(k),k=1,ngas),
     &   (yobs_ghost(k),k=1,nwin),
     &   col_out(:lnbc(col_out))
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
!            ssss=ssss(:lnbc(ssss))//','//cc(:lnbc(cc))
         enddo
         write(lun_outc,'(a)')ssss(:lnbc(ssss))
         ssss=''
      enddo ! do j=1,nrow
      write(*,'(a,x,i4,x,a)') "Completed  100 %"
c Close files
c      close(lun_oof)
c      close(lun_col)

      call freelun(inputlun)
      do j=1, mluns
         if (luns(j) .gt. 0) call freelun(luns(j))
      end do
      
      close(lun_tav)
      close(lun_vav)
      close(lun_ada)
      close(lun_aia)
c     close(lun_gaa)
      close(lun_out)
      close(lun_outc)
      end
