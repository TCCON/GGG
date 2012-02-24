c  write_aux.f
c  Writes out a nice version of the mav file for data users. 

      implicit none
      include "../ggg_int_params.f"

      integer*4 lun_mav,lun_out,nlev,nspeci,mvmr,lr,lnbc,lcolon,i,
     & j,ii,lspace,ii1,jj,ico2,ico,ich4,in2o,ih2o,ihdo,flag,nlhead
      parameter (lun_mav=12)      ! input file (.mav)
      parameter (lun_out=14)      ! output file (.map)
      parameter (mvmr=28000)
      real*4 z(mlev),t(mlev),p(mlev),d(mlev),vmr(mvmr)
      character version*62,mavfile*80,outfile*80,mavstring*64,
     & runlabmav*57,string*48,head*1800,hdr(mlev)*(11),string1*25,
     & mav_version*64

      version=
     &' write_aux              version 1.0.0   2012-01-10   DW'
      write(*,*) version

      write(*,'(a)')
     & 'Enter name of .mav file'
      read(*,'(a)')mavfile
c     lr=lnbc(mavfile)

c  Read the entire contents of the .mav file
      open(lun_mav,file=mavfile,status='old')
      read(lun_mav,'(a)')mav_version
      read(lun_mav,'(14x,a)')string
      lcolon=index(string,':')
      read(string(lcolon+1:),'(a)')runlabmav

      do i=1,9999
          call read_mav_head(lun_mav,nlev,nspeci,nlhead)
          call read_mav_aux(lun_mav,nlhead,nlev,nspeci,z,t,p,d,vmr,head)
c         write(*,*)runlabmav
c         write(*,*)head
c         write(*,*)mav_version
c         outfile=runlabmav(:lr-13)//'.map'
c         lr=lnbc(runlabmav)

c  Figure out where the date ends (usually has an 's' for TCCON spectrum names).
          lr=index(runlabmav,'s')
          if (lr.eq.0) lr=11
c         write(*,*)lr
          outfile=runlabmav(:lr-1)//'.map'

c  Figure out what is in the mav file and select out the appropriate columns
          ii1=1
          jj=1
          ico2=0
          ico =0
          ich4=0
          in2o=0
          ih2o=0
          ihdo=0
          flag=0
          do ii=1,99999
c             write(*,*)ii1
              if(head(ii1:ii1).eq.' ') then
                  ii1=ii1+1
              else
                  lspace=index(head(ii1:),' ')+ii1
                  if(index(head(ii1:),' ')+ii1.ge.lnbc(head)) then
                   hdr(jj)=head(ii1:)
                       flag=1
                  else
                       hdr(jj)=head(ii1:lspace)
                  endif
                  if(hdr(jj).eq.'1co2') then
                       ico2 = jj
                  elseif(hdr(jj).eq.'1co') then
                       ico  = jj
                  elseif(hdr(jj).eq.'1h2o') then
                       ih2o = jj
                  elseif(hdr(jj).eq.'1hdo') then
                       ihdo = jj
                  elseif(hdr(jj).eq.'1n2o') then
                       in2o = jj
                  elseif(hdr(jj).eq.'1ch4') then
                       ich4 = jj
                  endif
                  ii1=lspace+1
                  jj=jj+1
              endif
              if(flag.eq.1) then
c                 write(*,*)'flag is 1'
                  goto 77
              endif
          enddo ! do ii=1,9999 parse headers
77        continue

c         write(*,*)outfile
c         write(*,*)'ico2,ico,ih2o,ihdo,ich4,in2o=',
c    &               ico2,ico,ih2o,ihdo,ich4,in2o
          if(ico2.eq.0) then
              write(*,*)'Warning, 1co2 was not found.'
          endif
          if(ico .eq.0) then
              write(*,*)'Warning, 1co was not found.'
          endif
          if(ich4.eq.0) then
              write(*,*)'Warning, 1ch4 was not found.'
          endif
          if(in2o.eq.0) then
              write(*,*)'Warning, 1n2o was not found.'
          endif
          if(ih2o.eq.0) then
              write(*,*)'Warning, 1h2o was not found.'
          endif
          if(ihdo.eq.0) then
              write(*,*)'Warning, 1hdo was not found.'
          endif

c Determine which columns to output
          string1=hdr(ih2o)(2:4)//','//hdr(ihdo)(2:4)//','
     &    //hdr(ico2)(2:4)//','
     &    //hdr(in2o)(2:4)//','//hdr(ico)(2:3)//','//hdr(ich4)(2:4)
c         write(*,*)string1

c Write new .map file
          open(lun_out,file=outfile,status='unknown')
          write(lun_out,'(1x,a)')outfile
          write(lun_out,'(a)')mav_version
          write(lun_out,'(a)')version
          write(lun_out,'(1x,a,a)')'Height,Temp,Pressure,Density,',
     & string1

c Match units to string1
          write(lun_out,'(1x,a,a)')'km,K,hPa,molecules_cm3,',
     & 'parts,parts,ppm,ppb,ppb,ppb'

c    Note that the pe11.3 format without 0p affects the subsequent f11.3 format.
c    It multiplies it by 10!
10        format ((2(f7.2,","),1p,4(e10.3,","),0p,1(f8.3,","),
     & 1p,1(e10.3,","),0p,1(f8.3,","),f7.1))
          do j=1,nlev
             if(z(j).ge.0) then ! Skip cell levels
             write(lun_out,10)z(j),t(j),
     & p(j)*1013.25,d(j),
     & vmr(ih2o-4+(j-1)*nspeci),           ! H2O in parts
     & vmr(ihdo-4+(j-1)*nspeci),           ! HDO in parts
     & vmr(ico2-4+(j-1)*nspeci)*0.99*1e6,  ! CO2*0.99 in ppm
     & vmr(in2o-4+(j-1)*nspeci)*1e9,       ! N2O in ppb
     & vmr(ico -4+(j-1)*nspeci)*1e9,       ! CO  in ppb
     & vmr(ich4-4+(j-1)*nspeci)*1e9        ! CH4 in ppb
            endif
          enddo ! do j=1,nlev
          close(lun_out)

          read(lun_mav,'(a)',end=66) mavstring
          if(mavstring(1:14).eq.'Next Spectrum:') then
             read(mavstring(15:),'(a)') runlabmav
          else
             write(*,*) mavstring
             write(*,*)'Failed to find Next Spectrum string'
             stop
          endif ! Next Spectrum
      enddo ! do 1,9999
66    continue
      close(lun_mav)
      end
