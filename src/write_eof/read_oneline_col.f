c  Program COLLATE_ALL
c
c  Reads the .col output files (gas_1234.runlog.col) produced by GFIT and
c  writes it out to a string.
c
c  INPUT FILES:
c     multiggg.sh          batch file containing names of .col files
c     runlog.xrl           runlog file
c     gas_1234.runlog.col  files containing column amounts
c
      subroutine read_oneline_col(output_string,luns, ncol, noc,
     & grl_array_size, specname_grl_array, iyr_array, doy_array,
     & delta_t_array, zpdtim_array, grl_array_counter,gfit_version,
     & gsetup_version, atmsum, gctsum, fciasum, sciasum, solarsum)
      implicit none
      include "params.f"

      integer irow,jtg,lnbc,fnbc,fbc,
     & npp,ldot,l2,l3,l4,totnit,ntc,mlabel,jval,
     & ncol,icol,mval,
     & iyrwas,doywas,jj,noc,
     & nrow,nit,mit,i,iyr,doy,istat
      integer lun_col
      integer luns(mluns)
      parameter (mval=mrow*mcol) ! Max number of values (NROW * NCOL)
      parameter (mlabel=16000)  ! Max # of characters in column labels

      character
     & apf*2,
     & col1*1,
     & colfile*40,
     & csformat*90,
     & col_string*500,
     & collabel*(mlabel),
     & specname_grl*57,
     & specname_col*57,
c     & collate_version*64,
     & specname_gwas*57,
     & out_string*8000,
     & output_string*10000,
     & varr(5)*16,
     & gfit_version_number*10,
     & gsetup_version_number*10

      real*8 airmass(mval),asza,cl,tilt,cc,zlo,zobs,rmin,rmax,
     & fqshift,graw,obslat,obslon,opd,ovcol(mval),rmsfit,
     & trms,lasf,wavtkr,aipl,sia,fvsi,azim,wspd,wdir,osds

      real*4
     & vsf(mval),vsf_err(mval),ymiss
      parameter (ymiss=9.8765e+29)

      real*8 zpdtim,tout,pout,hout,tins,pins,hins,fovi,fovo,amal,
     & snr,zenoff,zoff,sg,zpdwas,max_delta_t,delta_t,zmin
      integer bytepw,ifirst,ilast,possp,nss

      logical append_qcflag
      logical append_spectrum_name

      integer grl_array_size
      character specname_grl_array(grl_array_size)*57
      integer iyr_array(grl_array_size), doy_array(grl_array_size)
      integer delta_t_array(grl_array_size)
      integer zpdtim_array(grl_array_size)
      integer grl_array_counter, gac
      logical specname_found


      out_string = ""
      nrow=1
      append_qcflag=.false.
      append_spectrum_name=.true.
c     append_spectrum_name=.false.
c This program now probably requires the spectrum be appended
    
c      collate_version=
c     &' collate_all_results          Version 1.3.0   2010-06-24   GCT'
c!      write(6,*) collate_version

c  Initialize character arrays (Necessary for the G77 compiler).
      do i=1,mlabel
         collabel(i:i)=' '
      end do

c  Read in the retrieved absorber amounts (YOBS+-YERR)
      jj=1
      do icol=1,ncol     !  main loop (over windows)
        npp=0
        trms=0.0d0
        totnit=0
        ntc=0
        jval=icol-ncol
c temporary set lun_col to read the correct file without changin
c too much code.

        lun_col = luns(icol)


        specname_col='='
        read(lun_col,'(a)',end=24) col_string
        if ( lnbc(col_string) .le. 2 ) go to 24 ! skip blank line at EOF
        l2=fbc(col_string(2:))+1 ! First space following spectrum name
        l3=fnbc(col_string(l2:))+l2-1 ! First character of NIT
        l4=fbc(col_string(l3:))+l3-1 ! First space following NIT

c       write(*,*)'gfit_version: ',gfit_version
        if (index(gfit_version,'2.40.2') .ne. 0) then !  old col file format 
           csformat='(1x,a21,i2,1x,f5.3,3(1x,f4.1),1x,f5.3,'
     &          //'1x,f6.4,f7.3,1x,9(0pf7.3,1pe10.3,0pf9.4,1pe8.1))'
        elseif (index(gfit_version,'4.8.') .ne. 0) then ! assume new col file format 
           write(csformat,'(a,i2.2,a)')'(1x,a',l4-5,
     &         ',i3,f6.3,4f5.1,f6.3,f7.4,f8.3,15(f7.3,e11.4,f9.4,e8.1))'
        else ! assume previous col file format
           write(csformat,'(a,i2.2,a)')'(1x,a',l4-5,
     &         ',i3,f6.3,3f5.1,f6.3,f7.4,f8.3,15(f7.3,e11.4,f9.4,e8.1))'
        endif

        if (index(gfit_version,'4.8.') .ne. 0) then ! assume new col file format 
           read(col_string,csformat) specname_col,nit,cl,tilt,cc,
     &     fqshift,
     &     sg,zlo,rmsfit,zmin,
     &     (airmass(jtg),ovcol(jtg),vsf(jtg),vsf_err(jtg),jtg=1,noc)
        else ! assume previous col file format
           read(col_string,csformat) specname_col,nit,cl,tilt,
     &     fqshift,
     &     sg,zlo,rmsfit,zmin,
     &     (airmass(jtg),ovcol(jtg),vsf(jtg),vsf_err(jtg),jtg=1,noc)
        endif
           totnit=totnit+nit
           if(nit.lt.mit) ntc=ntc+1  ! Number of Times Converged
           if(rmsfit.le.0.0) then
              write(*,*) 'rmsfit <= 0', colfile,irow
              write(*,*)specname_col,nit,cl,tilt,fqshift,sg,zlo,rmsfit
c              stop 'rmsfit <= 0'   ! Commented 2009-03-18
           endif
           if(rmsfit.gt.rmax) rmax=rmsfit
           if(rmsfit.lt.rmin) rmin=rmsfit
24         continue
           


           
!     loop until we find the specname, 
!     usually this will be the first test
           if (specname_col.ne.specname_grl) then
              specname_found=.false.              
              
              do while (.not. specname_found)
                 gac=grl_array_counter
                 if(gac.le.0) gac=grl_array_size
                 
                 do while(
     &              (gac.ne.(grl_array_counter+1)).and..not.
     &               ((grl_array_counter.eq.grl_array_size).and.
     &             (gac.eq.1)).and..not.specname_found
     &              )
!                  write(*,*) "L154: ", 
!     & grl_array_size, grl_array_counter,gac
!                  write(*,*) specname_col, specname_grl_array(gac)
                    if(specname_col.eq.
     &                      specname_grl_array(gac)) then
                       specname_found=.true.

                       specname_grl = specname_grl_array(gac)
                       iyr = iyr_array(gac)
                       doy = doy_array(gac)
                       zpdtim = zpdtim_array(gac)

                    else
                       gac=gac-1
                       if(gac.eq.0) gac=grl_array_size
                    endif
                 end do
!                  write(*,*) "l155"

                 do while(.not.specname_found)
                   !not in the array, read some more
!                     write(*,*) "l158: ", grl_array_counter
                    
                 call read_runlog(lun_rlg,col1,specname_grl,iyr,doy,
     &                zpdtim,
     &                obslat,obslon,zobs,asza,zenoff,azim,osds,opd,
     &                fovi,fovo,amal,ifirst,ilast,
     &                graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &                tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,
     &                aipl,istat)
              
                 specname_grl_array(grl_array_counter) = specname_grl
                 iyr_array(grl_array_counter) = iyr
                 doy_array(grl_array_counter) = doy
                 zpdtim_array(grl_array_counter) = zpdtim
                 

!                  write(*,*) "l178", grl_array_counter,
!     & specname_grl_array(grl_array_counter)
                 grl_array_counter = grl_array_counter+1
                 if(grl_array_counter.gt.grl_array_size) 
     &                grl_array_counter=1

                 if (specname_col.eq.specname_grl) then
                    specname_found=.true.              
                 endif

                 end do
              end do
!               write(*,*) "L181", specname_col, specname_grl
!!!!!!!!!!!!!!!
!              if(istat.ne.0) go to 14     ! Exit Loop  irow=1,mrow

              ldot=index(specname_grl,'.')
              delta_t=zpdtim-zpdwas
c
c  Create speparate entries in the .vsw file for spectra whose
c  ZPD times differ by more than MAX_DELTA_T
              if(iyr.ne.iyrwas .or. 
     &         doy.ne.doywas .or.
     &        dabs(delta_t).ge.max_delta_t ) then  ! New observation
c
c  The following if-statement shouldn't be necessary. But occasionally
c  you get simultaneous InGaAs/Si scans with very different ZPD times.
c  You don't want them to have separate entries in the .vsw file. So....
                 if(specname_grl(4:ldot-2).ne.specname_gwas(4:ldot-2)
     &           .or.  specname_grl(ldot:).ne.specname_gwas(ldot:)) then
                    irow=irow+1
                 endif
              endif
              iyrwas=iyr
              doywas=doy
              zpdwas=zpdtim
              specname_gwas=specname_grl
           endif     ! while(specname_col.ne.specname_grl)
        
           out_string = out_string(:lnbc(out_string))
     &     //' '//col_string(l2:lnbc(col_string))
           npp=npp+1
           trms=trms+1.0d0/rmsfit**2

        if(ncol*nrow.gt.mval) then
           write(*,*)ncol,nrow
           write(*,*)'Increase parameter MVAL to ',nrow*ncol
           stop 'collate'
        endif
      end do  ! icol=1,mcol     !  main loop (over windows)

c Select out number from gfit_version string
      call substr(gfit_version(:lnbc(gfit_version)),varr,5,nss)
      gfit_version_number = varr(3)
      call substr(gsetup_version(:lnbc(gsetup_version)),varr,5,nss)
      gsetup_version_number = varr(3)
c Write output to a string to pass along
      output_string=""
      write(output_string,*) out_string(:lnbc(out_string))
     & //' '//gfit_version_number(:lnbc(gfit_version_number))
     & //' '//gsetup_version_number(:lnbc(gsetup_version_number))
     & //' '//atmsum//' '//gctsum//' '//fciasum//' '//sciasum
     & //' '//solarsum
 
      end
