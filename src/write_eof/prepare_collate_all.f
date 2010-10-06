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
      subroutine prepare_collate_all(header,luns, ncol, noc,
     & grl_array_size, specname_grl_array, iyr_array, doy_array,
     & delta_t_array, zpdtim_array, grl_array_counter, gfit_version,
     & gsetup_version)
      implicit none
      include "params.f"
      integer i1,idot,irow,k,ktg,lnbc,nn,
     & mlabel,mval,
     & lr,ncol,icol,
     & nlhead,iyrwas,doywas,jj,noc,
     & nss,mit,i,j,flag,lspace
      integer lun_col
      parameter (mval=mrow*mcol) ! Max number of values (NROW * NCOL)
      parameter (mlabel=16000)  ! Max # of characters in column labels

      integer*4 luns(mluns), temp_lun
      character cdum*20,colabel*15000,
     & gfit_version*80,gsetup_version*80,
     & collabel*(mlabel),
     & specname_grl*38,runlog*80,
     & tabel*80,
     & header*20000,
     & colfile*40,
c     & collate_version*64,
     & window(mcol)*10,specname_gwas*38,
     & hdri(mval)*(26)

      real*8 fcen,width,rmin,rmax,
     & r8was,trms

      real*8 zpdwas,max_delta_t

      logical append_qcflag
      logical append_spectrum_name

      integer grl_array_size
      character specname_grl_array(grl_array_size)*38
      integer iyr_array(grl_array_size), doy_array(grl_array_size)
      integer delta_t_array(grl_array_size)
      integer zpdtim_array(grl_array_size)
      integer grl_array_counter

      grl_array_counter=1
      specname_grl=""
      do grl_array_counter=1,grl_array_size
         specname_grl_array(grl_array_counter)=""
         iyr_array(grl_array_counter)=0
         doy_array(grl_array_counter)=0
         zpdtim_array(grl_array_counter)=0
      end do
      grl_array_counter=1

      append_qcflag=.false.
      append_spectrum_name=.true.
c     append_spectrum_name=.false.
c This program now probably requires the spectrum be appended
    
c      collate_version=
c     &' collate_all_results          Version 1.3.0   2010-06-24   GCT'
c      write(6,*) collate_version
      lr=0

c  Initialize character arrays (Necessary for the G77 compiler).
      do i=1,mlabel
         collabel(i:i)=' '
      end do

      do icol=1, mluns
         luns(icol)=0
      end do

c  Find the number of windows/columns (NCOL)
      open(lun_mul,file='multiggg.sh',status='old')
      do icol=1,mcol  
 1351    read(lun_mul,'(a)',end=99) tabel
        if(tabel(1:1).eq.':') go to 1351
         temp_lun=0
         call getlun(temp_lun)
         if (temp_lun .ne. 0) then 
             luns(icol) = temp_lun
         else    
            write(*,*) "Cannot find free lun to read in: ", tabel
            stop "Error preparing col files"
         endif
      end do  ! icol=1,mcol     !  main loop (over windows)

      read(lun_mul,*,end=99) tabel
      stop 'Increase parameter mcol'
 99   close(lun_mul)
      ncol=icol-1

c  Read in the retrieved absorber amounts (YOBS+-YERR)
      open(lun_mul,file='multiggg.sh',status='old')
      jj=1

      do icol=1,ncol     !  main loop (over windows)

        trms=0.0d0
135     read(lun_mul,'(a)') tabel
        if(tabel(1:1).eq.':') go to 135
        lun_col = luns(icol)

        colfile=tabel(index(tabel,'<')+1:index(tabel,'.ggg'))//'col'
        open(lun_col,file=colfile,status='old') ! .col file
        idot=index(colfile,'.')
        i1=index(colfile(:idot),'_')
        if(i1.eq.0) i1=index(colfile(:idot),'^')   !  new format
        if(i1.eq.0) i1=index(colfile(:idot),'-')   !  old format
        collabel=collabel(:lnbc(collabel)+2)//colfile(:idot-1)//' '//
     $  colfile(:idot-1)//'_error'
        window(icol)=colfile(:idot-1)
c
c  Read header lines of .col file and locate column containing "OVC_gas".
c  in order to read data from appropriate target gas.
        read(lun_col,*) nlhead
        read(lun_col,'(a)') gfit_version
        read(lun_col,'(a)') gsetup_version
        do k=4,nlhead-2
           if(nlhead.ge.23) read(lun_col,'(34x,a)')colabel
           if(nlhead.eq.20) read(lun_col,'(34x,a)')colabel
           if(nlhead.eq.21) read(lun_col,'(a)')colabel
           if(k.eq.6) runlog=colabel    ! GCT 2009-03-04
           if(index(colabel,'runlogs').gt.0) runlog=colabel
        end do

        read(lun_col,'(a)') colabel
        read(colabel,*) fcen, width, mit
        read(lun_col,'(a)')colabel
        k=1
        flag=0

        do j=1,99999
           if(colabel(k:k).eq.' ') then
               k=k+1
           else
               lspace=index(colabel(k:),' ')+k
!if I find the end of the line, then set flag=1
!this is used to end the loop later
               if(index(colabel(k:),' ')+k.ge.lnbc(colabel)) then
                   hdri(jj)=window(icol)(:lnbc(window(icol)))
     &                     //'_'//colabel(k:lnbc(colabel))
                   flag=1
                   jj=jj+1
               else
!ignore the spectrum column, otherwise append 
!the window name to the column headers
                   if (index(colabel(k:lspace),'Spectrum').ge.1) then
                   else
                   hdri(jj)=window(icol)(:lnbc(window(icol)))
     &                     //'_'//colabel(k:lspace)
                   jj=jj+1
                   endif
! 
              endif
               k=lspace+1
           endif
           if(flag.eq.1) then
              goto 56
           endif
        enddo ! j=1,99999

56      continue
        ktg=1+index(colabel,' OVC_'//colfile(:i1-1))
        if ( ktg .gt. 1) then
          call substr(colabel(:ktg-1),cdum,1,nss)
          call substr(colabel(ktg:lnbc(colabel)),cdum,1,noc)
          noc=(noc+1)/4
          ktg=(nss-4)/4
        endif
        iyrwas=-99999
        doywas=-99999
        zpdwas=-99999.9d0
        irow=0
        rmin=9999999.0
        rmax=0.0
        specname_grl=' '
        specname_gwas='x'

        lr=lnbc(runlog)
        if(runlog(lr-2:lr-2).eq.'o') then
           max_delta_t=0.0004  ! 1.44s (ACE)
        else
           max_delta_t=0.0025  ! 9.0s (ground-based TCCON)
        endif
      end do
c
c  Read auxilliary measurements from runlog
        open(lun_rlg,file=runlog(:lnbc(runlog)), status='old')   !DG000906
        read(lun_rlg,*) nlhead,nn
        do i=2,nlhead
           read(lun_rlg,*)
        end do
        r8was=-9999999.9d0

!this is the end of the header

c Write header to a string to pass along
      write(header,*) (hdri(k)(:lnbc(hdri(k))+1),k=1,jj-1)
      header=header(:lnbc(header))//' GFIT_Version GSETUP_Version'

      end
