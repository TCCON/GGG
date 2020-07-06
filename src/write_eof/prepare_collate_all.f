c  Subroutine PREPARE_COLLATE_ALL
c
c  Reads the .col output file headers. Requires multiggg.sh to run.
c
      subroutine prepare_collate_all(header,luns, ncol,ktg,
     & gfit_version,
     & gsetup_version,tllsum,solarsum,lnit,
     & csformat)
      implicit none
      include "../gfit/ggg_int_params.f"
      include "params.f"

      integer i1,idot,idum,
     & k,lnbc,
     & mlabel,mval,
     & lr,ncol,icol,lnit,
     & nlhead,
     & jj,
c     & noc,
     & mit,i,j,flag,lspace
      integer lun_col
      parameter (mval=mrow*mcol) ! Max number of values (NROW * NCOL)
      parameter (mlabel=16000)  ! Max # of characters in column labels

      integer*4 luns(mluns), temp_lun,ktg(mcol),nss
      character
     & cdum*20,
     & rlgfile*(mfilepath),
     & colfile*64,
     & colabel*15000,
     & collabel*(mlabel),
     & tabel*500,
     & windows(mcol)*10,
     & hdri(mval)*(26),
     & cksum*32

      real*8 fcen,width,trms

c     logical append_qcflag
c     logical append_spectrum_name

c     integer grl_array_size
c     character specname_grl_array(grl_array_size)*(nchar)
c     integer iyr_array(grl_array_size), doy_array(grl_array_size)
c     integer zpdtim_array(grl_array_size)
c     integer grl_array_counter

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

      idum=lun_qc   ! Avoid compiler warning (unused parameter)
      idum=lun_rlg  ! Avoid compiler warning (unused parameter)

c     grl_array_counter=1
c     do grl_array_counter=1,grl_array_size
c        specname_grl_array(grl_array_counter)=""
c        iyr_array(grl_array_counter)=0
c        doy_array(grl_array_counter)=0
c        zpdtim_array(grl_array_counter)=0
c     end do
c     grl_array_counter=1

c     append_qcflag=.false.
c     append_spectrum_name=.true.
c This program now requires the spectrum be appended
    
      lr=0
c initialize the checksum values to zero, in case the linelists aren't found
      tllsum  ='00000000000000000000000000000000'
      solarsum='00000000000000000000000000000000'

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
         if(lnbc(tabel).le.0) go to 1351  ! blank line
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
135      read(lun_mul,'(a)') tabel
         if(tabel(1:1).eq.':') go to 135
         if(lnbc(tabel).le.0) go to 135  ! blank line
         lun_col = luns(icol)

         colfile=tabel(index(tabel,' ')+1:index(tabel,'.ggg'))//'col'
         open(lun_col,file=colfile,status='old') ! .col file
         idot=index(colfile,'.')
         i1=index(colfile(:idot),'_')
         if(i1.eq.0) i1=index(colfile(:idot),'^')   !  new format
         if(i1.eq.0) i1=index(colfile(:idot),'-')   !  old format
         collabel=collabel(:lnbc(collabel)+2)//colfile(:idot-1)//' '//
     $   colfile(:idot-1)//'_error'
         windows(icol)=colfile(:idot-1)
c
c  Read header lines of .col file and locate column containing "OVC_gas".
c  in order to read data from appropriate target gas.
         read(lun_col,*) nlhead
         read(lun_col,'(a)') gfit_version
         read(lun_col,'(a)') gsetup_version
         do k=4,nlhead-2
            read(lun_col,'(a32,2x,a)')cksum,colabel
            if(k.eq.6) rlgfile=colabel(:mfilepath)
            if(index(colabel,'runlogs').gt.0)rlgfile=colabel(:mfilepath)
            if(index(colabel,'telluric_linelists.md5').gt.0) then
               if(index(cksum,'  ').eq.0) tllsum=cksum
            endif
            if(index(colabel,'solar_merged.108').gt.0) then
               if(index(cksum,'  ').eq.0) solarsum=cksum
            endif
         end do

         csformat=colabel(:lnbc(colabel))

c         write(*,*)'csformat=',csformat
         read(lun_col,'(a)') colabel
         read(colabel,*) fcen, width, mit
         read(lun_col,'(a)')colabel
c         write(*,*)'colabel=',colabel
         lnit=index(colabel,'Nit')
c         write(*,*)'colabel=',colabel(:lnbc(colabel))
c         write(*,*)'colfile=',colfile(:i1-1)
         ktg(icol)=1+index(colabel,' AM_'//colfile(:i1-1))
         if ( ktg(icol) .gt. 1) then
            call substr(colabel(ktg(icol)-1:lnbc(colabel)),cdum,1,nss)
            ktg(icol)=nss/4
         endif
c         write(*,*)'ktg=',ktg

         k=1
         flag=0

         do j=1,99999
            if(colabel(k:k).eq.' ') then
               k=k+1
            else
               lspace=index(colabel(k:),' ')+k-1
c if I find the end of the line, then set flag=1
c this is used to end the loop later
               if(index(colabel(k:),' ')+k.ge.lnbc(colabel)) then
                  hdri(jj)=windows(icol)(:lnbc(windows(icol)))
     &            //'_'//colabel(k:lnbc(colabel))
                  flag=1
                  jj=jj+1
               else
c ignore the spectrum column, otherwise append 
c the window name to the column headers
                  if (index(colabel(k:lspace),'Spectrum').ge.1) then
                  else
                     hdri(jj)=windows(icol)(:lnbc(windows(icol)))
     &               //'_'//colabel(k:lspace)
                     jj=jj+1
                  endif
c 
               endif
               k=lspace+1
            endif
            if(flag.eq.1) then
               goto 56
            endif
         enddo ! j=1,99999

56       continue
      end do
!this is the end of the header

c Write header to a string to pass along
      write(header,*) (hdri(k)(:lnbc(hdri(k))+1),k=1,jj-1)
      header=header(:lnbc(header))//' GFIT_Version GSETUP_Version'
     & //' TELLURIC_checksum'
     & //' SOLAR_checksum'

      end
