c  subroutine READ_ONELINE_COL
c
c  Reads the .col output files (gas_1234.runlog.col) produced by GFIT and
c  writes it out to a string.
c
      subroutine read_oneline_col(output_string,luns, ncol,
     & csformat,lnit,ktg,specname_oof,
     & gfit_version,
     & gsetup_version, tllsum, solarsum)
      implicit none
      include "../gfit/ggg_int_params.f"
      include "params.f"

      integer idum,
     & lnbc,fnbc,fbc,
     & l2,l3,l4,dot,
     & ncol,icol,lnblnk
      integer lun_col,nss
      integer luns(mluns)

      character
     & col_string*500,
     & specname_col*(nchar),specname_oof*(nchar),
     & out_string*20000,
     & output_string*22000,
     & varr(5)*16,
c     & frmt*50,
     & gfit_version_number*16,
     & gsetup_version_number*16

      integer
     & jtg,mtg,lnit,nit,jj
      parameter(mtg=1000)
      integer*4
     & ktg(mcol)
      real*8
     & cl,tilt,cc,
     & fqshift,
     & sg,zlo,rmsfit,zmin,
     & airmass(mtg),ovcol(mtg)
      real*4
     & vsf(mtg),vsf_err(mtg)
      character
     & csform_array(mcol)*20,csform_condensed*128
      out_string = ""

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)
      idum=lun_qc   ! Avoid compiler warning (unused parameter)
      idum=lun_rlg  ! Avoid compiler warning (unused parameter)
      idum=lun_mul  ! Avoid compiler warning (unused parameter)
      idum=lnblnk(header)

c  Initialize character arrays (Necessary for the G77 compiler).
c     do i=1,mlabel
c        collabel(i:i)=' '
c     end do

c  Read in the retrieved absorber amounts (YOBS+-YERR)
      do icol=1,ncol     !  main loop (over windows)
c temporary set lun_col to read the correct file without changing
c too much code.

         lun_col = luns(icol)
c         write(*,*) 'lun_col = ',lun_col

         specname_col='='
         read(lun_col,'(a)',end=24) col_string
         if ( lnbc(col_string) .le. 2 ) go to 24 ! skip blank line at EOF
         l2=fbc(col_string(2:))+1 ! First space following spectrum name
         l3=fnbc(col_string(l2:))+l2-1 ! First character of NIT
         l4=fbc(col_string(l3:))+l3-1 ! First space following NIT

c This is only necessary because some of the .col file fields over-run their format, 
c so there is no space between numbers


c The following is a kludge to add in space (because of the lack of spaces):
         call substr(csformat,csform_array,mcol,nss)
c         write(*,*)'lnbc(csformat)=',lnbc(csformat)
c         write(*,*)'csformat=',csformat(:lnbc(csformat))
         csform_condensed = ''
         do jj=1,nss
            if (index(csform_array(jj),'1x').gt.0) then
               if (jj.eq.1) then
                  csform_condensed = 
     &            csform_condensed(:lnbc(csform_condensed))//
     &            csform_array(jj)(:lnbc(csform_array(jj)))
               else
                  csform_condensed = 
     &            csform_condensed(:lnbc(csform_condensed))//','//
     &            csform_array(jj)(:lnbc(csform_array(jj)))
               endif
            elseif (index(csform_array(jj-1),'1x').gt.0) then
               csform_condensed = 
     &         csform_condensed(:lnbc(csform_condensed))//','//
     &         csform_array(jj)(:lnbc(csform_array(jj)))
            elseif (index(csform_array(jj),'(').gt.0) then
               csform_condensed = 
     &         csform_condensed(:lnbc(csform_condensed))//','//
     &         csform_array(jj)(:index(csform_array(jj),'('))//'1x,'//
     &         csform_array(jj)
     &         (index(csform_array(jj),'(')+1:lnbc(csform_array(jj)))
            else
               csform_condensed = 
     &         csform_condensed(:lnbc(csform_condensed))//',1x,'//
     &         csform_array(jj)(:lnbc(csform_array(jj)))
            endif
         enddo ! jj=1,nss
         read(col_string,csformat) specname_col(:lnit-3),nit,cl,
     &   tilt,cc,
     &   fqshift,
     &   sg,zlo,rmsfit,zmin,
     &   (airmass(jtg),ovcol(jtg),vsf(jtg),vsf_err(jtg),jtg=1,ktg(icol))
c Check to ensure that the .col spectrum name is the same as the .oof
c spectrum name. Stop the program if this is not the case.
         dot=index(specname_col,'.')
         if(specname_col(:dot-2)//specname_col(dot:).ne.
     &   specname_oof(:dot-2)//specname_oof(dot:)) then
            write(*,*)'The spectrum name in the .col file does not '//
     &      'match the spectrum name in the .oof file: '
            write(*,*)specname_col(:lnit-3),', ',specname_oof(:lnit-3)
            write(*,*)'Please fix this problem and rerun write_eof.'
            stop
         endif
         write(col_string,csform_condensed)
     &   specname_col(:lnit-3),nit,cl,tilt,cc,
     &   fqshift,
     &   sg,zlo,rmsfit,zmin,
     &   (airmass(jtg),ovcol(jtg),vsf(jtg),vsf_err(jtg),jtg=1,ktg(icol))
c         write(*,*)'col_string after=',col_string(:lnbc(col_string))
c       else ! assume previous col file format
c         read(col_string,csformat) specname_col(:lnit-3),nit,cl,tilt,
c    &    fqshift,
c    &    sg,zlo,rmsfit,zmin,
c    &   (airmass(jtg),ovcol(jtg),vsf(jtg),vsf_err(jtg),jtg=1,ktg(icol))
c Check to ensure that the .col spectrum name is the same as the .oof
c spectrum name. Stop the program if this is not the case.
c          dot=index(specname_col,'.')
c          if(specname_col(:dot-2)//specname_col(dot:).ne.
c    &        specname_oof(:dot-2)//specname_oof(dot:)) then
c             write(*,*)'The spectrum name in the .col file does not '//
c    &        'match the spectrum name in the .oof file: '
c             write(*,*)specname_col(:lnit-3),', ',specname_oof(:lnit-3)
c             write(*,*)'Please fix this problem and rerun write_eof.'
c             stop
c          endif
c          write(col_string,csform_condensed)
c    &     specname_col,nit,cl,tilt,
c    &     fqshift,
c    &     sg,zlo,rmsfit,zmin,
c    &   (airmass(jtg),ovcol(jtg),vsf(jtg),vsf_err(jtg),jtg=1,ktg(icol))
c       endif


         out_string = out_string(:lnbc(out_string))
     &   //' '//col_string(l2:lnbc(col_string))

24       continue
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
     & //' '//tllsum//' '//solarsum
 
      end
