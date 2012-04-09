c  subroutine READ_ONELINE_COL
c
c  Reads the .col output files (gas_1234.runlog.col) produced by GFIT and
c  writes it out to a string.
c
      subroutine read_oneline_col(output_string,luns, ncol,
     & csformat,lnit,ktg, 
     & gfit_version,
     & gsetup_version, atmsum, gctsum, fciasum, sciasum, solarsum)
      implicit none
      include "../ggg_int_params.f"
      include "params.f"

      integer 
     & lnbc,fnbc,fbc,
     & l2,l3,l4,
     & ncol,icol,mval
      integer lun_col,nss
      integer luns(mluns)
      parameter (mval=mrow*mcol) ! Max number of values (NROW * NCOL)

      character
     & col_string*500,
     & specname_col*(nchar),
     & out_string*9000,
     & output_string*10000,
     & varr(5)*16,
     & frmt*50,
     & gfit_version_number*10,
     & gsetup_version_number*10

      integer
     & jtg,noc,mtg,lnit,nit,jj
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
     & csform_array(mcol)*20,csform_condensed*120
      out_string = ""

c  Initialize character arrays (Necessary for the G77 compiler).
c     do i=1,mlabel
c        collabel(i:i)=' '
c     end do

c  Read in the retrieved absorber amounts (YOBS+-YERR)
      do icol=1,ncol     !  main loop (over windows)
c temporary set lun_col to read the correct file without changing
c too much code.

        lun_col = luns(icol)

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
c      write(*,*)'csformat=',csformat(:lnbc(csformat))
       csform_condensed = ''
       do jj=1,nss
              if (index(csform_array(jj),'1x').gt.0) then
                 if (jj.eq.1) then
              csform_condensed = 
     & csform_condensed(:lnbc(csform_condensed))//
     & csform_array(jj)(:lnbc(csform_array(jj)))
                 else
              csform_condensed = 
     & csform_condensed(:lnbc(csform_condensed))//','//
     & csform_array(jj)(:lnbc(csform_array(jj)))
                 endif
              elseif (index(csform_array(jj-1),'1x').gt.0) then
              csform_condensed = 
     & csform_condensed(:lnbc(csform_condensed))//','//
     & csform_array(jj)(:lnbc(csform_array(jj)))
              elseif (index(csform_array(jj),'(').gt.0) then
              csform_condensed = 
     & csform_condensed(:lnbc(csform_condensed))//','//
     & csform_array(jj)(:index(csform_array(jj),'('))//'1x,'//
     & csform_array(jj)
     & (index(csform_array(jj),'(')+1:lnbc(csform_array(jj)))
              else
              csform_condensed = 
     & csform_condensed(:lnbc(csform_condensed))//',1x,'//
     & csform_array(jj)(:lnbc(csform_array(jj)))
              endif
       enddo ! jj=1,nss
c      write(*,*)'csform_condensed=',csform_condensed

        if (index(gfit_version,'4.8.') .ne. 0) then ! assume new col file format
           read(col_string,csformat) specname_col(:lnit-3),nit,cl,
     &     tilt,cc,
     &     fqshift,
     &     sg,zlo,rmsfit,zmin,
     &   (airmass(jtg),ovcol(jtg),vsf(jtg),vsf_err(jtg),jtg=1,ktg(icol))
           write(col_string,csform_condensed)
     &     specname_col(:lnit-3),nit,cl,tilt,cc,
     &     fqshift,
     &     sg,zlo,rmsfit,zmin,
     &   (airmass(jtg),ovcol(jtg),vsf(jtg),vsf_err(jtg),jtg=1,ktg(icol))
c          write(*,*)'col_string after=',col_string(:lnbc(col_string))
        else ! assume previous col file format
           read(col_string,csformat) specname_col(:lnit-3),nit,cl,tilt,
     &     fqshift,
     &     sg,zlo,rmsfit,zmin,
     &   (airmass(jtg),ovcol(jtg),vsf(jtg),vsf_err(jtg),jtg=1,ktg(icol))
           write(col_string,csform_condensed)
     &     specname_col,nit,cl,tilt,
     &     fqshift,
     &     sg,zlo,rmsfit,zmin,
     &   (airmass(jtg),ovcol(jtg),vsf(jtg),vsf_err(jtg),jtg=1,ktg(icol))
        endif


           out_string = out_string(:lnbc(out_string))
     &     //' '//col_string(l2:lnbc(col_string))

24      continue
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
