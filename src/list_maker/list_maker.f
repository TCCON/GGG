c  Program list_maker.f
c
c  Generates a chronological list of files (spectra), given
c  user-supplied starting and ending file name templates.
c
c  The starting and ending filenames don't have to exist.
c  They merely represents the lowest and highest values of date,
c  detector, and extension.
c
c  Motivated by the failure of the Unix "ls" command, on many
c  computers, to generate a list when the directory contains
c  too many files.
c
c  Uses the gindfile subroutine to search through the list of
c  data partitions where spectra may be stored.
c
      include "../gfit/ggg_int_params.f"

      integer*4 lunw,fnbc,lnbc,lf,iyyyy,imm,idd,j1,j2,jul,
     & idet,idet_initial, idet_final,
     & ncount, nquery,lfi,lff,iext_initial,iext_final,iext
      parameter (lunw=14)
      character path*(mpath),filename*57,initial_file*57,final_file*57,
     & det_initial*1,det_final*1,gggdir*(mpath),dl*1,fmt*10,version*50
      integer*4 ldot, lloc, nextdigit, idum

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

      version= ' list_maker    Version 1.11     2016-04-02   GCT '
      write(*,*) version
      call get_ggg_environment(gggdir, dl)
      gggdir=gggdir(:lnbc(gggdir))//'config'//dl//
     &  'data_part_list_maker.lst'
      write(*,*) gggdir

      if (iargc() == 0) then
         write(*,*) 'Enter initial file (e.g. pa20061101saaaaa.001)'
         read(*,'(a)') initial_file

         write(*,*) 'Enter final file (e.g. pa20061231saaaab.500)'
         read(*,'(a)') final_file
      elseif (iargc() == 2) then
         call getarg(1, initial_file)
         call getarg(2, final_file)
      else
         stop 'Usage: $gggpath/bin/list_maker initial_file final_file'
      endif
      lfi=fnbc(initial_file)
      if(lfi.gt.1) initial_file=initial_file(lfi:)
      lff=fnbc(final_file)
      if(lff.gt.1) final_file=final_file(lff:)
      
c      read(initial_file,*) firstline
c      ldot = lloc(firstline, ".")
c      lf=lnbc(firstline)
c      nextdigit = lf - ldot
c      write(fmt, '("(I", I1, ")")') nextdigit
c      write(*,*) 'fmt=',fmt
c      read(firstline,'(2x,i4,i2,i2,5x,a1)') iyyyy,imm,idd,det_initial
c      read(firstline(ldot+1:lf), fmt) iext_initial
c      idet_initial=ichar(det_initial)
c      call julian(iyyyy,imm,idd,j1)

      ldot = lloc(initial_file, ".")
      lf=lnbc(initial_file)
      nextdigit = lf - ldot
      write(fmt, '("(I", I1, ")")') nextdigit
      write(*,*) 'fmt=',fmt
      read(initial_file,'(2x,i4,i2,i2,5x,a1)') iyyyy,imm,idd,det_initial
      read(initial_file(ldot+1:lf), fmt) iext_initial
      idet_initial=ichar(det_initial)
      call julian(iyyyy,imm,idd,j1)

c      read(final_file,*) firstline
c      ldot = lloc(firstline, ".")
c      lf=lnbc(firstline)
c      read(firstline,'(2x,i4,i2,i2,5x,a1)') iyyyy,imm,idd,det_final
c      read(firstline(ldot+1:lf), fmt) iext_final
c      idet_final=ichar(det_final)
c      call julian(iyyyy,imm,idd,j2)

      ldot = lloc(final_file, ".")
      lf=lnbc(final_file)
      read(final_file,'(2x,i4,i2,i2,5x,a1)') iyyyy,imm,idd,det_final
      read(final_file(ldot+1:lf), fmt) iext_final
      idet_final=ichar(det_final)
      call julian(iyyyy,imm,idd,j2)

      open(lunw,file='list_maker.out',status='unknown')
      filename=initial_file
      lf=lnbc(filename)
      ncount=0
      nquery=0
c      write(fmt, '(i1)') nextdigit
      do jul=j1,j2
         call caldat(jul,iyyyy,imm,idd)
         write(filename(3:10),'(i4,i2.2,i2.2)')iyyyy,imm,idd
         do iext=iext_initial,iext_final
            do idet=idet_initial,idet_final
               write(filename(lf-(nextdigit+1):lf-nextdigit),'(2a1)')
     &           char(idet),'.'
               write(filename(lf-(nextdigit-1):lf),
     &          "(I" // char(48+nextdigit) // "." //
     &         char(48+nextdigit) // ")") iext
c     &          "(I" // ADJUSTL(fmt) // "." // ADJUSTL(fmt) // ")") iext
               call gindfile(gggdir,filename,path)
               nquery=nquery+1
               if(lnbc(path).gt.0) then
                  ncount=ncount+1
                  write(lunw,'(a)')filename
               endif
c              write(*,*)file_exists,filename
               if( mod(nquery,10000).eq.0 ) write(*,*) ncount,
     &         ' files listed   ',nquery,' queries performed ',
     &         filename
            end do
         end do
      end do
      close(lunw)
      write(*,*) ncount,' files listed   ',nquery,' queries performed '
      stop
      end
