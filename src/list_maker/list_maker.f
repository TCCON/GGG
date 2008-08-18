c  Program list_maker.f
c
c  Generates a chronological list of files (spectra), given
c  user-supplied starting and ending file names.
c
c  The starting and ending filenames don't have to exist.
c  They merely represents the lowest and highest values of date,
c  detector, and extension.
c
c  Motivated by the failure of the Unix "ls" command, on many
c  computers, to generate a list when the directory contains
c  too many files.
c
c  Uses the gindfile command to search through the list of
c  data partitions where spectra may be stored.
c
      integer*4 lunw,fnbc,lnbc,lf,iyyyy,imm,idd,j1,j2,jul,
     & idet,idet_initial, idet_final,
     & ncount, nquery,lfi,lff,iext_initial,iext_final,iext
      parameter (lunw=14)
      character path*80,filename*40,initial_file*40,final_file*40,
     & det_initial*1,det_final*1,dplist*80

      call getenv('GGGPATH',dplist)
      dplist=dplist(:lnbc(dplist))//'/config/data_part_list_maker.lst'
      write(*,*) dplist

      write(*,*) 'Enter initial file (e.g. pa20040520saaaaa.001)'
      read(*,'(a)') initial_file
      lfi=fnbc(initial_file)
      if(lfi.gt.1) initial_file=initial_file(lfi:)

      write(*,*) 'Enter final file (e.g. pa20061231saaaab.500)'
      read(*,'(a)') final_file
      lff=fnbc(final_file)
      if(lff.gt.1) final_file=final_file(lff:)
      
      read(initial_file,'(2x,i4,i2,i2,5x,a1,1x,i3)')iyyyy,imm,idd,
     & det_initial,iext_initial
      idet_initial=ichar(det_initial)
      call julian(iyyyy,imm,idd,j1)

      read(final_file,'(2x,i4,i2,i2,5x,a1,1x,i3)')iyyyy,imm,idd,
     & det_final,iext_final
      idet_final=ichar(det_final)
      call julian(iyyyy,imm,idd,j2)

      open(lunw,file='list_maker.out',status='unknown')
      filename=initial_file
      lf=lnbc(filename)
      ncount=0
      nquery=0
      do jul=j1,j2
         call caldat(jul,iyyyy,imm,idd)
         write(filename(3:10),'(i4,i2.2,i2.2)')iyyyy,imm,idd
         do iext=iext_initial,iext_final
            do idet=idet_initial,idet_final
               write(filename(lf-4:lf),'(2a1,i3.3)') char(idet),'.',iext
               call gindfile(dplist,filename,path)
               nquery=nquery+1
               if(lnbc(path).gt.0) then
                    ncount=ncount+1
                    write(lunw,'(a)')filename
               endif
c              write(*,*)file_exists,filename
               if( mod(nquery,10000).eq.0 ) write(*,*) ncount,
     &         ' files listed   ',nquery,' queries performed',
     &         filename
            end do
         end do
      end do
      close(lunw)
      write(*,*) ncount,' files listed   ',nquery,' queries performed'
      stop
      end
