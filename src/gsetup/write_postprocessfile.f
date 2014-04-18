      subroutine write_postprocessfile(ext,rlgfile,oblat)
c
c write the post-processing file based on geometry 'ext'
c
      implicit none
      include "../ggg_int_params.f"

      character
     &  ppfilename*19,     !post-processing batch file name
     &  ext*3,             !observation geometry
     &  rlgfile*40,        !name of occultation file
     &  gggdir*(mpath),    !root directory
     &  dl*1               !delimiter (='/' Unix, ='\' DOS)

      integer*4 lunw_pp, lunw, lrt, lr, lnbc
      real*8 oblat
      parameter (lunw=52)     ! for writing (general purpose)
      parameter (lunw_pp=76)  ! for post_processing.sh/post_processing.bat

      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)     !Length of root directory
      lr=lnbc(rlgfile)

      if (dl.eq.char(92)) then       !Back-Slash (\)
         ppfilename = 'post_processing.bat'
      else
         ppfilename = 'post_processing.sh'
      endif

      open(lunw_pp,file=ppfilename,status='unknown')

c     'gnd' geometry
      if (ext(1:3).eq.'gnd') then
        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results t'

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results v'

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  rlgfile(:lr-3)//'tsw'

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  rlgfile(:lr-3)//'vsw'

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'apply_airmass_correction '//
     &  rlgfile(:lr-3)//'vav'

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'apply_insitu_correction '//
     &  rlgfile(:lr-3)//'vav.ada'

c apply_ghost_correction is no longer needed if the LSE resampling IPP
c version is used!
c        write(lunw_pp,'(a)')
c     &  gggdir(:lrt)//'bin'//dl//'apply_ghost_correction '//
c     &  rlgfile(:lr-3)//'vav.ada.aia'

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'error_scale_factor '//
     &  rlgfile(:lr-3)//'vav.ada.aia'

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'write_official_output_file '//
     &  rlgfile(:lr-3)//'vav.ada.aia'

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'write_eof '//rlgfile(:lr-3)//'tav'

        write(lunw_pp,'(a,f8.4)')
     &  gggdir(:lrt)//'bin'//dl//'write_aux '//
     &  rlgfile(:lr-3)//'mav ',oblat

c     'bal' or 'orb' geometry
      elseif (ext(1:3).eq.'bal' .or. ext(1:3).eq.'orb') then
        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results t'

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results l'

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  rlgfile(:lr-3)//'tsw'

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  rlgfile(:lr-3)//'lsw'

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'zenang '//
     &  rlgfile(:lr-3)//'tav'

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'zencorr '//
     &  rlgfile(:lr)

        write(lunw_pp,'(a)')'# cp '//gggdir(:lrt)//
     & 'runlogs/'//ext(1:3)//dl//rlgfile(:lr)//' '//
     &  gggdir(:lrt)//'runlogs/'//ext(1:3)//dl//rlgfile(:lr)//'.bak'

        write(lunw_pp,'(a)')'# cp new_rl.out '//gggdir(:lrt)//
     & 'runlogs/'//ext(1:3)//dl//rlgfile(:lr)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'diurnret '//
     &  rlgfile(:lr-3)//'lav'
     
        write(lunw_pp,'(a)')'# cp '//gggdir(:lrt)//
     & 'vmrs/'//ext(1:3)//dl//rlgfile(:lr-3)//'vmr '//
     &  gggdir(:lrt)//'vmrs/'//ext(1:3)//dl//rlgfile(:lr-3)//'vmr.bak'

        write(lunw_pp,'(a)')'# cp '//gggdir(:lrt)//
     & 'vmrs/'//ext(1:3)//dl//rlgfile(:lr-3)//'vmr.new '//
     &  gggdir(:lrt)//'vmrs/'//ext(1:3)//dl//rlgfile(:lr-3)//'vmr'

        write(lunw_pp,'(a)')'# '//
     &  gggdir(:lrt)//'bin'//dl//'gsetup < '//'gsetup.input'

        write(lunw_pp,'(a)') '# cp mmm multiggg.sh'


c     'lab' geometry
      elseif (ext(1:3).eq.'lab') then
        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results t'

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results l'

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  rlgfile(:lr-3)//'tsw'

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  rlgfile(:lr-3)//'lsw'

      endif
      close(lunw_pp)

      return
      end
