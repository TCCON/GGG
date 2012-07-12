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
     &  gggdir(:lrt)//'bin'//dl//'collate_results < collate_t.input'
        open(lunw,file='collate_t.input', status='unknown')
        write(lunw,'(a)') 't'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results < collate_v.input'
        open(lunw,file='collate_v.input', status='unknown')
        write(lunw,'(a)') 'v'
        close(lunw)

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  '< average_results_t.input'
        open(lunw,file='average_results_t.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'tsw'
        close(lunw)

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  '< average_results_v.input'
        open(lunw,file='average_results_v.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'vsw'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'apply_airmass_correction < '//
     &  'apply_airmass_correction.input'
        open(lunw,file='apply_airmass_correction.input',
     &  status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'vav'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'apply_insitu_correction < '//
     &  'apply_insitu_correction.input'
        open(lunw,file='apply_insitu_correction.input',
     &  status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'vav.ada'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'apply_ghost_correction < '//
     &  'apply_ghost_correction.input'
        open(lunw,file='apply_ghost_correction.input',
     &  status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'vav.ada.aia'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'write_official_output_file < '//
     &  'write_official_output_file.input'
        open(lunw,file='write_official_output_file.input',
     &  status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'vav.ada.aia.gaa'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'write_eof < '//'write_eof.input'
        open(lunw,file='write_eof.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'tav'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'write_aux < '//'write_aux.input'
        open(lunw,file='write_aux.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'mav'
        write(lunw,'(f8.4)') oblat
        close(lunw)

c     'bal' or 'orb' geometry
      elseif (ext(1:3).eq.'bal' .or. ext(1:3).eq.'orb') then
        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results < collate_t.input'
        open(lunw,file='collate_t.input', status='unknown')
        write(lunw,'(a)') 't'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results < collate_l.input'
        open(lunw,file='collate_l.input', status='unknown')
        write(lunw,'(a)') 'l'
        close(lunw)

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  '< average_results_t.input'
        open(lunw,file='average_results_t.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'tsw'
        close(lunw)

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  '< average_results_l.input'
        open(lunw,file='average_results_l.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'lsw'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'zenang < '//'zenang.input'
        open(lunw,file='zenang.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'tav'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'zencorr < '//'zencorr.input'
        open(lunw,file='zencorr.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr)
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'diurnret < '//'diurnret.input'
        open(lunw,file='diurnret.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'lav'
        close(lunw)

c     'lab' geometry
      elseif (ext(1:3).eq.'lab') then
        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results < collate_t.input'
        open(lunw,file='collate_t.input', status='unknown')
        write(lunw,'(a)') 't'
        close(lunw)

        write(lunw_pp,'(a)')
     &  gggdir(:lrt)//'bin'//dl//'collate_results < collate_l.input'
        open(lunw,file='collate_l.input', status='unknown')
        write(lunw,'(a)') 'l'
        close(lunw)

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  '< average_results_t.input'
        open(lunw,file='average_results_t.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'tsw'
        close(lunw)

        write(lunw_pp,'(a)')gggdir(:lrt)//'bin'//dl//
     &  'average_results '//
     &  '< average_results_l.input'
        open(lunw,file='average_results_l.input',status='unknown')
        write(lunw,'(a)') rlgfile(:lr-3)//'lsw'
        close(lunw)

      endif
      close(lunw_pp)

      return
      end
