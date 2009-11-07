      subroutine get_ggg_environment(root, dl)
!   Returns the gggpath environment variable and delimiter (\ or /)

      integer*4      lrt 

      character*(*)   dl,       !filename path delimiter
     &              root        !the path to ggg

      call getenv('GGGPATH',root)
      lrt = lnbc(root)
      if(lrt<=0)then
          stop "GGGPATH environment variable not set"
      elseif(index(root,char(92)).ge.1)then
          dl = char(92)          !  \  PC
      else
          dl = char(47)          !  /  UNIX (default)
      endif

      if(root(lrt:lrt).ne.dl) root=root(1:lrt)//dl

      return
      end
