      subroutine get_ggg_environment(root, dl)
!   Returns the gggpath environment variable and delimiter (\ or /)

      integer*4      lrt,lnblnk 

      character*(*)   dl,       !filename path delimiter
     &              root        !the path to ggg

      call getenv('GGGPATH',root)
      lrt = lnblnk(root)
      if(lrt.le.0)then
         stop "GGGPATH environment variable not set"
      elseif(index(root,char(92)).ge.1)then
         dl = char(92)          !  \  PC
      else
         dl = char(47)          !  /  UNIX (default)
      endif

      if(root(lrt:lrt).ne.dl) root=root(1:lrt)//dl

      return
      end
