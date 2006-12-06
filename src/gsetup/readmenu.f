      subroutine readmenu(path,filename)
c  Reads and displays menu files and prompts user for their selection
c  Returns the name of the file (filename) chosen by the user.
      character header*80,string*76,path*(*),filename*(*)
      integer max,lunr,fbc,lnbc,k,kt
      parameter (max=999,lunr=19)
      write(*,*)'Readmenu:',path
      open(lunr,file=path,status='old')
 119  read(lunr,88)header
 88   format(a)
      write(6,*)
      write(6,88) '  #  '//header(:lnbc(header))
      do k=1,max
        read(lunr,'(a)',end=113) string
        write(6,'(i3,1x,a)') k,string
      end do
      write(6,*) 'warning! too many files: increase parameter max'
 113  rewind(lunr)
      write(6,9915)
 9915 format(' Enter # of file to be used ',$)
      read(5,*,err=119) kt
      if(kt.ge.k .or. kt.le.0) go to 119
      read(lunr,88) header
      do k=1,kt
        read(lunr,'(a)') string
      end do
      filename=string(:fbc(string)-1)
      close(lunr)
      return
      end
