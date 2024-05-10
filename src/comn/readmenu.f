      subroutine readmenu(lunr,menupath,lun_stdin,filename)
c  Reads and displays menu files and prompts user for their selection
c  Returns the name of the file (filename) chosen by the user.
      implicit none
      character header*80,string*76,menupath*(*),filename*(*)
      integer nmax,lunr,fbc,lnbc,k,kt,ls,lun_stdin
      parameter (nmax=999)
      write(*,*)'Readmenu:',menupath
      open(lunr,file=menupath,status='old')
 119  read(lunr,'(a)')header
      write(*,'(a)') '  #  '//header(:lnbc(header))
      do k=1,nmax
         read(lunr,'(a)',end=113) string
         ls=lnbc(string)
         write(*,'(i3,1x,a)') k,string(1:ls)
      end do
      write(*,*) 'Warning! too many files: increase parameter nmax'
 113  rewind(lunr)
      write(*,'(a),$') ' Enter # of file to be used '
      read(lun_stdin,*,err=119) kt
      if(kt.ge.k .or. kt.le.0) go to 119
      read(lunr,'(a)') header
      do k=1,kt
         read(lunr,'(a)') string
      end do
      filename=string(:fbc(string)-1)
      close(lunr)
      return
      end
