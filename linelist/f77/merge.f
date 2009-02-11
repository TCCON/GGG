c  Program to merge two HITRAN-format linelists
      character list1*40,list2*40,string1*100,string2*100
      integer lunr1,lunr2,lunw
      parameter (lunr1=14,lunr2=15,lunw=16)
      real*8 nu1,nu2
c
      write(*,*)' Enter linelist #1'
      read(*,88)list1
      write(*,88)list1
      write(*,*)' Enter linelist #2'
      read(*,88)list2
      write(*,88)list2
88    format(a)
      open(lunr1,file=list1,status='old')
      read(lunr1,88)string1
      read(string1,'(3x,f12.6)')nu1
      open(lunr2,file=list2,status='old')
      read(lunr2,88)string2
      read(string2,'(3x,f12.6)')nu2
      open(lunw,file='merge.out',status='unknown')
1     if(nu1.lt.nu2) then
          write(lunw,88)string1
          read(lunr1,88,end=76)string1
          read(string1,'(3x,f12.6)')nu1
      else
          write(lunw,88)string2
          read(lunr2,88,end=66)string2
          read(string2,'(3x,f12.6)')nu2
      endif
      go to 1
c
c Reached end of file2 (lunr2): Finish up file 1
66    close(lunr2)
      write(lunw,88)string1
67    read(lunr1,88,end=99)string1
      write(lunw,88)string1
      go to 67
c
c Reached end of file1 (lunw): Finish up file 2
76    close(lunr1)
      write(lunw,88)string2
77    read(lunr2,88,end=99)string2
      write(lunw,88)string2
      go to 77
c
99    close(lunr2)
      stop
      end
