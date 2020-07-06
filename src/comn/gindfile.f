      subroutine gindfile(dplist,filnam,path)

c   GINDFILE     Version 3.1.2     GCT    12-Nov-99
c   Locates a file and returns the path to it (including the filename).
c   Systematically searches through the list of partitions, starting at
c   the partition where the previous spectrum was found.
c
c   Input:
c          DPLIST:  C*(*)  the path to the data partition list
c          FILNAM:  C*(*)  the name of the file to locate
c
c   Output:
c          PATH:    C*(*)  the full pathname to the file (including FILNAM)
c
c   Error Return:
c          PATH = ' '  (lnbc(path)=0)  if FILNAM was not found
c
c  Assumptions:
c       1)    LUN 19 is available at time of call
c       3)    Partition names are 80 characters or less.
c       4)    Partition names beginning with a ':' are skipped
c       5)    There are no more than MPART (=50) uncommented partitions
c
      implicit none
      integer*4
     $ lf,       ! length of filnam string
     & fq, lq,   ! Positions of first and last ? in partition
     $ lnbc,     ! external function (Last Non Blank Character)
     $ lloc,     ! external function (Reverse-Index)
     $ lunr,     ! logical unit number of M4 partition list
     $ j,        ! dummy partition index
     $ ipart,    ! partition index
     $ npart,    ! number of partitions to be searched
     $ mpart     ! maximum supported number of partitions
      parameter (lunr=19,mpart=256)
      character dplist*(*),filnam*(*),path*(*),
     & partition(mpart)*(256)
      logical*4 flexst
      save ipart,npart,partition
      data ipart,npart/1,0/
c======================================================================
c
      lf=lnbc(filnam)
      if(lf.gt.0) then           ! skip empty file names
c
c  If first call, read file containing list of data partitions to be searched.
         if(npart.lt.1) then
            open(lunr,file=dplist,status='old')
            do j=1,mpart
 1             read(lunr,'(a)',end=88) partition(j)
               if(partition(j)(1:1).eq.':') go to 1
            end do
            write(*,*)'GINDFILE warning: Increase MPART=',mpart
88          npart=j-1
            close(lunr)
c            ipart=1
         endif
c
c  Systematically search over the NPART partitions starting at IPART
         do j=1,npart   !  loop over NPART partitions
            fq=index(partition(ipart),'?')  ! location of first ?
            if(fq.gt.0) then
               lq=lloc(partition(ipart),'?')  ! location of last ?
               path=partition(ipart)(:fq-1)//filnam(1:1+lq-fq)//
     &         partition(ipart)(lq+1:lnbc(partition(ipart)))//
     &         filnam(:lf)
            else
               path=partition(ipart)(:lnbc(partition(ipart)))//
     &         filnam(:lf)
            endif
            inquire(file=path,exist=flexst)
c            if(flexst) write(*,*) path(:lnbc(path)),lf
            if(flexst) return              !  success exit  !
            ipart=mod(ipart,npart)+1       !  increment partition index
         end do
      endif  !  (lf.gt.0) 
c
c   Note that after an unsucessful exit, IPART has the same value that it had
c   on entry, having been incremented NPART times and wrapped around once.
      path=' '
      return                         !  failure exit  !
      end
