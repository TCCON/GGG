      subroutine skiprec(lun,nrec)
c  Skips NREC records of the already open logical unit number LUN
c  This avoids having multiple read(lun,*) statements in the main program
c  or avoids having to set up do loops to perform the multiple reads.
      implicit none
      integer*4 lun,irec,nrec
      do irec=1,nrec
         read(lun,*)
      end do
      return
      end
