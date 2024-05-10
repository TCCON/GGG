      subroutine lowercase(ss)
c  Performs an in-place conversion of string SS to lower-case.

      implicit none
      integer*4 i,ic
      character ss*(*)
c
      do i=1,len(ss)
         ic=ichar(ss(i:i))
         if(ic.ge.65 .and. ic.le.90) ss(i:i)=char(ic+32)
      end do
      return
      end
