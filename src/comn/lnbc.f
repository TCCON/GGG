      integer function lnbc(string)
c  Returns location (LNBC) of Last Non-Blank Character in STRING
c  Recognizes the following characters as "blanks":
c     nul            (ASCII character # 0)
c     horizontal tab (ASCII character # 9) 
c     carriage return (ASCII character # 13)
c     space          (ASCII character # 32)
c     comma          (ASCII character # 44)
c
c  GCT  4-Feb-1993
c  GCT  9-Jul-1993
c  GCT 17-Oct-1997  
c  GCT 11-May-1998   Added comma and horizontal tab  
c  GCT 11-Jun-2006   Added carriage return
c
c
      implicit none
      integer ic
      character string*(*)
      do lnbc=len(string),1,-1
         ic=ichar(string(lnbc:lnbc))
         if(     ic.ne.0
     &     .and. ic.ne.9
     &     .and. ic.ne.13
     &     .and. ic.ne.32
     &     .and. ic.ne.44) return   ! Successful return
      end do
      lnbc=0   ! Entire string was blanks.
      return   ! Abnormal return; no non-blanks found.
      end
