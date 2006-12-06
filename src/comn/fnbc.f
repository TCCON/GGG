      integer function fnbc(string)
c  Returns position of First Non-Blank Character in STRING
c  LNBC=0 indicates that the entire string was blank
c  GCT  9-Jul-1993
c
c  Blanks include:
c     nul             (ASCII character # 0)
c     horizontal tab  (ASCII character # 9)
c     carriage return (ASCII character # 13)
c     space           (ASCII character # 32)
c     comma           (ASCII character # 44)

      implicit none
      integer ic
      character string*(*)
      do fnbc=1,len(string)
      ic=ichar(string(fnbc:fnbc))
      if(     ic.ne.0
     &  .and. ic.ne.9
     &  .and. ic.ne.13
     &  .and. ic.ne.32
     &  .and. ic.ne.44 ) return   !  It's not a blank
      end do
      fnbc=0
      return
      end
