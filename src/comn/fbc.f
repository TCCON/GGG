      integer function fbc(string)

c  Returns the location (FBC) of the First Blank Character in STRING
c  Recognizes the following characters as "blanks":
c     nul             (ASCII character # 0)
c     horizontal tab  (ASCII character # 9) 
c     carriage return (ASCII character # 13)
c     space           (ASCII character # 32)
c     comma           (ASCII character # 44)
c
c  GCT 11-May-1998 
      implicit none
      integer ic
      character string*(*)
      do fbc=1,len(string)
         ic=ichar(string(fbc:fbc))
         if(    ic.eq.0
     &     .or. ic.eq.9
     &     .or. ic.eq.13
     &     .or. ic.eq.32
     &     .or. ic.eq.44) return   ! Successful return
      end do
      return   ! Abnormal return; no blanks found; fbc=LEN+1
      end
