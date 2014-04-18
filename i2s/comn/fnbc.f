      integer function fnbc(string)
c  Returns position of First Non-Blank Character in STRING
c  FNBC=0 indicates that the entire string was blank
c  GCT  9-Jul-1993
      implicit none
      integer ic
      character string*(*)
      do fnbc=1,len(string)
        ic=ichar(string(fnbc:fnbc))
        if( ic.ne.32 .and. ic.ne.0 ) return   ! if it's not NULL or SPACE
      enddo
      fnbc=0
      return
      end
