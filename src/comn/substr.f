      subroutine substr(inputstring,outputarray,mss,nss)
c  Converts a space/tab/comma-delimited character string (INPUTSTRING),
c  into an array of its NSS sub-string components (OUTPUTARRAY).
c
c  Inputs:
c       INPUTSTRING   The character string to be parsed.
c       MSS           Declared dimension of OUTPUTARRAY in calling prog.
c
c   Outputs:
c       NSS           Number of sub-strings found in INPUTSTRING.
c       OUTPUTARRAY   Array of the NSS sub-strings that were located.
c
c  Special Notes:
c  1) If the actual number of sub-strings (NSS) exceeds the delared dimension
c     of OUTPUTARRAY (MSS), only the first MSS sub-strings will be copied to
c     OUTPUTARRAY, avoiding the possibility of array-bound violations. We
c     therefore reccommend that each call of substr be followed by an if
c     statement such as:
c     if(nss.gt.mss) write(*,*) 'SUBSTR warning: Increase parameter MSS to',nss
c  2) If the length of any of the sub-strings exceeds LEN(OUTPUTARRAY) 
c     a warning is printed, showing the full oversized sub-string,
c     which is then truncated.
c  3) Currently recognized delimiters include
c                       nul            (ASCII character # 0)
c                       horizontal tab (ASCII character # 9)
c                       space          (ASCII character # 32)
c                       comma          (ASCII character # 44)
c  4) Calls external functions FNBC.F (First Non-Blank Character),
c     and FBC.F (First Blank Character)
c
c  18-May-98  GCT
      implicit none
      integer mss,nss,ibeg,iend,fnbc,fbc,lenout
      character inputstring*(*),outputarray(mss)*(*)
c
      lenout=len(outputarray(1))
      nss=0
      iend=0
      ibeg=fnbc(inputstring)
      do while (ibeg.gt.iend)
         nss=nss+1
         iend=ibeg+fbc(inputstring(ibeg+1:))
c         write(*,*)'substr: ',nss,ibeg,iend-1,inputstring(ibeg:iend-1)
         if(iend-ibeg.gt.lenout) then
           write(6,*) ' Warning from SUBSTR.F: sub-string too long:'
           write(6,*) inputstring(ibeg:iend-1)
         endif
         if(nss.le.mss) outputarray(nss)=inputstring(ibeg:iend-1)
         ibeg=iend+fnbc(inputstring(iend+1:))
      end do
      return
      end
