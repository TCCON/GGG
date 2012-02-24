      subroutine cdsubstr(inputstring,outputarray,mss,nss)
c  Converts a comma-delimited character string (INPUTSTRING),
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
c
c  17-Aug-2006  GCT
c
      implicit none
      integer mss,nss,idiff,iwas,inext,lenout
      character inputstring*(*),outputarray(mss)*(*)
c
      lenout=len(outputarray(1))
      nss=0
      idiff=1
      inext=0
      do while (idiff.gt.0)
         iwas =inext
         idiff=index(inputstring(iwas+1:),',')
         inext=iwas+idiff
         nss=nss+1
c         write(*,*)nss,iwas,inext,' ',inputstring(iwas+1:inext-1)
         if(idiff-1.gt.lenout) then
           write(6,*) ' Warning from SUBSTR.F: sub-string too long:'
           write(6,*) idiff-1,' ',inputstring(iwas+1:inext-1)
         endif
         if(nss.le.mss) outputarray(nss)=inputstring(iwas+1:inext-1)
      end do
      outputarray(nss)=inputstring(iwas+1:) ! Sub-string following last comma
      return
      end
