C   Returns a uniform random deviate between 0.0 and 1.0.
C   Set idum to any nagative value to (re)initialize the sequence.
C     Re-coded per changes to Numerical Recipes
C     Bill Irion, May 14, 2003

        FUNCTION ran1(idum)
        implicit none
        INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
        REAL ran1,AM,EPS,RNMX
        PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     &   NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER j,k,iv(NTAB),iy
        SAVE iv,iy
        DATA iv /NTAB*0/, iy /0/
        if (idum.le.0.or.iy.eq.0) then
          idum=max(-idum,1)
          do j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iv(j)=idum
          enddo
          iy=iv(1)
        endif
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        j=1+iy/NDIV
        iy=iv(j)
        iv(j)=idum
        ran1=min(AM*iy,RNMX)
        return
        END
