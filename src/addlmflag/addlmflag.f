c  Program addlmflag 
c  Reads the linelist $ggpath/linelist/addlmflag_official.161
c  which contains the 1177 lines for which LM flags should be
c  present in the main linelist.
c  Program then reads and re-writes the main linelist to
c  addlmflag.out. For each line, it searches for a near-
c  identical line in addlmflag_official.161 and if it
c  finds one it set the LMflag to 1 before writing.
c
c  Since in a new linelist the position, intensities, 
c  widths and even quantum numbers may be different
c  to previously, it is not trivial to find the line
c  in LL2 that matches the one in LL1.
c
c  There are two ways that this could be done
c  1) The memory-efficient method reads only the small
c    linelist (addlmflag_official.161) into memory.
c   Then read/adjust/write the main linelist.
c   This risks multiple lines being adjusted in the
c   main linelist on the basis of a single line of
c   addlmflag_official.161. Or no lines being adjusted
c   if they are too different.
c  2) Read only the large main linelist into memory.
c   Then read through addlmflag_official.161, adjusting
c   the lmflags of the main linelist in memory, only
c   allowing one primary line to be adjusted per line
c   in addlmflag_official.161. Then write out the
c   partly-adjusted main linelist.
c
c  Program identifies matching lines based on their gas#,
c  isotope#, quantum#, position, strength, E"
c  Assumes a linelist naming convention  xxxxx.nnn  
c  where nnn is the number of characters per line
c     nnn = 101 (HITRAN_2000, Unix)
c     nnn = 102 (HITRAN_2000, DOS)
c     nnn = 161 (HITRAN_2004, Unix)
c     nnn = 162 (HITRAN_2004, DOS)
c
      implicit none
      integer*4 mline,nline2,i1,i2,i2best,nlmflag,
     & ls,le,posnvec,lrt,lnbc,ciso2kiso,
     & lunr_ll1,lunr_ll2,lunw_out
      parameter (mline=1200,lunr_ll1=14,lunr_ll2=15,lunw_out=16)

      integer*4 
     & molno1,kisot1,
     & molno2(mline),kisot2(mline)

      real*8 diff,dbest,dbthresh,wid,
     & abhw1,pshift1,sbhw1,tdabhw1,
     & freq1,stren1,eaco1,eprime1,

     & abhw2(mline),pshift2(mline),sbhw2(mline),tdabhw2(mline),
     & freq2(mline),stren2(mline),eaco2(mline),eprime2(mline)

      character version*48,fullrec1*160,
     & llpath1*256,llpath2*256,ciso1*1,ciso2*1,
     & llrformat*80,gggdir*256,dl*1
c============================================================
      version=' ADDLMFLAG    Version 1.01     2019-11-03    GCT'
      write(*,*) version
      write(*,*)

c     Platform specification:      DG090519
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)     !Length of gggdir

      dbthresh=0.25d0
      wid=1.0d0

      if (iargc() == 0) then
       write(6,'(a)') 'Enter full paths to linelist:'
         read(*,'(a)') llpath1
      elseif (iargc() == 1) then
         call getarg(1, llpath1)
      else
         stop 'Usage: $gggpath/bin/addlmflag path1 '
      endif

c  Read entire contents of LL2: the addlm_official.161 linelist,
c  into memory
      llrformat='(i2,a1,f12.6,2e10.3,2f5.4,f10.4,f4.2,f8.6)'

      llpath2=gggdir(:lrt)//'linelist/addlmflag_official.161'
      write(*,*)' Reading LL2: '//llpath2(:lnbc(llpath2))
      open(lunr_ll2,file=llpath2,status='old')
      do i2=1,mline
         read(lunr_ll2,llrformat,end=88)
     &   molno2(i2),ciso2,freq2(i2),stren2(i2),eaco2(i2),
     &   abhw2(i2),sbhw2(i2),eprime2(i2),tdabhw2(i2),pshift2(i2)
         kisot2(i2)=ciso2kiso(molno2(i2),ciso2) ! Wretched kludge for HITRAN CO2
      end do
      stop 'Increase parameter MLINE'
88    nline2=i2-1

c============================================================
c  Read & re-write primary linelist, making adjustments to
c  lmflag if found to be set in afflmflag_official.161.
      open(lunr_ll1,file=llpath1,status='old')
      open(lunw_out,file='addlmflag.out',status='unknown')
c
      nlmflag=0
      do i1=1,9999999  ! Main loop over lines in linelist1
         read(lunr_ll1,'(a)',end=99)fullrec1
         read(fullrec1,llrformat)
     &   molno1,ciso1,freq1,stren1,eaco1,
     &   abhw1,sbhw1,eprime1,tdabhw1,pshift1
         kisot1=ciso2kiso(molno1,ciso1) ! Wretched kludge for HITRAN CO2
c  Find best-matching line in secondary linelist.
c  Start/End search in linelist 1 at freq1 +/- wid cm-1
c         write(*,*)'posnvec=',posnvec(nline2,freq2,freq1)
         ls=max(     1,posnvec(nline2,freq2,freq1-wid))
         le=min(nline2,posnvec(nline2,freq2,freq1+wid))
c         write(*,*) 'freq1,ls,le = ',freq1,ls,le
         i2best=1
         dbest=1.E+32
c  Find best-matching line in linelist 2
         do i2=ls,le  ! Loop over linelist 2
            diff=sqrt(
     &      ((freq1-freq2(i2))/0.01)**2 +
     &      ((abhw1-abhw2(i2))/0.3)**2 +
     &      ((sbhw1-sbhw2(i2))/0.3)**2 +
     &      ((eprime1-eprime2(i2))/0.1)**2 +
     &      ((pshift1-pshift2(i2))/0.1)**2 +
     &      (dlog(stren1/stren2(i2)))**2 +
     &      100*(molno1-molno2(i2))**2 +
     &      100*(kisot1-kisot2(i2))**2)
            if(diff.lt.dbest) then
               i2best=i2
               dbest=diff
            endif
         enddo  ! i2=ls,le  ! Loop over linelist 2
c         if(i2best.le.0) stop 'i2best <= 0'
         if(dbest.lt.dbthresh) then
            write(*,*) 'Setting LM flag: ',i1,i2best,dbest
            fullrec1(146:146)='1'
            nlmflag=nlmflag+1
c            write(*,*) freq1,freq2(i2best)
c            write(*,*) stren1,stren2(i2best)
c            write(*,*) abhw1,abhw2(i2best)
c            write(*,*) sbhw1,sbhw2(i2best)
c            write(*,*) eprime1,eprime2(i2best)
c            write(*,*) molno1,molno2(i2best)
c            write(*,*) kisot1,kisot2(i2best)
c            write(*,*) 
         endif

         write(lunw_out,'(a160)')fullrec1
      enddo      !  do i1=1,nline1
99    close(lunw_out)
      write(*,*)'# lines found in LL2 = ',nline2
      write(*,*)'# times than lmflag was set = ',nlmflag
      if(nline2.ne.nlmflag) then
         write(*,*) 'Error: nline2 .ne. nlmflag: Adjust dbthresh'
      endif
      stop
      end
