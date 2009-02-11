c  Strips out the lines of a given gas & isotope from a prescribed
c  spectral interval from a HITRAN-format linelist.  Two output files
c  are created:
c  (i) the stripped-out lines (strippedNN_llname), 
c  (ii) the remaining lines (leftoverNN_llname).
c  where llname is the name of the input linelist and NN is the
c  target gas number.
c
      character string*100,inputfile*48,outputfile*52
      integer igas,ig,lunr,lunw,lunw2,lunw1,iso,is,iline,nstripped
      real*8 nus,nue,freq
      parameter (lunr=14,lunw2=15,lunw1=16)
c
      iline=0
      write(*,*)' Enter Input file name (e.g. /ggg/linelist/atm.h92)'
      read(*,*)inputfile
      outputfile='00_'//inputfile
      write(*,*)' Enter the gas and isotope numbers (e.g. 6 1)'
      write(*,*)' Note: -ve values strip ALL gases/isotopomers'
      read(*,*)igas,iso
      if (igas.lt.0) igas=0
88    format(a)
      write(*,*)' Starting & Ending frequencies (e.g. 500 7973.1):'
      read(*,*) nus,nue
      write(outputfile(1:2),'(i2.2)')igas
      open(lunr,file=inputfile,status='old')
      open(lunw1,file='leftover'//outputfile,status='unknown')
      open(lunw2,file='stripped'//outputfile,status='unknown')

      nstripped=0
      do iline=0,9999999
         read(lunr,88,end=99)string
         read(string,'(i2,i1,f12.6)')ig,is,freq
         lunw=lunw1
         if(freq.gt.nus .and. freq.lt.nue) then
            if(ig.eq.igas .or. igas.le.0) then
               if(is.eq.iso .or. iso.lt.0) then
                  nstripped=nstripped+1
                  lunw=lunw2
               endif
            endif
         endif
         write(lunw,88)string
         if(mod(iline,100000).eq.0) write(*,*)iline
      end do
c
99    close(lunr)
      close(lunw1)
      close(lunw2)
      write(*,*)'stripped'//outputfile,nstripped
      write(*,*)'leftover'//outputfile,iline-nstripped
      stop
      end
