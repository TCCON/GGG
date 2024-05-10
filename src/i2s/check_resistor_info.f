      subroutine check_resistor_info(code)

      implicit none

      integer*4
     & code,
     & mch,        ! Maximum number of data channels
     & npgn,       ! Number of gain settings in pre-amp
     & mrci,       ! Maximum R_Config_Index (revisions of PA gains)
     & lunw        ! Logical Unit Number for output file
      parameter (mch=2)   ! FIXME: this should not be repeated
      parameter (npgn=4)
      parameter (mrci=3)
      parameter (lunw=21)

      integer*4
     & resistors(npgn,mch,mrci), ! 3-D matrix to hold all R values
     & ichan,    ! Channel number (1=InGaAs=slave, 2=Si=master)
     & ipgn,       ! Index into Pre-amp GaiNs
     & lnbc        ! Integer function Last Non-Blank Character in string

      character
     & detname(mch)*6 ! Identifiers for each channel (detector name)

      data resistors /1500,1820,750000,2700,    ! InGaAs, revision 1
     &                2800,3320,4700000,20000,  ! Si, revision 1
     &                3090,3650,750000,2700,    ! InGaAs, revision 2
     &                7680,9090,4700000,20000,  ! Si, revision 2
     &                3090,3650,750000,2700,    ! InGaAs, revision 3
     &                7680,9090,1500000,20000/  ! Si, revision 3
      data detname /'InGaAs','Si    '/

      open(lunw,file='resistors.h',status='unknown')
c
c  If input error, produce a header file that will cause a compilation
c  error.  Also send a message to the standard output.
      if((code.le.0).or.(code.gt.mrci)) then
         write(lunw,'(a,i0)')'Error: bad resistor config index:',code
         write(*,'(a,i0)')'Error: bad resistor config index: ',code
c
c  If no error, produce header file.  We loop over the detector channels,
c  then over the resistors in each channel.
      else
         write(lunw,'(a)')'/* resistors.h Defines resistor strings */'
         write(lunw,'(a)')'#ifndef RESISTORS_H_INCLUDED'
         write(lunw,'(a)')'#define RESISTORS_H_INCLUDED'
         write(lunw,*)
         write(lunw,'(a,i0)')'#define R_Config_Index ',code
         write(lunw,*)

         do ichan=1,mch
            write(lunw,'(3a)')'static char *',
     &      detname(ichan)(1:lnbc(detname(ichan))),'_R[] = {'
            do ipgn=1,npgn
               write(lunw,'(a,i0,a)')'"',resistors(ipgn,ichan,code),'",'
            enddo
            write(lunw,'(a)')'};'
            write(lunw,*)
         enddo

         write(lunw,'(a)')'#endif'
      endif
      close(lunw)
      return
      end
