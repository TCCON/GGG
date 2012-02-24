c  Program computes zenith angle offsets from the CO2 VSF factors.
      integer*4 j,k,mgas,ngas,kgas,jgas,luns,lunq,
     $ mobs,iobs,naux,ncol,lnbc,lr,ncomm
      parameter (luns=13)
      parameter (lunq=14)
      parameter (mobs=1500)
      parameter (mgas=120)
      parameter (naux=19)
      character pabel*800,tavfile*12,clabel(naux+2*mgas)*32
      real*8 yval(naux),totcon,toterr,pobs,asza,del,dsbydt
c
      write(6,*)'ZENANG   Version 2.1.4   19 Jul 2011   GCT'
 1    write(6,101)
 101  format('Enter name of .tav file (e.g.  fts93rat.tav, etc) ',$)
      read(*,85) tavfile 
 85   format(a)
      lr=lnbc(tavfile)
      if(lr.le.0) go to 1
c=====================================================================
c  Read in TOTCONS from SUNRUN.TAV
      open(luns,file=tavfile(:lr-3)//'tav',status='old')
      read(luns,*) ncomm
      do j=3,ncomm       !  Skip comment lines
        read(luns,*)
      end do
      read(luns,'(a)') pabel
      call substr(pabel,clabel,naux+2*mgas,ncol)
      ngas=(ncol-naux)/2
      do kgas=1,ngas
        if(index(clabel(naux+2*kgas-1),'co2').eq.1)goto 301
      end do
      stop 'Did not find CO2'
 301  continue
c
      open(lunq,file='zenang.out',status='unknown')
      do iobs=1,mobs
        read(luns,'(f14.8,22(f13.5),240(e12.4))',end=14)
     &  (yval(k),k=1,naux),(totcon,toterr,jgas=1,kgas)
        pobs=yval(18)
        asza=yval(9) 
        del=(totcon-1)/dsbydt(asza)/(1+toterr**2)
        write(lunq,'(5f14.8)') yval(1),totcon,pobs,asza,del
      end do  ! iobs=1,mobs
      iobs=iobs+1
      read(luns,*,end=14)
      write(6,*) 'Warning: increase parameter MOBS'
14    close(luns)
      close(lunq)
      stop
      end
