      subroutine read_opus_header(path,iend,dtype,nsp,fxv,lxv,iy,im,
     & id,hh,mm,ss,ms,apt,dur,vel,apf,phr,res,lwn,foc,nip,dfr,
     & pkl,prl,gfw,gbw,lfl,hfl,possp,oblat,oblon,obalt,
     & tins,pins,hins,tout,pout,hout,sia,sis)
c
c  Subroutine to extract information from OPUS headers.
c
c INPUTS:
c       path  C*(*)  Path to spectrum
c       iend  I*4  Endianess of host computer
c       dtype I*4  Data Type (spectrum=1031; interferogram=2055)

c OUTPUTS:
c       nsp,   I*4  ! Number of Spectral Points
c       fxv,   R*8  ! First X-value
c       lxv,   R*8  ! Last X-value
c       iy,    I*4  ! Year (4 digit)
c       im,    I*4  ! Month
c       id,    I*4  ! Day
c       hh,    I*4  ! Hour
c       mm,    I*4  ! Minute
c       ss,    I*4  ! Seconds
c       ms,    I*4  ! Milli-Seconds
c       apt,   R*8  ! Aperture diameter (mm)
c       vel,   R*8  ! Bruker "Velocity" (laser fringe rate in kHz)
c       apf,   C*2  ! Apodization function (e.g. "BX")
c       phr,   R*8  ! Phase Resolution (cm-1) = 0.9/OPD_short
c       res,   R*8  ! Spectral Resolution (cm-1) = 0.9/OPD_long
c       lwn,   R*8  ! Laser WaveNumber (cm-1 ) = 15798.
c       foc,   R*8  ! Focal Length of telescope/collimating mirror (mm)
c       nip,   I*4  ! Number of Interferogram Points
c       dfr,   I*4  ! Digital Filter Reduction
c       pkl,   I*4  ! Interferogram Peak (ZPD) location: forward scan
c       prl,   I*4  ! Interferogram Peak (ZPD) location: reverse scan
c       gfw,   I*4  ! Number of Good ForWard scans
c       gbw,   I*4  ! Number of Good BackWard scans
c       lfl,   R*8  ! Low Folding Limit (cm-1)
c       hfl,   R*8  ! High Folding Limit (cm-1)
c       possp  I*4  ! Position Of Starting Spectral Point (bytes)
c
      implicit none
      character path*(*),apf*2,cval*36,instname*7
      integer*4 dtype,luns,mc,iend,bigendian
      integer*2 mrs,rs
      parameter (luns=19,mrs=26,mc=26,bigendian=1)
      integer*2 i2val(mrs)

      integer*4 nsp,iy,im,id,hh,mm,ss,ms,nip,pkl,prl,possp,i4val,lc,lnbc
      real*8 fxv,lxv,vel,phr,res,apt,lwn,foc,dur,r8val,sia,sis,lfl,hfl,
     & velocity(10),oblat,oblon,obalt,tins,pins,hins,tout,pout,hout
      integer*4 i,ndb,mdb,magic,prog,ip,ivel,ldot,btype,bpointer,
     & itype(mc),ilen(mc),ipoint(mc),reclen,dfr,gfw,gbw

c----------------------------------------------------------
      equivalence (i2val,i4val,r8val,cval)
      data velocity/0.d0,0.d0,4.d0,5.d0,7.5d0,10.d0,20.d0,40.d0,
     & 60.d0,80.d0/

c
c  Read Header Block.
        reclen=12
c        write(*,*)'read_opus_header:',path
        open(luns,file=path,form='unformatted',status='old',
     &  access='direct',recl=reclen)
        read(luns,rec=1)magic,prog        ! magic number, version number
        if(iend.eq.bigendian) call rbyte(magic,4,1)
        if(iend.eq.bigendian) call rbyte(prog,8,1)
c        write(*,*)magic,prog
        read(luns,rec=2)ip,mdb,ndb          ! pointer, max. size, cur. size
        if(iend.eq.bigendian) then
           call rbyte(ip,4,1)
           call rbyte(mdb,4,1)
           call rbyte(ndb,4,1)
        endif
c        write(*,*)'ip, mdb, ndb =', ip,mdb,ndb
        if(ndb.gt.mc) then
           write(*,*)'ndb,mc=',ndb,mc
           stop 'ndb > mc'
        endif
c
        ip=ip/reclen  ! convert from bytes to records
c  Read NDB Directory Blocks.
        do i=1,ndb
           read(luns,rec=i+ip) itype(i),ilen(i),ipoint(i)
        end do
        close(luns)
        if(iend.eq.bigendian) then
           call rbyte(itype,4,ndb)
           call rbyte(ilen,4,ndb)
           call rbyte(ipoint,4,ndb)
        endif
c
c  Initialize variables that might be missing from the header
c  and therefore not set inside rdopushead.
        gfw=0    ! Added Oct 8 2007 GCT
        gbw=0    ! Added Oct 8 2007 GCT
        pkl=0
        prl=0
        nip=0
        dfr=0
        reclen=2
        open(luns,file=path,form='unformatted',status='old',
     &  access='direct',recl=reclen)
c
c  Extract relevent parameter values from parameter blocks.
c  This is a bit messy because different version of OPUS
c  have different parameters in different blocks and have
c  a different order of the blocks. Also,  the same parameter
c  names appear in multiple blocks with different values.
        do i=1,ndb
           btype=mod(itype(i),2**30)  ! block type
           bpointer=ipoint(i)         ! block pointer
c
c  1047 (417 hex) is DSTAT parameters for Master spectrum (Si)
c  1047+32768 (8417 hex) is DSTAT parameters for Slave spectrum (InGaAs)
          if(mod(btype,32768).eq.1047) then  !  AMPL+SAMP+DSTAT+SPEC 
             call getopusparval(luns,bpointer,'NPT',iend,mrs,i2val,rs)
             if(rs.eq.2) nsp=i4val
             call getopusparval(luns,bpointer,'FXV',iend,mrs,i2val,rs)
             if(rs.eq.4) fxv=r8val
             call getopusparval(luns,bpointer,'LXV',iend,mrs,i2val,rs)
             if(rs.eq.4) lxv=r8val
c             write(*,*)'i,NSP,FXV,LXV=',i,nsp,fxv,lxv
             call getopusparval(luns,bpointer,'DAT',iend,mrs,i2val,rs)
             if(rs.gt.0) then
                if(index(cval,'/').gt.4) then  !  YYYY/MM/DD
                   read(cval,'(i4,1x,i2,1x,i2)')iy,im,id
                else                           !  DD/MM/YYYY
                   read(cval,'(i2,1x,i2,1x,i4)')id,im,iy
                endif
             endif
c             write(*,*)'iy,im,id=',iy,im,id
             call getopusparval(luns,bpointer,'TIM',iend,mrs,i2val,rs)
             if(rs.gt.10 .and. cval(9:9).eq.'.') then
                 read(cval,'(i2,1x,i2,1x,i2,1x,i3)')hh,mm,ss,ms
             elseif(rs.ge.6) then
                 ms=0
                 read(cval,'(i2,1x,i2,1x,i2)')hh,mm,ss
             endif

          elseif(btype.eq.160) then  ! ORGPAR block
             call getopusparval(luns,bpointer,'LAT',iend,mrs,i2val,rs)
             if(rs.eq.4) oblat=r8val
             call getopusparval(luns,bpointer,'LON',iend,mrs,i2val,rs)
             if(rs.eq.4) oblon=r8val
             call getopusparval(luns,bpointer,'ALT',iend,mrs,i2val,rs)
             if(rs.eq.4) obalt=r8val/1000.
             call getopusparval(luns,bpointer,'TOU',iend,mrs,i2val,rs)
             if(rs.eq.4) tout=r8val
             call getopusparval(luns,bpointer,'POU',iend,mrs,i2val,rs)
             if(rs.eq.4) pout=r8val
             call getopusparval(luns,bpointer,'HOU',iend,mrs,i2val,rs)
             if(rs.eq.4) hout=r8val
             call getopusparval(luns,bpointer,'SIA',iend,mrs,i2val,rs)
             if(rs.eq.4) sia=r8val
             call getopusparval(luns,bpointer,'SIS',iend,mrs,i2val,rs)
             if(rs.eq.4) sis=r8val
          elseif(btype.eq.96) then  ! OPTPAR block
             call getopusparval(luns,bpointer,'APT',iend,mrs,i2val,rs)
             if(rs.ge.4) read(cval(:index(cval,'mm')-1),*) apt
             call getopusparval(luns,bpointer,'VEL',iend,mrs,i2val,rs)
             ldot=index(cval(:4),'.')
c             write(*,*)cval(:2*rs)
             if(ldot.gt.0) then  ! in the new Brukers 
                read(cval(:4),*) vel  ! read velocity directly
             else              ! old Brukers
                lc=lnbc(cval(:4))
                read(cval(:lc),*) ivel ! read velocity code value
                vel=velocity(ivel)
             endif
c            write(*,*) 'rs, VEL=',rs,' ',cval(:4),vel
          elseif(mod(btype,32768).eq.2071) then  ! AMPL+SAMP+DSTAT+IGRM
             call getopusparval(luns,bpointer,'NPT',iend,mrs,i2val,rs)
             if(rs.eq.2) nip=i4val  ! number of Interferogram Points
             call getopusparval(luns,bpointer,'DAT',iend,mrs,i2val,rs)
             if(rs.gt.0) then
                if(index(cval,'/').gt.4) then  !  YYYY/MM/DD
                   read(cval,'(i4,1x,i2,1x,i2)')iy,im,id
                else                           !  DD/MM/YYYY
                   read(cval,'(i2,1x,i2,1x,i4)')id,im,iy
                endif
             endif
             call getopusparval(luns,bpointer,'TIM',iend,mrs,i2val,rs)
             if(rs.gt.10 .and. cval(9:9).eq.'.') then
                 read(cval,'(i2,1x,i2,1x,i2,1x,i3)')hh,mm,ss,ms
             elseif(rs.ge.6) then
                 ms=0
                 read(cval,'(i2,1x,i2,1x,i2)')hh,mm,ss
             endif

          elseif(btype.eq.64) then  !  FTPAR block
             call getopusparval(luns,bpointer,'APF',iend,mrs,i2val,rs)
             if(rs.gt.0) apf=cval(:2)
             call getopusparval(luns,bpointer,'PHR',iend,mrs,i2val,rs)
             if(rs.eq.4) phr=r8val
          elseif(btype.eq.48) then  !  AQPAR block
             call getopusparval(luns,bpointer,'RES',iend,mrs,i2val,rs)
             if(rs.eq.4) res=r8val
c             call getopusparval(luns,bpointer,'NSS',iend,mrs,i2val,rs)
c             if(rs.eq.2) nss=i4val
          elseif(btype.eq.32) then  !  INSTR block
             call getopusparval(luns,bpointer,'LFL',iend,mrs,i2val,rs)
             if(rs.eq.4) lfl=r8val
             call getopusparval(luns,bpointer,'HFL',iend,mrs,i2val,rs)
             if(rs.eq.4) hfl=r8val
             call getopusparval(luns,bpointer,'DUR',iend,mrs,i2val,rs)
             if(rs.eq.4) dur=r8val
             call getopusparval(luns,bpointer,'LWN',iend,mrs,i2val,rs)
             if(rs.eq.4) lwn=r8val
             call getopusparval(luns,bpointer,'DFR',iend,mrs,i2val,rs)
             if(rs.eq.2) dfr=i4val
             call getopusparval(luns,bpointer,'HUM',iend,mrs,i2val,rs)
             if(rs.eq.4) hins=r8val
             if(hins.gt.99.9) hins=99.9
             call getopusparval(luns,bpointer,'PIM',iend,mrs,i2val,rs)
             if(rs.eq.4) pins=r8val
             call getopusparval(luns,bpointer,'TSC',iend,mrs,i2val,rs)
             if(rs.eq.4) tins=r8val
             call getopusparval(luns,bpointer,'PKL',iend,mrs,i2val,rs)
             if(rs.eq.2) pkl=i4val
             call getopusparval(luns,bpointer,'PRL',iend,mrs,i2val,rs)
             if(rs.eq.2) prl=i4val
             call getopusparval(luns,bpointer,'GFW',iend,mrs,i2val,rs)
             if(rs.eq.2) gfw=i4val
             call getopusparval(luns,bpointer,'GBW',iend,mrs,i2val,rs)
             if(rs.eq.2) gbw=i4val
             call getopusparval(luns,bpointer,'INS',iend,mrs,i2val,rs)
             if(rs.gt.0) instname=cval(:7)
             call getopusparval(luns,bpointer,'FOC',iend,mrs,i2val,rs)
             if(rs.gt.0) then
                foc=r8val
             else
                if( instname(1:6) .eq. 'IFS120' ) then
                   if( instname(7:7) .eq. 'M' ) then
                      foc=220.  ! Bruker 12XM
                   else
                      foc=418.  ! Bruker 12XHR.
                   endif
                endif
             endif
c  If btype=1031 (407 hex) it's the master spectrum data bloak
c  If btype=1031+32768 (8407 hex) it's the slave spectrum data block
c  If btype=2055 (807 hex) it's the master interferogram data bloak
c  If btype=2055+32768 (8807 hex) it's the slave interferogram data block
          elseif(mod(btype,32768).eq.dtype) then  ! spectrum/igram data
             possp=bpointer
          endif
        end do  ! i=1,ndb
        close(luns)
        return
        end
