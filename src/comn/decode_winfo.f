      subroutine decode_winfo(winfo,mfp,
     & ntg,ncbf,nfp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,pars)

c  Inputs:
c     winfo (c*(*)  string defining window and fitted parameters
c     mfp    I*4    Max allowed number of fitted parameters
c
c  Outputs:
c     ntg    I*4    Number of Target gases
c     ncbf   I*4    Number of Continuum Basis Functions
c     nfp    I*4    Number of fitted parameters (ntg+ncbf+..)
c     iptg   I*4    Pointer to first target gas (usually =1)
c     ipcl   I*4    Pointer to first CBF (=iptg+ntg)
c     ipfs   I*4    Pointer to FS  (=ipcl+ncbf if fitted)
c     ipsg   I*4    Pointer to SG 
c     ipzo   I*4    Pointer to Zero Offset
c     ipcf   I*4    Pointer to Channel fringes (=0)
c     pars C(mfp)** Target gas names extracted from winfo

      integer*4 lnbc,lc,lw,isv,lbf,mfp,ntg,ncbf,nfp,
     & iptg,ipcl,ipfs,ipsg,ipzo,ipcf
      character winfo*(*),pars(mfp)*(*)

      iptg=0
      ipcl=0
      ipfs=0
      ipsg=0
      ipzo=0
      ipcf=0
      lc=index(winfo,':')
      call lowercase(winfo)
      call substr(winfo(lc+1:),pars,mfp,ntg)
c      if(ntg+nfp+3.gt.mfp) then
c         write(*,*)' gfit: Error: NTG > MTG ',ntg,mfp
c         stop 'Increase parameter MTG in gfit.f'
c      endif
      if(ntg.gt.0) iptg=1

      lw=lnbc(winfo)
      isv=ntg     !   ISV is index into the State Vector (SV)
      if(index(winfo(:lc),' ncbf=').gt.0) then
         lbf = index(winfo(:lc),' ncbf=')
         read(winfo(lbf+6:),*)ncbf
         if(ncbf+ntg.gt.mfp) then
            write(*,*) 'ncbf,ntg,mfp = ',ncbf,ntg,mfp
            stop 'ncbf+ntg > mfp'
         endif
         do j=1,ncbf
            write(pars(ntg+j),'(a3,i2.2,a1)') 'cbf',j,' '
         end do
      else
         ncbf=0
         if(index(winfo(:lc),' cl ').gt.0) then
            write(pars(ntg+1),'(a3,i2.2,a1)') 'cbf',1,' '
            ncbf=ncbf+1
         endif
         if(index(winfo(:lc),' ct ').gt.0) then
            write(pars(ntg+2),'(a3,i2.2,a1)') 'cbf',2,' '
            ncbf=ncbf+1
         endif
         if(index(winfo(:lc),' cc ').gt.0) then
            write(pars(ntg+3),'(a3,i2.2,a1)') 'cbf',3,' '
            ncbf=ncbf+1
         endif
      endif   ! index(winfo(:lc),' ncbf=').gt.0
      if(ncbf.gt.0) ipcl=isv+1
      isv=isv+ncbf
      if(index(winfo(:lc),' fs ').gt.0) then
         isv=isv+1
         if(isv.gt.mfp) stop 'fs: isv > mfp'
         ipfs=isv
         pars(isv)='fs  '
      endif
      if(index(winfo(:lc),' sg ').gt.0) then
         isv=isv+1
         if(isv.gt.mfp) stop 'sg: isv > mfp'
         ipsg=isv
         pars(isv)='sg  '
      endif
      if(index(winfo(:lc),' zo ').gt.0) then
         isv=isv+1
         if(isv.gt.mfp) stop 'zo: isv > mfp'
         ipzo=isv
         pars(isv)='zo  '
      endif
      nfp=isv
      if(nfp.gt.mfp) stop 'NFP > MFP  Increase parameter MFP'

      return
      end
