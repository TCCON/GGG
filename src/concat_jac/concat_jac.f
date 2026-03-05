c  Program concat_jac.f
c  Concatenates (in wavenumber) Jacobian files of individual windows
c  of a gas into a single pseudo-contiguous jac file covering all the
c  windows of the gas. This can then be fed to avg_ker.f to generate
c  a collective averaging kernel for each gas, for each spectrum.
c
c  This cannot be done from the kernels of the individual windows;
c  you have to go back to the Jacobians to get the correct weighting.
c  
c  Assumes equal weight for each spectral point in each individual
c  window, so wide windows with lots of target absorption lines have
c  more influence than narrow windows with a single absorption line.
c
c  Uses multiggg.sh as the driving input file. The Jac files of the
c  individual windows are opened and read in parallel, so there could
c  be alot of LUNs open simultaneously.
c
c  The main complication is that the retrieved parameters can differ
c  from window to window of the same gas. Even if the exact same gases
c  are fitted, but they are listed in a different order, it creates a
c  problem. So it is necessary to create a super-set of the retrieved
c  parameters and set the Jacobians to zero for any parameter in the
c  super-set that is NOT retrieved in a particular window. An index
c  is then created that maps the Jac's of each individual window onto
c  the superset, so that they are read into the correct row of the
c  concatenated Jac matrix. 
c
c  The ordering of the rows of the JAC matrix is as follows:
c     First NTG  rows are the super-set Target gas JACs
c     Next  NCBF rows are the super-set Continuum JACs
c     Next  row is FS  JACs (if fitted)
c     Next  row is SG  JACs (if fitted)
c     Next  row is ZO  JACs (if fitted)
c
c
c  Pseudo-code:
c
cc  Read multiggg.sh to create mapping (gwmap) between gases & windows.
cc  gwmap(igas) is the row in multiggg.sh of the last window of igas.
cc  So the windows of igas occupy rows windows(igas-1)+1 to windows(igas)
cc
cc  Loop over gases
c      kwin1=1
c      do igas=1,ngas
c         kwin2=gwmap(igas)
c         
cc  Create a super-set of the retrieved parameters.
cc  Loop over windows of a particular gas, reading just the first record (winfo).
c         do jwin=kwin1,kwin2
c            open(lunr_jac+jwin,file=gggdir//'jac/'//jac_infile(jwin),status='old') ! open Jac file
c            read(lunr_jac+jwin,'(a)') winfo(jwin)
c            call decode_winfo(winfo(jwin),mfp,ntg_sw,ncbf_sw,nfp_sw,
c     &      iptg,ipcl,ipfs,ipsg,ipzo,ipcf,svlabels_sw)
c            
c            do jfp=1,nfp_sw
c               svmap(jfp)=0
c               do i2=nfp_ss,1,-1
c                  if(svlabels_sw(jfp).eq.svlabels_ss(i2)) svmap(jfp)=i2
c               end do
c               if(svmap(jfp).eq.0) then  ! New entry for State Vector (SV)
cc                 Inset new entry into right place in svlabels_ss,
cc                 right-shifting entries already at and above this location.
c                  nfp_ss=nfp_ss+1
c               endif         !  if(svmap(jfp).eq.0) then  ! New entry to SV
c            end do        !   do jfp=1,nfp_sw
c         end do  !   do jwin=kwin1,kwin2
c
cc  Loop again over same windows, but this time reading the entire jac
cc  file and writing their collective Jacobians
c         open(lunw_jac,file='/home/toon/ddd/jac/'//jac_outfile,status='unknown')
c         write(lunw_jac,'(30a8)') (svlabels_ss(j),j=1,nfp_ss)
c         do ispec=1,mspec  ! Loop over spectra in the same Jac file.
c            nmp_ss=0
c            do jwin=kwin1,kwin2
c               call decode_winfo(winfo(jwin),mfp,ntg_sw,ncbf_sw,nfp_sw,
c     &         iptg,ipcl,ipfs,ipsg,ipzo,ipcf,svlabels_sw)
c               call clistindex(nfp_sw,svlabels_sw,nfp_ss,svlabels_ss,sw2ssmap)
c               read(lunr_jac+jwin,*) kmp,ktg,kfp,klev
c               do jtg=1,ntg_sw !  Target Gas Column Jacobians
c                  read(lunr_jac+jwin,*) (pd(nmp_ss+jmp,sw2ssmap(jtg)),
c     &            jmp=1,kmp)
c               end do
c               do ilev=1,klev  ! Single-Level Jacobians of 1st target gas
c                  read(lunr_jac+jwin,*)(slpd(nmp_ss+jmp,ilev),jmp=1,kmp)
c               end do
c               do jtg=ktg+1,kfp
c                  read(lunr_jac+jwin,*)(pd(nmp_ss+jmp,sw2ssmap(jtg)),
c     &            jmp=1,kmp)   ! Continuum, FS, SG, ZO Jacobians
c               end do
c               nmp_ss=nmp_ss+kmp
c               read(lunr_jac+jwin,*)(splosconc(ilev),ilev=1,klev)
c               read(lunr_jac+jwin,*) pout
c               read(lunr_jac+jwin,*)(z(ilev),ilev=1,klev)
c               read(lunr_jac+jwin,*)(p(ilev),ilev=1,klev)
c               read(lunr_jac+jwin,*)(ynoise(ifp),ifp=1,kfp)
c               read(lunr_jac+jwin,*)(aprxmcx(ifp),ifp=1,kfp)
c            end do  !   do jwin=kwin1,kwin2
c
cc  Write wavenumber-concatenated Jac file
c            write(lunw_jac,'(a)') specname
c            write(lunw_jac,*) zmin,asza,rmsocl
c            write(lunw_jac,*) nmp_ss,ntg_ss,nfp_ss,klev
c            do jtg=1,ntg_ss
c               write(lunw_jac,*) (pd(jmp,jtg),jmp=1,nmp_ss)
c            end do
c            do ilev=1,klev
c               write(lunw_jac,*)(slpd(jmp,ilev),jmp=1,nmp_ss)
c            end do
c            do jtg=ntg_ss+1,nfp_ss
c               write(lunw_jac,*)(pd(jmp,jtg),jmp=1,nmp_ss)   ! CL, CT, CC, FS, ZO Jacobians
c            end do
c            write(lunw_jac,*)(splosconc(ilev),ilev=1,klev)
c            write(lunw_jac,*) pout
c            write(lunw_jac,*)(z(ilev),ilev=1,klev)
c            write(lunw_jac,*)(p(ilev),ilev=1,klev)
c            write(lunw_jac,*)(ynoise(ifp),ifp=1,nfp_ss)
c            write(lunw_jac,*)(aprxmcx(ifp),ifp=1,nfp_ss)
c         end do     !   ispec=1,mspec
c         write(*,*) 'Warning: Increase parameter mspec = ', mspec
c99       nspec=ispec-1  ! hit end of Jac file
c         close(lunw_jac)
c         kwin1=kwin2+1
c      end do   !   do igas=1,ngas
c      stop
c      end

      implicit none
      integer*4 lunr_mggg,lunr_jac,lunw_jac,mgas,ngas,igas,
     & mlev,ilev,klev,mwin,nwin,jwin,kwin1,kwin2,kus,i2,
     & mfp,nfp_sw,nfp_ss,jfp,nmp_ss,mmp,jmp,ntg_sw,ntg_ss,jtg,
     & ktg,kfp,ifp,kmp,ksp,kgg,ls,fbc,j,jj,
     & l_,ldot,lrt,mspec,nspec,ispec,lnbc,
     & ipfs_ss,ipsg_ss,ipzo_ss,
     & iptg,ipcl,ipfs,ipsg,ipzo,ipcf,ncbf_sw,ncbf_ss

      parameter (lunr_mggg=17,lunr_jac=25,lunw_jac=25,mgas=18,
     & mmp=100000,mfp=30,mlev=150,mwin=200,mspec=1000000)
      integer*4 gwmap(mgas),svmap(mfp),sw2ssmap(mfp),index

      real*4 zmin,zminwas,asza,aszawas,
     & airmass,amwas,
     & rmsocl,trmsocl,pd(mmp,mfp),
     & pout,z(mlev),p(mlev),
     & ynoise(mfp),aprxmcx(mfp),
     & tynoise(mfp),taprxmcx(mfp),
     & splosconc(mlev),slpd(mmp,mlev)

      character string*128,jac_infile(mwin)*256,jac_outfile*256,
     & gasname*8,gaswas*8,version*60,
     & gggdir*256,dl*1,
     & specnamewas*64,specname*64,winfo(mwin)*128,
     & svlabels_sw(mfp)*8,svlabels_ss(mfp)*8

      version= 
     &' CONCAT_JAC.F             Version 0.02        2021-02-23 GCT'

      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)     !Length of gggdir

c  Read multiggg.sh to create mapping (gwmap) between gases & windows.
c  gwmap(igas) is the row in multiggg.sh of the last window of igas.
c  So the windows of igas occupy rows windows(igas-1)+1 to windows(igas)
      gaswas='        '
      igas=0
      open(lunr_mggg,file='multiggg.sh',status='old')
      do jwin=1,mwin     ! Loop over windows in multiggg file
         read(lunr_mggg,'(a)',end=88) string
         ksp=fbc(string)
c  By indexing the part of the string after the space and adding the
c  space index back, we handle cases where "_" or ".ggg>" appear in the
c  path to gfit
         kus=index(string(ksp:),'_')+ksp-1
         kgg=index(string(ksp:),'.ggg>')+ksp-1
         gasname=string(ksp+1:kus-1)
         jac_infile(jwin)='j_'//string(ksp+1:kgg-1)
         ls=kgg-ksp+1
         if(gasname.ne.gaswas) then
            igas=igas+1
         endif
         write(*,*) igas,jwin,ksp,kus,kgg,ls,jac_infile(jwin)(:ls)
         gwmap(igas)=jwin
         gaswas=gasname
      enddo ! jwin=1,mwin
      write(*,*)'Warning: Increase parameter mwin = ',mwin
88    ngas=igas
      close(lunr_mggg)
      nwin=jwin-1
 
c      do igas=1,ngas
c         write(*,*) 'igas,gwmap(igas) = ',igas,gwmap(igas)
c      end do

c  Loop over gases
      kwin1=1
      do igas=1,ngas
         kwin2=gwmap(igas)
         
c  Create a super-set of the retrieved parameters.
         ntg_ss=0
         nfp_ss=0
         ipfs_ss=0
         ipsg_ss=0
         ipzo_ss=0
         ncbf_ss=0
         nmp_ss=0
c         write(*,*)'igas,kwin1,kwin2 = ',igas,kwin1,kwin2
c  Loop over windows of a particular gas reading just the first record (winfo).
         do jwin=kwin1,kwin2
            open(lunr_jac+jwin,file=gggdir(:lrt)//dl//'jac'//dl//
     &      jac_infile(jwin),status='old') ! open .jac file
            read(lunr_jac+jwin,*)   !  Version info
            read(lunr_jac+jwin,'(a)') winfo(jwin)
            write(*,*)jwin,'  winfo = '//winfo(jwin)
            call decode_winfo(winfo(jwin),mfp,ntg_sw,ncbf_sw,nfp_sw,
     &      iptg,ipcl,ipfs,ipsg,ipzo,ipcf,svlabels_sw)
            
            do jfp=1,nfp_sw
               svmap(jfp)=0
               do i2=nfp_ss,1,-1
                  if(svlabels_sw(jfp).eq.svlabels_ss(i2)) svmap(jfp)=i2
               end do
               if(svmap(jfp).eq.0) then  ! New entry for State Vector (SV)
c               write(*,*) 'new entry to SV',jfp,svlabels_sw(jfp)
c   Find the correct place in svlabels_ss to insert the new SV entry,
c   right-shifting the entries already at and above this location.
                  if(jfp.le.ntg_sw) then   ! New target gas entry for state vector
                     do jj=nfp_ss,ntg_ss+1,-1
                        svlabels_ss(jj+1)=svlabels_ss(jj)
                     end do
                     svlabels_ss(ntg_ss+1)=svlabels_sw(jfp)
                     ntg_ss=ntg_ss+1
                  elseif(jfp.le.ntg_sw+ncbf_sw) then  ! New CBF entry for state vector
                     do jj=nfp_ss,ntg_ss+ncbf_ss+1,-1
                        svlabels_ss(jj+1)=svlabels_ss(jj)
                     end do
                     svlabels_ss(ntg_ss+ncbf_ss+1)=svlabels_sw(jfp)
                     ncbf_ss=ncbf_ss+1
                  else                ! New FS SG or ZO entry for state vector
                     svlabels_ss(nfp_ss+1)=svlabels_sw(jfp)
                  endif
                  nfp_ss=nfp_ss+1
               endif         !  if(svmap(jfp).eq.0) then  ! New entry to SV
            end do        !   do jfp=1,nfp_sw
 
         end do  !   do jwin=kwin1,kwin2
         write(*,*)'svlabels_ss = ',(svlabels_ss(j),j=1,nfp_ss)

c  Loop again over same windows, but this time reading the entire jac
c  files and writing their collective Jacobians
         l_=index(jac_infile(kwin1)(3:),'_')+2
         ldot=index(jac_infile(kwin1)(l_+1:),'.')+l_
         jac_outfile=jac_infile(kwin1)(1:l_)//
     &               jac_infile(kwin1)(ldot+1:)
         open(lunw_jac,file=gggdir(:lrt)//dl//'jac'//dl//jac_outfile,
     &   status='unknown')
         write(lunw_jac,'(a)') version
         write(lunw_jac,'(30a8)') (svlabels_ss(j),j=1,nfp_ss)
         do ispec=1,mspec  ! Loop over spectra in the same Jac file.
            nmp_ss=0
            trmsocl=0.0
            tynoise=0.0
            taprxmcx=0.0
            do jwin=kwin1,kwin2
               call decode_winfo(winfo(jwin),mfp,ntg_sw,ncbf_sw,nfp_sw,
     &         iptg,ipcl,ipfs,ipsg,ipzo,ipcf,svlabels_sw)
c               if(ispec.eq.1) then
c                write(*,*)'svlabels_sw=',(svlabels_sw(jfp),jfp=1,nfp_sw)
c               endif
               sw2ssmap=0
               call clistindex(nfp_sw,svlabels_sw,nfp_ss,svlabels_ss,
     &         sw2ssmap)
c               if(ispec.eq.1) then
c               do jfp=1,nfp_sw
c                  write(*,*) jfp,svlabels_sw(jfp),sw2ssmap(jfp),
c     &            svlabels_ss(sw2ssmap(jfp))
c               end do
c               endif
c  Check that the different windows have the same spectrum names,
c  and zmin/asza values
               read(lunr_jac+jwin,'(a)',end=99) specname
               read(lunr_jac+jwin,*) zmin,asza,rmsocl,airmass
               trmsocl=trmsocl+rmsocl
               if(jwin.eq.kwin1) then
                  specnamewas=specname
                  zminwas=zmin
                  aszawas=asza
                  amwas=airmass
               else
                  if(specname(1:15).ne.specnamewas(1:15) .or.
     &            specname(17:).ne.specnamewas(17:)) then
                     write(*,*)'       jwin,       kwin1,      specname(
     &jwin),       specnamewas(kwin1)'
                     write(*,*) jwin,kwin1,
     &               '  '//specname(:24)//'  '//specnamewas(:24)
                     stop 'mismatched spectra'
                  endif
                  if(abs(zmin-zminwas).gt.0.001) stop 'zmin mismatch'
                  if(abs(asza-aszawas).gt.0.001) stop 'asza mismatch'
                  if(abs(airmass-amwas).gt.0.001)stop 'airmass mismatch'
               endif
               read(lunr_jac+jwin,*) kmp,ktg,kfp,klev
               if(nmp_ss+kmp.gt.mmp) stop 'nmp_ss > mmp'
               if(kfp.gt.mfp) stop 'kfp > mfp'
               if(ntg_sw .ne. ktg) stop 'ktg mismatch'
               if(nfp_sw .ne. kfp) stop 'kfp mismatch'
               do jtg=1,ntg_ss
                  do jmp=1,kmp
                     pd(nmp_ss+jmp,jtg)=0.0
                  end do
               end do
               do jtg=1,ntg_sw !  Target Gas Column Jacobians
                  read(lunr_jac+jwin,*) (pd(nmp_ss+jmp,sw2ssmap(jtg)),
     &            jmp=1,kmp)
               end do
               do ilev=1,klev  ! Single-Level Jacobians of 1st target gas
                  read(lunr_jac+jwin,*)(slpd(nmp_ss+jmp,ilev),jmp=1,kmp)
               end do
               do jtg=ntg_ss+1,nfp_ss
                  do jmp=1,kmp
                     pd(nmp_ss+jmp,jtg)=0.0
                  end do
               end do
               do jtg=ktg+1,kfp
                  read(lunr_jac+jwin,*)(pd(nmp_ss+jmp,sw2ssmap(jtg)),
     &            jmp=1,kmp)   ! Continuum, FS, SG, ZO Jacobians
               end do
               nmp_ss=nmp_ss+kmp
               read(lunr_jac+jwin,*)(splosconc(ilev),ilev=1,klev)
               read(lunr_jac+jwin,*) pout
               read(lunr_jac+jwin,*)(z(ilev),ilev=1,klev)
               read(lunr_jac+jwin,*)(p(ilev),ilev=1,klev)
               ynoise=0.0
               aprxmcx=0.0
               read(lunr_jac+jwin,*)(ynoise(sw2ssmap(ifp)),ifp=1,kfp)
               read(lunr_jac+jwin,*)(aprxmcx(sw2ssmap(ifp)),ifp=1,kfp)

c  Add together a priori cobstraints from individual windows.
               do ifp=1,kfp
                  tynoise(sw2ssmap(ifp)) =tynoise(sw2ssmap(ifp)) +
     &            ynoise(sw2ssmap(ifp))
                  taprxmcx(sw2ssmap(ifp))=taprxmcx(sw2ssmap(ifp))+
     &            aprxmcx(sw2ssmap(ifp))
               end do

            end do  !   do jwin=kwin1,kwin2
            write(*,*) ispec,jwin,'  '//specname

c  Write wavenumber-concatenated Jac file
            write(lunw_jac,'(a)') specname
            write(lunw_jac,*) zmin,asza,trmsocl/(kwin2-kwin1+1),airmass
            write(lunw_jac,*) nmp_ss,ntg_ss,nfp_ss,klev
            do jtg=1,ntg_ss
               write(lunw_jac,*) (pd(jmp,jtg),jmp=1,nmp_ss)
            end do
            do ilev=1,klev
               write(lunw_jac,*)(slpd(jmp,ilev),jmp=1,nmp_ss)
            end do
            do jtg=ntg_ss+1,nfp_ss
               write(lunw_jac,*)(pd(jmp,jtg),jmp=1,nmp_ss)   ! CL, CT, CC, FS, ZO Jacobians
            end do
            write(lunw_jac,*)(splosconc(ilev),ilev=1,klev)
            write(lunw_jac,*) pout
            write(lunw_jac,*)(z(ilev),ilev=1,klev)
            write(lunw_jac,*)(p(ilev),ilev=1,klev)
            write(lunw_jac,*)(tynoise(ifp),ifp=1,nfp_ss)
            write(lunw_jac,*)(taprxmcx(ifp),ifp=1,nfp_ss)
         end do     !   ispec=1,mspec
         write(*,*) 'Warning: Increase parameter mspec = ', mspec
99       nspec=ispec-1  ! hit end of Jac file
         do jwin=kwin1,kwin2
            close(lunr_jac+jwin)
         end do
         close(lunw_jac)
         write(*,*)
         
         kwin1=kwin2+1
      end do   !   do igas=1,ngas

      stop
      end
