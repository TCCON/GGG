      subroutine generate_qc_flag(nrow,spectrum,qc)

      implicit none

      integer*4 nrow                ! # spectra (runlog)    intent in
      real*4 qc(nrow)               ! qc flag            intent inout
      character spectrum(nrow)*38   ! runlog spectrum names intent in

      integer i,i1,idot,icol,jtg,k,ktg,mcol,lunm,lunc,nit,nss,
     & mit, l2,l3,l4, fnbc, fbc, nlhead, date(nrow), selday,
     & istart, istop, j, jj  
      parameter (lunm=11)      ! multiggg.sh
      parameter (lunc=14)      ! .col file
      parameter (mcol=50)      ! max number of microwindows  

      character colabel*500,
     & gfit_version*80,gsetup_version*80,col_string*500,
     & csformat*90, cdum*20,
     & specname_col*38,runlog*80,tabel*80,
     & colfile*40,window(mcol)*10

      real*8 airmass,cl,tilt,zlo,fcen,width,fqshift,
     & ovcol,rmsfit,sg,zmin, cl_arr(nrow)

      real*4
     & vsf,vsf_err, cl_threshold, cl_limit, cl_max
      parameter(cl_threshold=0.9)

c initialise qc flag 

      do i=1, nrow
         qc(i)=0.0
      enddo
 
c parse o2_7885 col file and derive qc information

c Find the o2_7885 col file in multiggg.sh 
      open(lunm,file='multiggg.sh',status='old')
      do icol=1,mcol     !  main loop (over windows)
135     read(lunm,'(a)') tabel
        if(tabel(1:1).eq.':') go to 135
        colfile=tabel(index(tabel,'<')+1:index(tabel,'.ggg'))//'col'
        idot=index(colfile,'.')
        i1=index(colfile(:idot),'_')
        if(i1.eq.0) i1=index(colfile(:idot),'^')   !  new format
        if(i1.eq.0) i1=index(colfile(:idot),'-')   !  old format
        window(icol)=colfile(:idot-1)
        if (index(window(icol),"o2_7885") .gt. 0) then
           exit 
        endif
      enddo
      close(lunm)

      write(*,*) 'performing QC based on ', colfile
      open(lunc,file=colfile,status='old') ! .col file
c
c  Read header lines of .col file and locate column containing "OVC_gas".
c  in order to read data from appropriate target gas.
      read(lunc,*) nlhead
      read(lunc,'(a)') gfit_version
      read(lunc,'(a)') gsetup_version
      do k=4,nlhead-1
        read(lunc,'(a)')colabel
        if(index(colabel,'runlogs').gt.0) runlog=colabel(:80)
      end do
      read(colabel,*) fcen, width, mit
      read(lunc,'(a)')colabel
      ktg=1+index(colabel,' OVC_'//colfile(:i1-1))
      if ( ktg .gt. 1) then
          call substr(colabel(:ktg-1),cdum,1,nss)
          ktg=(nss-4)/4
      endif

      do i=1, nrow
           specname_col='='
           read(lunc,'(a)', end=25) col_string
           l2=fbc(col_string(2:))+1       ! First space following spectrum name
           l3=fnbc(col_string(l2:))+l2-1  ! First character of NIT
           l4=fbc(col_string(l3:))+l3-1   ! First space following NIT

           if (index(gfit_version,'2.40.2') .ne. 0) then !  old col file format 
               csformat='(1x,a21,i2,1x,f5.3,3(1x,f4.1),1x,f5.3,'
     &         //'1x,f6.4,f7.3,1x,9(0pf7.3,1pe10.3,0pf9.4,1pe8.1))'
           else                                  ! assume new col file format 
               write(csformat,'(a,i2.2,a)')'(1x,a',l4-5,
     &     ',i3,f6.3,3f5.1,f6.3,f7.4,f8.3,15(f7.3,e11.4,f9.4,e8.1))'
           endif

           read(col_string,csformat) specname_col,nit,cl,tilt,fqshift,
     &     sg,zlo,rmsfit,zmin,(airmass,ovcol,vsf,vsf_err,jtg=1,ktg)

           if (index(specname_col,spectrum(i)) .eq. 0) then
              write(*,*) "qc code ordering assumption violated"
              write(*,*) "Terminating processing"
              stop 
           endif

           cl_arr(i)=cl

c get the date from the filename
           read(specname_col,'(a2,i8,a)'), cdum, date(i), cdum

c converged retrievals
           if (nit .lt. mit) qc(i)=qc(i)+1. 
c O2 vf criteria 
           if (vsf .lt. 0.9 .or. vsf .gt. 1.1) qc(i)=0.
c O2 rms criteria 
           if (date(i) .lt. 20060630 .and. rmsfit .lt. 0.5) then 
              qc(i)=qc(i)+0.1
           endif
           if (date(i) .ge. 20060630 .and. rmsfit .lt. 1.0) then
              qc(i)=qc(i)+0.1
           endif
c shift criteria
           if (ABS(fqshift) .lt. 5. .and. ABS(sg) .lt. 3.) then 
              qc(i)=qc(i)+0.01 
           endif

      enddo   ! loop on i=1, nrow
      close(lunc)

c derive qc information dependent on data characteristics for day
c assumes records from the same day are contiguous in col/runlog
c cl criteria
      selday=99999999
      istart=1
      write(*,*) "deriving global daily qc info for dates:"
      do i=1, nrow
           if (date(i) .ne. selday .or. i .eq. nrow) then
              if (i .ne. nrow) istop=i-1
              if (i .eq. nrow) istop=i
              cl_max=0.0d0
              do jj=istart,istop
                 if(cl_arr(jj).gt.cl_max) cl_max=cl_arr(jj)
              end do
              cl_limit=cl_threshold*cl_max
c              cl_limit=cl_threshold*MAXVAL(cl_arr(istart:istop))
              do j=istart, istop
                 if (cl_arr(j) .gt. cl_limit) then
                    qc(j)=qc(j)+0.001
                 endif
              enddo
              selday=date(i) 
              istart=i
              if (i .ne. nrow) write(*,*) selday
           endif
      enddo   

c translate the qc value to a simple qc index
c 1 = S0T0
c 2 = S2T2 without 3-sigma o2 vf criteria
      do i=1, nrow
         if (qc(i) .lt. 1.1) qc(i)=0.
         if (qc(i) .ge. 1.1 .and. qc(i) .lt. 1.111) qc(i)=1.
         if (qc(i) .ge. 1.111) qc(i)=2.
      enddo

25    if (i-1 .ne. nrow) then 
          write(*,*) "End of file reached after ", i-1, " records"
          write(*,*) nrow, "records expected from runlog"
          write(*,*) "Terminating processing"
          stop 
      endif

      return 
      end
