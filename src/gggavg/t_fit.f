c  Subroutine T_FIT
c  Takes a 2-D matrix Y(i,j) of observations (e.g. slant columns from a
c  series of different spectra and windows) and fits each point as the
c  product of two piece-wise functions: (i) the average row values v(i),
c  and (ii) the the column scale factors s(j).
c  It does this by minimizing the quantity
c      SUM_ij([Y(i,j)-v(i).s(j)]/Y_err(i,j))**2 + SUM_j(S(j)-1)**2
c  with respect to each element of v(i) and s(j). Thus, i+j quantities
c  are derived from the i.j observations. Therefore, calling this
c  subroutine with either i.le.1 or j.le.1 is a do-nothing operation.
c
c  The first term represents our desire to fit the measurements.
c  The second term represents out desire to keep S(j) close to 1.0,
c  unless the measurements say otherwise. Essentially, S(j)=1 +/-1
c  is our a priori estimate.
c
c  Differentiating the penalty function equation wrt v(i) yields
c    SUM_j(-2.s(j)*[Y(i,j)-v(i).s(j)]/Y_err(i,j)**2 = 0  for all i
c  which can be reorganized as
c    v(i) = SUM_j(s(j)*Y(i,j)/Y_err(i,j)**2) / SUM_j(s(j)/Y_err(i,j))**2
c
c  Differentiating the penalty function equation wrt s(j) yields
c    SUM_i(-2.v(i)*[Y(i,j)-v(i).s(j)]/Y_err(i,j)**2 + 2.(s(j)-1) = 0  for all j
c  which can be reorganized as
c    s(j) = [1+SUM_i(v(i)*Y(i,j)/Y_err(i,j)**2)]/[1+SUM_i(v(i)/Y_err(i,j))**2]
c
c  Since the equations for v(i) depends on s(j), and the equations for s(j)
c  depend on v(i), these sets of simultaneous equations are solved iteratively


      subroutine t_fit(valmiss,nspe,mwin,nwin,yobs,yerr,
     $ybar,eybar,scale,escale,sew,error_weight)
c
c Inputs: 
c    valmiss           R*8  Missing Value (-999.0)
c    nspe              I*4  Number of i-values (i.e. spectra)
c    mwin              I*4  Declared dimension of j-values)
c    nwin              I*4  Actual number of windows
c    yobs(nspe,nwin)   R*8  Data values (e.g. slant columns)
c    yerr(nspe,nwin)   R*8  Uncertainties (e.g. slant column errors)
c
c  Outputs:
c    ybar(nspe)        R*8  Mean y-value (averaged over windows) for each spectrum
c    eybar(nspe)       R*8  Uncertainty in ybar
c    scale(nwin)       R*8  Multiplicative scale factor
c    escale(nwin)      R*8  Uncertainty in Multiplicative scale factor
c    sew
c    error_weight
c    
      implicit none
      integer nspe,nwin,mwin,iwin,jspe,loop
      real*8 yobs(mwin,nspe),yerr(mwin,nspe),eps,ss,vv,se,ve,
     $ybar(nspe),eybar(nspe),scale(mwin),escale(mwin),
     $sew(mwin),error_weight,total_ybar,total_scale,vnew,snew,res,
     $valmiss,ybar_num,ybar_denom,scale_num,scale_denom,grr,gwt,twt,trr
      parameter (eps=1.0d-5)
c
      do iwin=1,nwin
        scale(iwin)=1.0d0
      end do
c
      do jspe=1,nspe
        ybar(jspe)=1.0d0
      end do
c
      do loop=1,49
c
c  Re-evaluate ybars (or column abundances)
         vv=0.0d0
         total_ybar=0.0d0
         do jspe=1,nspe       ! loop over spectra
            ybar_num=0.d0
            ybar_denom=0.d0
            do iwin=1,nwin   !  loop over windows
               if( yerr(iwin,jspe).ne.valmiss .and.
     &         yerr(iwin,jspe) .ne. 0.0 ) then
                  se=scale(iwin)/yerr(iwin,jspe)
                  ybar_num=ybar_num+se*yobs(iwin,jspe)/yerr(iwin,jspe)
                  ybar_denom=ybar_denom+se*se
               endif
            end do
            if(ybar_denom.eq.0.0d0) then
               eybar(jspe)=valmiss
               vnew=valmiss
             else 
               eybar(jspe)=1.d0/dsqrt(ybar_denom)
               vnew=ybar_num/ybar_denom
            endif
            vv=vv+abs(vnew-ybar(jspe))
            ybar(jspe)=vnew
            total_ybar=total_ybar+abs(ybar(jspe))
         end do   !  jspe=1,nspe   end of ybar re-evaluation
c
c  Re-evaluate scale factors (for each window)
         ss=0.0d0
         total_scale=0.0d0
         do iwin=1,nwin                   ! loop over windows
            scale_num=1.0d0    ! A priori estimate
            scale_denom=1.0d0  ! A priori estimate
            do jspe=1,nspe        !  loop over spectra
               if(yerr(iwin,jspe).ne.valmiss .and.
     &         yerr(iwin,jspe).ne.0.0) then
                  ve=ybar(jspe)/yerr(iwin,jspe)
                  scale_num=scale_num+ve*yobs(iwin,jspe)/yerr(iwin,jspe)
                  scale_denom=scale_denom+ve*ve
               endif
            end do
            if(scale_denom.eq.0.0d0) then
               write(6,*)' No spectra for window'
               escale(iwin)=0.0d0
               snew=1.0d0
            else
               escale(iwin)=1.d0/dsqrt(scale_denom)
               snew=scale_num/scale_denom
c      write(*,*)loop,scale_num,scale_denom
            endif
            ss=ss+abs(snew-scale(iwin))
            scale(iwin)=snew
            total_scale=total_scale+abs(scale(iwin))
         end do    ! iwin=1,nwin    end of scale re-evaluation
c
c      write(6,*)loop,vv/total_ybar,ss/total_scale,scale(1)
      if(ss.le.eps*total_scale .and. vv.le.eps*total_ybar) go to 88
      end do                           !end of main iteration
      write(6,*)'Warning: t_fit did not converge'
c
c  Error_weight is the factor by which the given standard errors (STERR)
c  must be scaled to be consistent with the observed scatter in the points.
c  In a perfect world it would be unity.
88    grr=0.0d0
      gwt=0.0d0
      do iwin=1,nwin                   !calculate the rms fit
        trr=0.0d0
        twt=0.0d0
        do jspe=1,nspe
          if( yerr(iwin,jspe).ne.valmiss .and.
     &    yerr(iwin,jspe).ne.0.0) then
          twt=twt+1.0d0
          res=yobs(iwin,jspe)-ybar(jspe)*scale(iwin)
          trr=trr+( res/yerr(iwin,jspe) )**2
          endif
        end do
        sew(iwin)=dsqrt(trr/twt)   ! Scale Error Weight
        grr=grr+trr
        gwt=gwt+twt
      end do                         !end of rms calculation
      error_weight=dsqrt(grr/gwt)
      return
      end
