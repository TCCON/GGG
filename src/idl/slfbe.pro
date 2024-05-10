pro slfbe,nobs,x,ux,y,uy,a,ea,b,eb,chi2on,pcc

;  Fits a straight line of the form
;                   y = A + B.x
;  to the data in vectors x and y, which have uncertainties ux and uy
;
; INPUTS:
;    NMP       I*4  The number of measured points/rows
;    X(NMP)    R*4  Vector of X-Values
;    UX(NMP)   R*4  Vector of X-Errors
;    Y(NMP)    R*4  Vector of Y-Values
;    UY(NMP)   R*4  Vector of Y-Errors
;
; OUTPUTS:
;    A         R*4  Y-Offset of fitted line (Y-value at x=0)
;    UA        R*4  Uncertainty in OFF
;    B         R*4  Gradient of fitted line
;    UB        R*4  Uncertainty in B
;    CHI2ON    R*4  Factor by which deviations of the data points from
;                   the fitted line exceed their stated uncertainties.
;    PCC       R*4  Pearson Correlation Coefficient
;
;  The fitted line is determined in four distinct steps:
;
;  1) Determine the Chi2 for any line (y=a+b.x)
;
;  This is done by determining the closest approach of the line
;  to each point,  by minimizing the Chi2 of each point [x(i),y(i)]
;  from the line  y = a + b.x
;
;  Assuming that ux and uy are uncorrelated, the contours of
;  Chi2 around the point i form the ellipse:
;     Chi2(i)  = [x(i)-x]^2/ux(i)^2 + [y(i)-y]^2/uy(i)^2           (1)
;
;  Along the line y=a+b.x  Chi2(i) varies as
;     Chi2(i)  = [x(i)-x]^2/ux(i)^2 + [y(i)-a-b.x]^2/uy(i)^2       (2)
;  and will have a minimum at 'closest approach'
;
;  Differentiate equation (2) wrt x:
;    dChi2(i)/dx  = -2(x(i)-x)/ux(i)^2 - 2b(y(i)-a-b.x)/uy(i)^2
;
;  where this is zero, set x=x'
;    [x(i)-x']/ux(i)^2 = -b[y(i)-a-b.x']/uy(i)^2
;    x'[1/ux(i)^2+b^2/uy(i)^2] = x(i)/ux(i)^2+b(y(i)-a)/uy(i)^2
;    x' = [x(i)/ux(i)^2+b(y(i)-a)/uy(i)^2] / [1/ux(i)^2+b^2/uy(i)^2]  (3)
;
;  (x',y') = (x',a+bx') is the point on the line that is closest
;  to the measured point (x(i),y(i)) in terms of the Chi2 ellipse.
;  This is not the physically closest point (i.e., perpendicular
;  to line) unless ux=uy, when the ellipse becomes a circle.
;
;  Note that x' is a weighted average of  x(i) and (y(i)-a)/b
;  with weights of  1/ux(i)^2 and b^2/uy(i)^2, respectively.
;
;  Substituting x=x' into the Chi2 equations yields:
;    Chi2(i)  = (x(i)-x')^2/ux(i)^2 + (y(i)-a-b.x')^2/uy(i)^2
;    Chi2(i)  =  (y(i)-a-b.x(i))^2/(uy(i)^2+b^2.ux(i)^2)        (4)
;   
;  Sum the CHi2(i) over all points i
;    Chi2  =  SUM {(y(i)-a-b.x(i))^2/(uy(i)^2+b^2.ux(i)^2)}     (5)

;  In the case of ux(i)=0, equation (5) simplifies to the more
;  common equation for a least squares straight line fit.
;
;  Note that equation (5) is symmetrical wrt x and y, in the 
;  sense that switching x and y, ux and uy, and substituting
;  1/b for b  and a/b for a, results in the exact same equation
;
; 2) Determine A and B that minimizes Chi2
;
;  Having derived out the Chi2 function, minimize it wrt a and b: 
;     dChi2/da = SUM {-2(y(i)-a-b.x(i)).w(i)}          (6)
;     dChi2/db = SUM {-2[y(i)-a-b.x(i)].[(y(i)-a).b.ux(i)^2+x(i).uy(i)^2].w(i)^2}  (7)
;  where  w(i) = 1/[uy(i)**2+b^2*ux(i)^2]
;  The presence of b in the denominator makes the differentials
;  non-linear in b, necessitating an iterative approach.
;  Unless ux=0, in which case the b in the denominator doesn't
;  matter and so both differentials are linear in a and b.
;
;  Setting  dChi2/da = 0
;       a.SUM {w(i)} = SUM {[y(i)-b.x(i)].w(i)}
;       a.tw + b.twx = twy                                               (8)
;   So a is the weighted average value of (y-b.x), with weights of w
;   Note that a depends on b, which we don't know.
;
;  Setting  dChi2/db = 0
;       b.SUM {x(i).v(i) } = Sum {[y(i)-a].v(i)}
;       b.tvx + a.tv = tvy                                               (9)
;  So b is the weighted average of (y(i)-a)/x(i)  with weights
;     v(i) = [(y(i)-a).b.ux(i)^2+x(i).uy(i)^2].w(i)^2
;
;  So we have the matrix equation:
;     |tw twx| |a| = |twy|                                              (10)
;     |tv tvx| |b| = |tvy|
;
;  
;  Numerical Recipes claim that it is better to retrieve the angle
;  of the fitted line, c, rather than the gradient, b=tan(c)
;     dChi2/db = SUM {-2[y(i)-a-b.x(i)].[(y(i)-a).b.ux(i)^2+x(i).uy(i)^2].w(i)^2} 
;  Let b=tan(c)
;  db/dc = 1+tan(c)^2 = 1+b^2
;    dChi2/dc = dChi2/db . db/dc 
;             = SUM {-2[y(i)-a-b.x(i)](1+b^2).v(i)}
;  Setting to zero yields
;    SUM {(y(i)-a)(1+b^2).v(i)}  = b.SUM {x(i)(1+b^2).v(i)}
;  Since the (1+b^2) term is independent of i, it cancels.
;  So I don't see any difference in doing the minimization
;  in angle versus gradient space.
;
; 3) Iteration scheme
;    Solve equations (8) and (9) simultaneously give an improved
;    estimate of a and b.  It is not perfect because tw, tv, twx, tvx
;    all depend on a and b.
;
; 4) Compute uncertainties in a and b
;   This is awkward to do properly due to the b in the denominator.
;   So instead, assume that the denominator is fixed at the value
;   [uy(i)**2+b^2*ux(i)^2]  and use the standard equations for the
;   uncertainties from solving eqtn (10) on the last iteration.
;   The resulting uncertainties are independent of the y-values,
;   but depend on the y-errors.
;   
;
;  Notes:
;   The point representing the weighted average of the x- and
;   y-values [xbar,ybar] does not, in general, lie on the fitted
;   line. If it did, the line fitting would be much simpler. 
;   To convince yourself of this consider two points:
;     [xa,ya] has a small x-error and a large y-error
;     [xb,yb] has a small y-error and a large x-error
;   The point [xbar,ybar] won't fall on the line connecting the a & b.


;      implicit none
;      integer*4 nobs,i,iter,mit,iflag
;      parameter (mit=30)
;      real*4 x(nobs),ux(nobs),y(nobs),uy(nobs)
;      real*4 a,ea,b,eb,awas,bwas,
;     & xdif,xbar,ydif,ybar,
;     & xmin,xmax,ymin,ymax,
;     & res,chi2on,pcc 
;      real*8 tiny,tw,tdd,tx,ty,txx,tyy,txy,
;     & twx,twvx,twy,twvy,tu,tvx,tvy,
;     & ux2,uy2,denom,dd,u,w,tot

;  Transform variables so that they all lie in the range -1 to +1
;  Do this about a weighted average [xbar, ybar] to minimize
;  numerical errors associated with computing the correlation
;  coefficient later.

      mit=30
      xmin=+1e+36
      ymin=+1e+36
      xmax=-1e+36
      ymax=-1e+36
      twx=0.0d0
      twvx=0.0d0
      twy=0.0d0
      twvy=0.0d0
      tiny=1.d-36
      for i=0,nobs-1 do begin
         if(x(i) gt xmax) then xmax=x(i)
         if(x(i) lt xmin) then xmin=x(i)
         if(y(i) gt ymax) then ymax=y(i)
         if(y(i) lt ymin) then ymin=y(i)
         ux2=ux(i)*1.0d0*ux(i)
         if(ux2 lt tiny) then ux2=tiny
         twx=twx+1/ux2
         twvx=twvx+x(i)/ux2
         uy2=uy(i)*1.0d0*uy(i)
         if(uy2 lt tiny) then uy2=tiny
         twy=twy+1/uy2
         twvy=twvy+y(i)/uy2
      endfor
      xbar=twvx/twx
      ybar=twvy/twy
;      print,' xbar, ybar = ',xbar,ybar
;      print,' xmin, xmax = ', xmin,xmax
;      print,' ymin, ymax = ', ymin,ymax
      xdif=max([xmax-xbar,xbar-xmin])+tiny
      ydif=max([ymax-ybar,ybar-ymin])+tiny

      xdif=1.
      ydif=1.
      xbar=0.
      ybar=0.
      for i=0,nobs-1 do begin
         x(i)=(x(i)-xbar)/xdif
         ux(i)=ux(i)/xdif
         y(i)=(y(i)-ybar)/ydif
         uy(i)=uy(i)/ydif
      endfor

;  Iteration loop
      a = 0.0
      b = 0.0  ; first iteration ignores the x-uncertainties
      iflag=long(0)
      iter=long(0)
      while (iter lt mit and iflag lt 4 ) do begin
;        print, iter,iflag,ybar+ydif*(a-b*xbar/xdif),b*ydif/xdif
         tw  = tiny
         tv  = tiny
         twx = tiny
         tvx = tiny
         twy = 0.0d0
         tvy = 0.0d0
         for i=long(0),nobs-1  do begin
            uy2 = uy(i)*1.0d0*uy(i)
            ux2 = ux(i)*1.0d0*ux(i)
            dd=uy2+b*ux2*b
            if(dd le tiny) then begin
              w = 1/tiny
              v = w*x(i)
            endif else begin
              w = 1./dd
              v = w*(x(i)*uy2+b*(y(i)-a)*ux2)*w
            endelse
            tw  = tw + w
            tv  = tv + v
            twx = twx + w*x(i)
            tvx = tvx + v*x(i)
            twy = twy + w*y(i)
            tvy = tvy + v*y(i)
         endfor
         awas=a
         bwas=b
         denom=tw*tvx-tv*twx
         if (denom le 0.) then begin
            print,' SLFBE: denom =< 0 ',tw,tvx,tv,twx
            b=tvy/tvx
            a=twy/tw
         endif else begin
            b=(tw*tvy-tv*twy)/denom
            a=(tvx*twy-twx*tvy)/denom
         endelse
         if(abs(a-awas) lt 0.0000001*(1.+abs(a))) then iflag=iflag+1
         if(abs(b-bwas) lt 0.0000001*(1.+abs(b))) then iflag=iflag+1
;         print,'iter,iflag,awas,a,bwas,b = ',iter,iflag,awas,a,bwas,b
         iter=iter+1
      endwhile
      if(iflag lt 4) then print,' SLFBE Failed to converge'
;      print,'slfbe: iter,nobs,tw,tv,twx,tvx,twy,tvy,denom = ',iter,nobs,tw,tv,twx,tvx,twy,tvy,denom

      ea=sqrt(tvx/denom)
      eb=sqrt(tw/denom)
;  Comment:  The uncertainties are independent of the y-values since the
;  terms used (tvx,tw,denom) don't involve twy or tvy. Essentially, I'm
;  assuming that there are no x-uncertainties, as in a normal SLF, for
;  the computation of ea and eb.  But for a and b I'm doing it properly.
 
;  Compute Pearson Correlation Coefficient (PCC).
;  Since weighting here is different, xbar and ybar
;  cannot be assumed zero, but can be assumed small.
      tdd=0.0d0
      tx= 0.0d0
      ty= 0.0d0
      txx=0.0d0
      tyy=0.0d0
      txy=0.0d0
      for i=0,nobs-1 do begin
         uy2 = uy(i)*1.0d0*uy(i)
         ux2 = ux(i)*1.0d0*ux(i)
         dd=uy2+b*ux2*b
         if(dd lt tiny) then dd=tiny
         tdd=tdd+1/dd
         tx =tx +x(i)/dd
         ty =ty +y(i)/dd
         txx=txx+(x(i)/dd)*x(i)
         tyy=tyy+(y(i)/dd)*y(i)
         txy=txy+(x(i)/dd)*y(i)
      endfor
      pcc = (tdd*txy-tx*ty)/sqrt((tdd*txx-tx^2)*(tdd*tyy-ty^2))

;  Compute Chi2/N (in transformed co-ordinates)
      tot=0.0d0
      for i=0,nobs-1 do begin
         dd=uy(i)^2+(b*ux(i))^2
         if(dd le tiny) then dd=tiny
         res = y(i)-a-b*x(i)
         tot = tot + (res/dd)*res
      endfor
      chi2on=sqrt(tot/(nobs-2))
;      print,'chi2on1 = ',chi2on

;  Convert A, EA, B, EB back to original co-ordinates.
      eb=eb*ydif/xdif
      ea=sqrt((ydif*ea)^2+(xbar*eb)^2)
      a=ybar+ydif*(a-b*xbar/xdif)
      b=b*ydif/xdif

;  Optional: Convert VX, UX and VY, UY back to original co-ordinates
      for i=0,nobs-1 do begin
         x(i)=xbar+xdif*x(i)
         ux(i)=ux(i)*xdif
         y(i)=ybar+ydif*y(i)
         uy(i)=uy(i)*ydif
      endfor
      return
      end
