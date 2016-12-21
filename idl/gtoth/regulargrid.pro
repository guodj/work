
;===========================================================================
pro regulargrid,x_old,nxreg_old,xreglimits_old,x,xreg,nxreg,xreglimits,$
                w,wreg,nw,wregpad,triangles,symmtri
;
;    Regularize grid and interpolate "w" via triangulate and trigrid.
;    The original "w" data is interpolated into "wreg", for points outside
;    the convex hull of the irregular grid the "wregpad(nw)" array is used.
;
;    If "x_old" and "x" or "nxreg_old" and "nxreg" are different
;    a triangulization is done first and a regular coordinate array
;    "xreg" is created. The size of the "xreg" array is in "nxreg(2)",
;    "xreglimits(4)" gives the limits. The triangles are saved in "triangles".
;
;    "q" can be interpolated from the irregular grid to the regular one by:
;
;    qreg(*,*)=trigrid(x(*,*,0),x(*,*,1),q,triangles,[0,0],xreglimits)
;
;===========================================================================
   on_error,2

   ;Floating underflow is not a real error, the message is suppressed
   err=check_math(1,1)

   xx=x(*,*,0)
   yy=x(*,*,1)

   ; Test distribution. If you discomment the next lines you can
   ; take a look at the different "shock wave" representation
   ; on your grid for the 0-th variable (usually rho)
   ; for i=0L,n_elements(xx)-1 do $
   ;   if abs(xx(i))-0.2*abs(yy(i)) gt 0.01 then $
   ;       w(i,*,0)=2. else w(i,*,0)=1.

   ; Check if nxreg==nxreg_old and xreglimits==xreglimits_old and x==x_old
   newx=1
   nrectan=0
   if symmtri ne 1 and symmtri ne 2 then $
   if n_elements(nxreg_old) eq n_elements(nxreg) then $
   if max(abs(nxreg_old-nxreg)) eq 0 then $
   if n_elements(xreglimits) eq n_elements(xreglimits_old) then $
   if max(abs(xreglimits-xreglimits_old)) eq 0 then $
   if n_elements(x_old) eq n_elements(x) then $
   if max(abs(x_old-x)) eq 0 then newx=0

   if xreglimits(0) eq xreglimits(2) then begin
      xreglimits(0)=min(xx) & xreglimits(2)=max(xx)
   endif
   if xreglimits(1) eq xreglimits(3) then begin
      xreglimits(1)=min(yy) & xreglimits(3)=max(yy)
   endif

   if newx then begin
      print,'Triangulating...'
      x_old=x
      nxreg_old=nxreg
      xreglimits_old=xreglimits

      triangulate,xx,yy,triangles

      ; calculate conjugate triangulation and rectangles if required
      if symmtri eq 1 or symmtri eq 2 then $
          symm_triangles,xx,yy,triangles,$
                         triangconj,ntriangles,rectangles,nrectan

      xreg=dblarr(nxreg(0),nxreg(1),2)
      dx=(xreglimits(2)-xreglimits(0))/(nxreg(0)-1)
      dy=(xreglimits(3)-xreglimits(1))/(nxreg(1)-1)
      for i=0,nxreg(1)-1 do xreg(*,i,0)=dx*indgen(nxreg(0))+xreglimits(0)
      for i=0,nxreg(0)-1 do xreg(i,*,1)=dy*indgen(nxreg(1))+xreglimits(1)

      wreg=dblarr(nxreg(0),nxreg(1),nw)
   endif
   if not keyword_set(wregpad) then begin
      wregpad=dblarr(nw)
      for iw=0,nw-1 do begin
         wmax=max(w(*,*,iw))
         wmin=min(w(*,*,iw))
         if wmax*wmin lt 0 then wregpad(iw)=0 $
         else                   wregpad(iw)=wmin-0.1*(wmax-wmin)
      endfor
   endif

   case 1 of

   symmtri eq 3: for iw=0,nw-1 do $
      wreg(*,*,iw)=grid_data(xx,yy,reform(w(*,*,iw)),nxreg,xreglimits,$
                             triangles,wregpad(iw))

   symmtri eq 0 or (symmtri lt 3 and nrectan eq 0): for iw=0,nw-1 do $
      wreg(*,*,iw)=trigrid(xx,yy,w(*,*,iw),triangles, $
           [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))

   symmtri eq 1 and nrectan gt 0: $
      fit_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
           triangles,ntriangles,rectangles

   symmtri eq 2 and nrectan gt 0: $
      average_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
           triangles,triangconj

   endcase

   err=check_math(0,0)
   ;Floating underflow is not a real error, the message is suppressed
   if err ne 32 and err ne 0 then print,'Math error in regulargrid:',err

end

