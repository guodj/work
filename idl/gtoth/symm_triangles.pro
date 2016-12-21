;==============================================================================
pro symm_triangles,xx,yy,triangles,$
         triangconj,ntriangles,rectangles,nrectan

   ntriangles=n_elements(triangles(0,*))
   print,'Triangulation includes ',ntriangles, '  triangles'
   print,'Checking triangulation ...'

   npoints=n_elements(xx)

   dist=dblarr(npoints-1)
   i = 0L
   for i=0L,npoints-2 do $
      dist(i)=(xx(i+1)-xx(i))^2+(yy(i+1)-yy(i))^2
   dist2=min(dist)
   rectangles=lonarr(3,ntriangles)
   ;Structure of the rectangles array:
   ;If(rectangles(0,i)=1 then the Ith triangle from the triangles array
   ;is rectangular one
   tmp1=lonarr(3) & nrec_tri=0
   for i=0L,ntriangles-1 do begin
      if abs((xx(triangles(0,i))-xx(triangles(1,i)))*$
         (xx(triangles(1,i))-xx(triangles(2,i)))+$
         (yy(triangles(0,i))-yy(triangles(1,i)))*$
         (yy(triangles(1,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
      then begin
	 rectangles(0,i)=1
         tmp1(0)=triangles(1,i)
         if xx(triangles(0,i)) lt xx(triangles(2,i)) then begin
            tmp1(1)=triangles(0,i)
            tmp1(2)=triangles(2,i)
         endif else begin
            tmp1(1)=triangles(2,i)
            tmp1(2)=triangles(0,i)
         endelse
         for j=0,2 do triangles(j,i)=tmp1(j)
      endif

      if abs((xx(triangles(0,i))-xx(triangles(1,i)))*$
         (xx(triangles(0,i))-xx(triangles(2,i)))+$
         (yy(triangles(0,i))-yy(triangles(1,i)))*$
         (yy(triangles(0,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
      then begin
         rectangles(0,i)=1
         tmp1(0)=triangles(0,i)
         if xx(triangles(1,i)) lt xx(triangles(2,i)) then begin
            tmp1(1)=triangles(1,i)
            tmp1(2)=triangles(2,i)
         endif else begin
            tmp1(1)=triangles(2,i)
            tmp1(2)=triangles(1,i)
         endelse
         for j=0,2 do triangles(j,i)=tmp1(j)
      endif

      if abs((xx(triangles(0,i))-xx(triangles(2,i)))*$
         (xx(triangles(1,i))-xx(triangles(2,i)))+$
         (yy(triangles(0,i))-yy(triangles(2,i)))*$
         (yy(triangles(1,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
      then begin
         rectangles(0,i)=1
         tmp1(0)=triangles(2,i)
         if xx(triangles(0,i)) lt xx(triangles(1,i)) then begin
            tmp1(1)=triangles(0,i)
            tmp1(2)=triangles(1,i)
         endif else begin
            tmp1(1)=triangles(1,i)
            tmp1(2)=triangles(0,i)
         endelse
         for j=0,2 do triangles(j,i)=tmp1(j)
      endif
   endfor
   ;Rectangles(1,i) is not equal to zero if the ith rectangular triandgle
   ;has a common long side with the jth rectangular triangle. In this case
   ;rectangles(2,i)=j
   nrectan=0L
   for i=0L,ntriangles-1 do begin
     if rectangles(0,i) gt 0 then begin
        nrec_tri=nrec_tri+1
        if rectangles(1,i) eq 0 then begin
        for j=i+1L,ntriangles-1 do begin
           if rectangles(0,j) gt 0 then $
           if triangles(1,i) eq triangles(1,j) then $
           if triangles(2,i) eq triangles(2,j) then begin
              nrectan=nrectan+1
              rectangles(1,i)=1
              rectangles(2,i)=j
              rectangles(1,j)=1
              rectangles(2,j)=i
              goto,out
           endif
        endfor
        out:
        endif
     endif
   endfor

   if nrectan ne 0  then begin
      print,'Among    ',nrec_tri, '  rectangular triangles'
      print,'there are',nrectan, '   pairs which have common circumcircle'
      tmp2=lonarr(4)
      ndiag1=0
      ndiag2=0
      triangconj=lonarr(3,ntriangles)
      for i=0L,ntriangles-1 do begin
         if rectangles(1,i) eq 1 then begin
            if rectangles(2,i) gt i then begin
               for j=0,2 do tmp2(j)=triangles(j,i)
               tmp2(3)=triangles(0,rectangles(2,i))
               if yy(tmp2(1)) lt yy(tmp2(2)) then ndiag1=ndiag1+1 else $
                  ndiag2=ndiag2+1
               triangconj(0,i)=tmp2(1)
               triangconj(1,i)=tmp2(0)
               triangconj(2,i)=tmp2(3)
               triangconj(0,rectangles(2,i))=tmp2(2)
               triangconj(1,rectangles(2,i))=tmp2(0)
               triangconj(2,rectangles(2,i))=tmp2(3)
            endif
         endif else for j=0,2 do triangconj(j,i)=triangles(j,i)
      endfor
      print,' Among them ',ndiag1, ' are formed by the triangles,'
      print,' having the common side which is oriented as /////'
      print,' and ',ndiag2, ' have the triangles common side,'
      print,' oriented as \\\\\'
      print,' Calculating the conjugated triangulation ...'
   endif

end
