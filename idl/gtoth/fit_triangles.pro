;==============================================================================

pro fit_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
        triangles,ntriangles,rectangles

   tmp2=lonarr(4)

   for iw=0,nw-1 do begin
      data=reform(w(*,*,iw))
      print,'Calculating the fitting triangulization for iw=',iw
      for i=0L,ntriangles-1 do begin
         if rectangles(1,i) eq 1 then begin
            if rectangles(2,i) gt i then begin
               tmp2(0:2)=triangles(0:2,i)
               tmp2(3)  =triangles(0,rectangles(2,i))
               if abs(data(tmp2(0))-data(tmp2(3))) lt $
                  abs(data(tmp2(1))-data(tmp2(2))) then begin
                  triangles(0,i)=tmp2(1)
                  triangles(1,i)=tmp2(0)
                  triangles(2,i)=tmp2(3)
                  triangles(0,rectangles(2,i))=tmp2(2)
                  triangles(1,rectangles(2,i))=tmp2(0)
                  triangles(2,rectangles(2,i))=tmp2(3)
               endif
            endif
         endif
      endfor
      wreg(*,*,iw)=trigrid(xx,yy,data,triangles, $
         [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))
   endfor
   print,'Using fitted triangulation'
end

