;===================================================================
pro plotgrid,x,y,lines=lines,xstyle=xstyle,ystyle=ystyle,polar=polar
;===================================================================

on_error,2

if not keyword_set(x) then begin
   print,'Usage: plotgrid, x [,y] [,/lines], [/polar], [,xstyle=3] [,ystyle=1]
   retall
endif

xx=reform2(x)
sizx=size(xx)

if (n_elements(polar) eq 0) then polar = 0

if not keyword_set(y) then begin
   case sizx(0) of
     3:begin
         if sizx(3) ne 2 then goto, ERROR1
         yy=xx(*,*,1)
         xx=xx(*,*,0)
       end
     2:begin
         if sizx(2) ne 2 then goto, ERROR1
         yy=xx(*,1)
         xx=xx(*,0)
         lines=0
       end
   else: goto, ERROR1
   endcase
endif else begin
   yy=reform2(y)
   sizy=size(yy)
   if sizx(0) ne sizy(0)            then goto, ERROR2
   if max(abs(sizx-sizy)) ne 0      then goto, ERROR2
   if sizx(0) ne 2 and sizx(0) ne 1 then goto, ERROR2
   if sizx(0) eq 1 then lines=0
endelse

if keyword_set(lines) then begin
   plot, xx, yy, XSTYLE=xstyle, YSTYLE=ystyle, POLAR=polar, /NOERASE, /NODATA

   for ix=0,sizx(1)-1 do begin
      oplot,xx(ix,*),yy(ix,*),POLAR=polar
   endfor
   for iy=0,sizx(2)-1 do begin
      oplot,xx(*,iy),yy(*,iy),POLAR=polar
   endfor
endif else begin
   if polar then $
      plot, xx, yy, PSYM=3, SYMSIZE=!p.symsize, $
         XSTYLE=xstyle, YSTYLE=ystyle, /NOERASE, /POLAR  $
   else $   
      plot, xx, yy, PSYM=1, SYMSIZE=!p.symsize, $
         XSTYLE=xstyle, YSTYLE=ystyle, /NOERASE
endelse

return

ERROR1:
   print,'size(x)=',sizx
   print,'Error: plotgrid,x  requires x(nx,ny,2) array'
   retall

ERROR2:
   print,'size(x)=',sizx,' size(y)=',sizy
   print,'Error: plotgrid,x,y requires x(nx,ny) y(nx,ny) arrays'
   retall


end

