;===========================================================================
pro getaxes,ndim,x,xx,yy,zz,cut,cut0,rSlice,plotdim,variables
;===========================================================================
on_error,2
case ndim of
  1: xx=x
  2: begin
        xx=x(*,*,0)
        yy=x(*,*,1)
     end
  3: begin
       xx=x(*,*,*,0)
       yy=x(*,*,*,1)
       zz=x(*,*,*,2)
     end
endcase

if keyword_set(cut0) then begin
   xx=xx(cut0)
   if ndim gt 1 then yy=yy(cut0)
   if ndim gt 2 then zz=zz(cut0)
endif

!x.title=variables(0)
if plotdim gt 1 then !y.title=variables(1)
if plotdim gt 2 then !z.title=variables(2)

if ndim eq 3 and plotdim eq 2 then begin
   siz=size(cut)
   case 1 of
     siz(0) eq 2: rSlice=zz(0)
     siz(1) eq 1: rSlice=xx(0)
     siz(2) eq 1: rSlice=yy(0)
   endcase
   print,'Normal coordinate of 2D slice:',rSlice
endif else        rSlice=0.0

; Cut with fixed X value?
siz=size(cut)
; in 2D
if siz(0) eq 2 and siz(1) eq 1 then begin
   xx=yy
   !x.title=variables(1)
endif
; in 3D
if siz(0) eq 3 then begin
   case 1 of
   plotdim eq 1: begin
      xx=zz
      !x.title=variables(2)
   end
   siz(1) eq 1: begin
         xx=yy
         yy=zz
         !x.title=variables(1)
         !y.title=variables(2)
   end
   siz(2) eq 1: begin
      yy=zz
      !y.title=variables(2)
   end
   else: print,'internal error in getaxes'
   endcase
endif

end

