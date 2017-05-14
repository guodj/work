;===========================================================================
function diff3,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using IDL's 1D deriv() function
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff3 is intended for 2D arrays only'
   retall
endif

dadx=a

if direction eq 1 then for i2=0,siz(2)-1 do dadx(*,i2)=deriv(x(*,i2),a(*,i2))
if direction eq 2 then for i1=0,siz(1)-1 do dadx(i1,*)=deriv(x(i1,*),a(i1,*))

return,dadx

end

