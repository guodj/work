;===========================================================================
function diff2,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using 2nd order centered differencing
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff2 is intended for 2D arrays only'
   retall
endif

n1=siz(1)
n2=siz(2)

if direction eq 1 then begin
   ind1=indgen(n1)
   jnd1=ind1+1
   jnd1(n1-1)=n1
   hnd1=ind1-1
   hnd1(0)=0
   dadx=(a(jnd1,*)-a(hnd1,*))/(x(jnd1,*)-x(hnd1,*))
endif
if direction eq 2 then begin
   ind2=indgen(n2)
   jnd2=ind2+1
   jnd2(n2-1)=n2
   hnd2=ind2-1
   hnd2(0)=0
   dadx=(a(*,jnd2)-a(*,hnd2))/(x(*,jnd2)-x(*,hnd2))
endif

return,dadx

end

