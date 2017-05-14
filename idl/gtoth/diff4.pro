;===========================================================================
function diff4,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using 4th order centered differencing
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff4 is intended for 2D arrays only'
   retall
endif

n1=siz(1)
n2=siz(2)

dadx=a

if direction eq 1 then begin
   if n1 lt 5 then begin
      print,'Cannot take 4th order X gradient of grid with less than 5 columns'
      retall
   endif
   dadx(2:n1-3,*)=(a(4:n1-1,*)-8*a(3:n1-2,*)+8*a(1:n1-4,*)-a(0:n1-5,*)) $
                 /(x(3:n1-2,*)-x(1:n1-4,*))/6;
   dadx(0,*)   =dadx(2,*)
   dadx(1,*)   =dadx(2,*)
   dadx(n1-2,*)=dadx(n1-3,*)
   dadx(n1-1,*)=dadx(n1-3,*)
endif
if direction eq 2 then begin
   if n2 lt 5 then begin
      print,'Cannot take 4th order Y gradient of grid with less than 5 rows'
      retall
   endif
   dadx(*,2:n2-3)=(a(*,4:n2-1)-8*a(*,3:n2-2)+8*a(*,1:n2-4)-a(*,0:n2-5)) $
                 /(x(*,3:n2-2)-x(*,1:n2-4))/6;
   dadx(*,0)   =dadx(*,2)
   dadx(*,1)   =dadx(*,2)
   dadx(*,n2-2)=dadx(*,n2-3)
   dadx(*,n2-1)=dadx(*,n2-3)
endif

return,dadx

end

