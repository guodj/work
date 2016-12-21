;===========================================================================
function symmdiff,direction,a
;
; Take symmetric difference of "a" with respect to a mirror plane in direction
; "direction"
;
;===========================================================================
on_error,2

siz=size(a)
dim=siz(0)
nx=siz(1)

diff=a

case dim of
1: for i=0,nx-1 do diff(i)=a(i)-a(nx-1-i)
2: begin
     ny=siz(2)
     case direction of
     1: for i=0,nx-1 do diff(i,*)=a(i,*)-a(nx-1-i,*)
     2: for i=0,ny-1 do diff(*,i)=a(*,i)-a(*,ny-1-i)
     endcase
   end
3: begin
     ny=siz(2)
     nz=siz(3)
     case direction of
     1: for i=0,nx-1 do diff(i,*,*)=a(i,*,*)-a(nx-1-i,*,*)
     2: for i=0,ny-1 do diff(*,i,*)=a(*,i,*)-a(*,ny-1-i,*)
     3: for i=0,nz-1 do diff(*,*,i)=a(*,*,i)-a(*,*,nz-1-i)
     endcase
   end
4: begin
     ny=siz(2)
     nz=siz(3)
     case direction of
     1: for i=0,nx-1 do diff(i,*,*,*)=a(i,*,*,*)-a(nx-1-i,*,*,*)
     2: for i=0,ny-1 do diff(*,i,*,*)=a(*,i,*,*)-a(*,ny-1-i,*,*)
     3: for i=0,nz-1 do diff(*,*,i,*)=a(*,*,i,*)-a(*,*,nz-1-i,*)
     endcase
   end
endcase

return,diff
end

