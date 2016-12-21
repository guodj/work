;==========================================
function coarse,a,boxsize
;
; Produce a coarser array from "a" by averaging out cells in a box.
; The box size can be defined by a scalar (n long interval, n*n squarle,
; or ,n*n*n cube) or as an
; array of the same dimension as "a" (n1*n2 rectangle or n1*n2*n3 brick)

;on_error,2

if(n_elements(a) eq 0 or n_elements(boxsize) eq 0)then begin
   print,'Calling sequence is: array_co=coarse(array, boxsize)'
   retall
endif

siz=size(a)
ndim=siz(0)

if(ndim eq 0 or ndim gt 3)then begin
   print,'coarse requires a 1,2 or 3D array for the 1st argument'
   retall
endif
nx=siz(1:ndim)

siz=size(box)
if(siz(0) eq 0)then begin
   n=intarr(ndim)+boxsize
endif else if siz(0) eq ndim then begin
   n=boxsize
endif else begin
   print,'boxsize should either be a scalar, or an array '
   print,'of the same dimension as the number of dimensions of the array'
   retall
endelse

case ndim of
   1: begin
      result=dblarr(nx(0)/n(0))
      for ix=0,(nx(0)-1)/n(0) do $
        for i=0,n(0)-1 do $
           result(ix)=result(ix)+a(ix*n(0)+i)
      result=result/n(0)
   end
   2: begin
      result=dblarr(nx(0)/n(0),nx(1)/n(1))
      for ix=0,(nx(0)-1)/n(0) do $
      for iy=0,(nx(1)-1)/n(1) do $
        for i=0,n(0)-1 do $
        for j=0,n(1)-1 do $
           result(ix,iy)=result(ix,iy)+a(ix*n(0)+i,iy*n(1)+j)
      result=result/n(0)/n(1)
   end
   3: begin
      result=dblarr(nx(0)/n(0),nx(1)/n(1),nx(2)/n2)
      for ix=0,(nx(0)-1)/n(0) do $
      for iy=0,(nx(1)-1)/n(1) do $
      for iz=0,(nx(2)-1)/n(2) do $
        for i=0,n(0)-1 do $
        for j=0,n(1)-1 do $
        for k=0,n(2)-1 do $
           result(ix,iy,iz)=result(ix,iy,iz)+a(ix*n(0)+i,iy*n(1)+j,iz*n(2)+k)
      result=result/n(0)/n(1)/n(2)
   end
endcase
return,result
end

