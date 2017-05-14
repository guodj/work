;===========================================================================
function triplet,x0,x1,dx,y0,y1,dy,z0,z1,dz,w0,w1,dw
;
; Produce an index array corresponding to the Fortran 90 triplet notation
;
; Usage: cut=triplet(0,99,2,0,99,2)
;
;        velvector=25*25  &  velpos=dblarr(velvector,2)
;        velpos(*,*)=x(triplet(0,99,4,0,99,4,0,1,1))
;
; Note: the resulting indices are valid for an array of size
;
;      (x1+1)*(y1+1)*(z1+1)*(w1+1)
;
;===========================================================================

if keyword_set(dw) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy,z1+1,z0,z1,dz,w1+1,w0,w1,dw)

if keyword_set(dz) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy,z1+1,z0,z1,dz)

if keyword_set(dy) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy)

if keyword_set(dx) then $
   return,quadruplet(x1+1,x0,x1,dx)

print,'Error in TRIPLET: All strides are 0!'
retall

end

