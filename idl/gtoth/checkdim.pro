;===========================================================================
pro checkdim,idim,nx,x0,x1,dx

; Check quadruplet for conditions nx>x1>=x0>=0 and dx>0
;===========================================================================

   if nx le 0 then begin
      print,'Size must be positive for dimension',idim
      retall
   endif
   if x1 ge nx then begin
      print,'Maximum index must be less than size for dimension',idim
      retall
   endif
   if x0 lt 0 then begin
      print,'Minimum index must be greater than 0 for dimension',idim
      retall
   endif
   if x0 gt x1 then begin
      print,'Minimum index must be less than maximum index for dimension',idim
      retall
   endif
   if dx le 0 then begin
      print,'Stride must be a positive integer for dimension',idim
      retall
   endif

return
end

