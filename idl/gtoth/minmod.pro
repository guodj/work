;===========================================================================
function minmod,a,b
;
; Calculate minmod limited slope of a and b slopes

on_error,2

; get sign of a
if a gt 0 then s=1 else s=-1

; calculate limited slope
c = s*max([0,min([abs(a),s*b])])

return,c
end
