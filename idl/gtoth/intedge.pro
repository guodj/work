;===========================================================================
function intedge,f,xc
;
; Integrate the neighbouring values of "f" for the four edges described by "xc"
; The size of "f", "xc", and the result are n1*n2, (n1-1)*(n2-1), and n1*n2
; respectively, but only the inner (n1-2)*(n2-2) points are calculated, the
; edge values are 0-s.
;===========================================================================

siz=size(f)
n1=siz(1)
n2=siz(2)

intf=dblarr(n1,n2)
intf(1:n1-2,1:n2-2)=-(xc(1:n1-2,1:n2-2)-xc(0:n1-3,1:n2-2))*f(1:n1-2,2:n2-1) $
                    -(xc(1:n1-2,0:n2-3)-xc(1:n1-2,1:n2-2))*f(2:n1-1,1:n2-2) $
                    -(xc(0:n1-3,0:n2-3)-xc(1:n1-2,0:n2-3))*f(1:n1-2,0:n2-3) $
                    -(xc(0:n1-3,1:n2-2)-xc(0:n1-3,0:n2-3))*f(0:n1-3,1:n2-2)

return,intf

end

