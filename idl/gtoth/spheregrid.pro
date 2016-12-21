;===========================================================================
pro spheregrid,nvector,vectors,x,w,xreg,wreg
;
;    Transform vector variables from x,y,z to radial,phi,z components
;
;===========================================================================
  on_error,2

  xreg=x
  xreg(*,*,*,0)=sqrt(x(*,*,*,0)^2+x(*,*,*,1)^2+x(*,*,*,2)^2)
  xreg(*,*,*,2)=-atan(x(*,*,*,2),x(*,*,*,0))
  xreg(*,*,*,1)=atan(x(*,*,*,1),sqrt(x(*,*,*,0)^2+x(*,*,*,2)^2))
  phi=xreg(*,*,*,2)
  theta=xreg(*,*,*,1)
  wreg=w
  sinphi=sin(phi)
  cosphi=cos(phi)
  sintheta=sin(theta)
  costheta=cos(theta)
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     ivz=ivy+1
     wreg(*,*,*,ivx)=(w(*,*,*,ivx)*cosphi-w(*,*,*,ivz)*sinphi)*costheta $
                     +w(*,*,*,ivy)*sintheta
     wreg(*,*,*,ivz)=-w(*,*,*,ivx)*sinphi-w(*,*,*,ivz)*cosphi
     wreg(*,*,*,ivy)=(-w(*,*,*,ivx)*cosphi+w(*,*,*,ivz)*sinphi)*sintheta $
                     +w(*,*,*,ivy)*costheta
  endfor

  ;Remove 2*pi jumps from phi
  pi=4*atan(1) & pi2=2*pi & sz=size(phi) & nx2=sz(2) & nx3=sz(3)
  for ix3=1,nx3-1 do while phi(1,1,ix3-1) gt phi(1,1,ix3) do $
     phi(*,*,ix3)=phi(*,*,ix3)+pi2

  ;Remove turn over from theta
  for ix2=1,nx2-1 do $
  if theta(1,ix2-1,1) ge theta(1,ix2,1) then begin
     if theta(1,ix2,1) lt 0 then $
          theta(*,ix2-1,*)=-pi-theta(*,ix2-1,*) $
     else $
          theta(*,ix2,*)=pi-theta(*,ix2,*)
  endif

  xreg(*,*,*,2)=phi
  xreg(*,*,*,1)=theta
end

