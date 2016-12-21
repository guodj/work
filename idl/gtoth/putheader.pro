;===========================================================================
pro putheader,multix,multiy,ix,iy,ninfo,headline,nx

on_error,2

if ninfo lt 1 then return
info=strtrim(headline,2)
if ninfo gt 1 then info=info+' (nx='+string(nx,format='(i6,2(i4))')+')'
xyouts,5+(ix*!d.x_size)/multix,-12+((iy+1)*!d.y_size)/multiy,/DEV,info

end
