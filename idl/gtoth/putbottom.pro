;===========================================================================
pro putbottom,multix,multiy,ix,iy,ninfo,nx,it,time

on_error,2

if ninfo lt 1 then return
info=''
if ninfo gt 2 then info='nx='+string(nx,format='(i6,2(",",i4))')+' '
if ninfo gt 1 then info=info+'it='+string(it,format='(i6)')+', '
info=info+'time='+string(time,format='(g12.5)')
xyouts,5+(ix*!d.x_size)/multix,8+(iy*!d.y_size)/multiy,/DEV,info

end
