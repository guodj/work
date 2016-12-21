;==========================================
function reform2,x
;==========================================
  ;Remove all degenerate dimensions from x

  if n_elements(x) lt 2 then return,x

  siz=size(x)
  siz=siz(1:siz(0))
  return,reform(x,siz(where(siz gt 1)))

end

