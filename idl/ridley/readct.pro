;*****************************************************************************
 
pro readct, ncolors, ctname
 
;*****************************************************************************
 
  openr,11,ctname
  ncolors=0
  readf,11,ncolors
  color = fltarr(3,255)
  readf,11,color
  close,11

  tvlct,color(0,*),color(1,*),color(2,*)
 
  return
 
end
