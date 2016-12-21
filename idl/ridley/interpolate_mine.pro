
; outx is the time you want to interpolate TO
; inx is the time you want to interpolate FROM
; iny is the thing you want to interpolate

function interpolate_mine, outx, iny, inx

  nx = n_elements(outx)
  outy = fltarr(nx)

  for i=0L,nx-1 do begin

      if (outx(i) lt min(inx)) then begin
          outy(i) = iny(0)
      endif else begin
          if (outx(i) gt max(inx)) then begin
              outy(i) = iny(n_elements(iny)-1)
          endif else begin
              d = outx(i) - inx
              if (min(abs(d)) eq 0) then begin
                  l = where(d eq 0)
                  outy(i) = iny(l(0))
              endif else begin
                  lp = where(d gt 0,cp)
                  lp = lp(cp-1)
                  lm = where(d lt 0,cm)
                  lm = lm(0)
                  x = (outx(i) - inx(lm))/ (inx(lp) - inx(lm))
                  outy(i) = x*iny(lp) + (1-x)*iny(lm)
              endelse
          endelse
      endelse

  endfor

  return, outy

end
