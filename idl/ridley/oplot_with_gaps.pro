
pro oplot_with_gaps, x, y, gap, linestyle = linestyle, thick = thick, $
                     color = color

  if n_elements(linestyle) eq 0 then linestyle = 0
  if n_elements(thick) eq 0 then thick = 4
  if n_elements(color) eq 0 then color = 0

  n = n_elements(x)
  dx = x(1:n-1) - x(0:n-2)

  l = where(dx gt gap, c)

  if c eq 0 then begin
     l = n
     c = 1
  endif

  iS = 0
  iE = l(0)
  if (iE eq n) then iE = n-1
  oplot, x(iS:iE), y(iS:iE), linestyle = linestyle, thick = thick, $
         color = color

  for i = 1, c-1 do begin
     iS = l(i-1)+1
     iE = l(i)
     oplot, x(iS:iE), y(iS:iE), linestyle = linestyle, thick = thick, $
            color = color
     if (iS eq iE) then $
        oplot, x(iS:iE), y(iS:iE), linestyle = linestyle, thick = thick, $
               color = color, psym = 4

  endfor

end
