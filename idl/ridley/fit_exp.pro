
pro fit_exp, datax, datay, y0, ymax, tau, rms, linear = linear, $
             absolute=absolute, average = average

  if (n_elements(ymax) eq 0) then ymax = -1.0e32
  if (n_elements(tau)  eq 0) then tau  = -1.0e32
  if (n_elements(y0)  eq 0)  then tau  = -1.0e32
  if (n_elements(linear) eq 0) then linear = 0 else linear = 1
  if (n_elements(absolute) eq 0) then absolute = 0 else absolute = 1
  if (n_elements(average) eq 0) then average = 0

  newx = datax
  newy = datay

  if average gt 0 then bin_data, datax, datay, nPts

  if (y0 lt -1.0e31) then y0 = min(datay)

  dy = max(newy)-min(newy)
  dx = max(newx)-min(newx)

  if (ymax lt -1.0e31) then ymax = dy*3

  if (tau lt -1.0e31) then begin
      if (linear) then tau = dy/dx else tau = dx
  endif

  if (linear) then fit = y0 + tau*newx $
  else fit = y0 + ymax*(1.0-exp(-newx/tau))
  if absolute then rms = mean(abs(fit-newy)) $
  else rms = sqrt(mean((fit-newy)^2))

  tau_save  = tau
  ymax_save = ymax
  y0_save   = y0
  rms_save  = rms

  print, "Init Value ", tau_save, y0_save, ymax_save, rms_save

  for i=1,10000 do begin

      v = randomu(s,3)-0.5

      tau   = tau_save   + 0.1*tau_save*(v(0))
      y0    = y0_save    + 0.1*y0_save*(v(1))
      ymax  = ymax_save  + 0.1*ymax_save *(v(2))

      if (linear) then fit = y0 + tau*newx $
      else fit = y0 + ymax*(1.0-exp(-newx/tau))

      if absolute then rms = mean(abs(fit-newy)) $
      else rms = sqrt(mean((fit-newy)^2))

      if (rms lt rms_save) then begin

          tau_save = tau
          ymax_save = ymax
          y0_save = y0
          rms_save = rms

          print, "New Value ", i, tau_save, y0_save, ymax_save, rms_save

;;          plot, newx, newy, psym = 3
;;          x = findgen(1001)
;;
;;          if linear then y = y0 + tau * x $
;;          else y = y0 + ymax*(1.0-exp(-x/tau))
;;          oplot, x, y, linestyle = 1
;;;          wait,1.0
          
      endif

  endfor

  tau = tau_save
  ymax  = ymax_save
  rms = rms_save

  if (linear) then fit = y0 + tau*datax $
  else fit = y0 + ymax*(1.0-exp(-datax/tau))

  if absolute then rms = mean(abs(fit-datay)) $
  else rms = sqrt(mean((fit-datay)^2))

end
