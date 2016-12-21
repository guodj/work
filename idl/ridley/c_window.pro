pro c_window, npts, wind, window_shape

; window shape :
;   1 - no window (rectangular)
;   2 - trianglular
;   3 - blackman
;   4 - hamming
;   5 - hanning
;   6 - kaiser
;   7 - lanczos
;   8 - tukey
;   these shapes are given in Introduction to Digital Signal Processing
;   by Proakis and Manolakis

  M = npts + 1

  case (window_shape) of

    1 : wind = fltarr(M)+1.0
    2 : wind = 1.0 - abs(2.0*(findgen(M)-((M-1)/2.0))/(M-1))
    3 : wind = 0.42 - 0.5*cos(2.0*!pi*findgen(M)/(M-1)) + 		$
		     0.08*cos(4.0*!pi*findgen(M)/(M-1))
    4 : wind = 0.54 -0.46*cos(2.0*!pi*findgen(M)/(M-1))
    5 : wind = 0.5*(1.0-cos(2.0*!pi*findgen(M)/(M-1)))

  endcase

  return

end

