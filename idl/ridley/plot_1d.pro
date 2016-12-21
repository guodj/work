
filelist = findfile("b0*_t*.*ALL")

nfiles = n_elements(filelist)

;iVar = 7
iVar = 12
iVar = 3
iVar = 4

for ifile = 0,nfiles-1 do begin

  file  = filelist(ifile)

  read_thermosphere_file, file, nvars, nalts, nlats, nlons, $
    vars, data, nBLKlat, nBLKlon, nBLK

  if (nBlk eq 1 and nlats eq 1 and nlons eq 1) then begin

      data = reform(data)

      altitude = data(2,*) / 1000.0
      latitude  = data(1,0)
      longitude = data(0,0)

      print, "rho : ",mm(alog10(data(3,*)))
      print, "temp : ",mm(data(4,*))
      print, "vel : ",mm(data(11,*))

;      plot, data(7,*), altitude

if (iFile eq 0) then begin
      plot, data(iVar,2:nAlts-3), altitude, xrange=[0,500],$
        yrange = [0,max(altitude)], title = file, ystyle = 1
  endif else oplot, data(iVar,2:nAlts-3), altitude

;if (iFile eq 0) then begin
;      plot, data(iVar,2:nAlts-3), altitude, xrange=[200,2000],$
;        yrange = [0,max(altitude)], title = file, ystyle = 1
;  endif else oplot, data(iVar,2:nAlts-3), altitude

;      oplot, data(4,*), altitude, linestyle = 1

;      oplot, data(10,*)*20, altitude, linestyle = 2
;      oplot, data(11,*), altitude, linestyle = 1

      wait,5.0/nfiles

  endif

endfor

stop

;plot_oi, data(iVar,*), altitude, xtitle = "NDensity", ytitle = 'Altitude', $
;  xrange = mm(data(3:nVars-3,*))
;
;iAlt = 50
;
;for iVar = 3, nVars-3 do begin
;
;    iAlt = nalts * randomu(s,1)
;
;    oplot, data(iVar,*), altitude, linestyle = iVar-3
;    oplot, [data(iVar,iAlt)], [altitude(iAlt)], psym = 4
;    xyouts, data(iVar,iAlt), altitude(iAlt), vars(iVar)
;
;endfor

plot_oi, data(7,*), altitude, xrange = [1.0e10,1.0e13]
oplot, data(12,*), altitude, linestyle = 1

oplot, data(6,*), altitude
oplot, data(11,*), altitude, linestyle = 2

;plot_oi, data(13,*), altitude
;oplot, data(14,*), altitude, linestyle = 1

end
