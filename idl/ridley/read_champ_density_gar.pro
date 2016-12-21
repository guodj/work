
pro read_champ_density_gar, file, time, data, vars, note

  f = findfile(file)
  if (strlen(f(0)) le 1) then begin

     time = dblarr(1)
     time(0) = -999.0

  endif else begin

     spawn, 'wc '+file,wc
     nLines = fix(wc)
     nLines = nLines(0)-2

     close,1
     openr,1,file
     line = ''
     note = ''
     readf,1,note
     readf,1,line
     vars = strsplit(line, ';', /extract)
     nVars = n_elements(vars)
     data = fltarr(nVars,nLines)
     readf,1,data
     close,1

     time = dblarr(nLines)
     for i=0,nLines-1 do begin
        itime = [data(0,i),1,data(1,i),0,0,data(2,i)]
        c_a_to_r, itime, rTime
        time(i) = rTime
     endfor
     p=plot, time,data(12,*),TITLE='Time vs Neutral density at 400 km'
     p.save, "neutral density.png"
     p=plot, time, data(06,*), TITLE='Time vs Satellite Height'
     p.save, "height.png"
     ;plot, data(10,*), data(12,*), TITLE='Satellite M Local Time vs Neutral Density '
     ;plot, data(08,*), data(12,*),TITLE='Satellite QD Latitude vs Neutral Density'
     ;plot, data(09,*), data(12,*), TITLE='Satellite Mlon vs Neutral density'
     ;plot, data(06,*), data(12,*), TITLE='Satellite Height vs Neutral Density'
  endelse

end
