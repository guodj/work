
pro read_madrigal, file, time, data, vars

  spawn, 'wc '+file, wc
  nLines = fix(wc(0))-1

  close,1
  openr,1,file

  line = ''

  readf,1,line
  vars = strsplit(line,' ',/extract)

  nVars = n_elements(vars)-6
  Vars = Vars(6:nVars+6-1)

  data = fltarr(nVars, nLines)
  time = dblarr(nLines)

  for iLine = 0, nLines-1 do begin

     readf,1,line
     d = strsplit(line,' ',/extract)
     itime = fix(d(0:5))
     c_a_to_r, itime, t
     time(iLine) = t

     for i = 0, nVars-1 do begin
        if (strpos(d(6+i),'missing') eq 0) then d(6+i)='-1.0e32'
        if (strpos(d(6+i),'NaN') eq 0) then d(6+i)='-1.0e32'
        if (strpos(d(6+i),'assumed') eq 0) then d(6+i)='1.0'
        data(i,iLine) = float(d(6+i))
     endfor

  endfor

  close,1

end
