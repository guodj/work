

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;              Make ascii file
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO make_ascii, single, t, basetime,     $
		numofvar, titlenam, unitnam, placement,       $
		numofpts, datafile, fileout

  openw,1,fileout

  for i=0,n_elements(datafile)-1 do		$
    if strlen(datafile(i)) gt 0 then 		$
      printf,1, 'File number ',tostr(i+1),' : ',datafile(i)

  printf,1,'Variable 1 : Year'
  printf,1,'Variable 2 : Month'
  printf,1,'Variable 3 : Day'
  printf,1,'Variable 4 : Hour'
  printf,1,'Variable 5 : Minute'
  printf,1,'Variable 6 : Second'

  for curvar = 0,numofvar do begin

    ytlr =  ' (' + strcompress(unitnam(placement(curvar,0))) + ')'

    ytl  = strcompress(titlenam(placement(curvar,0)))

    printf,1,'Variable ',tostr(curvar+7),' : ',ytl,' ',ytlr

  endfor

  printf,1,' '
  printf,1,'Format for data is : 6I3,',tostr(numofvar+1),'E11.3'

  format = '(6I3,'+tostr(numofvar+1)+'E11.3)'

  printf,1,' '

  done = 0

  currow = intarr(numofvar+1)
  dumdata = fltarr(numofvar+1)
  missing = -1.0e32
  ti = double(t)
  for i=0,numofvar do ti(i,*) = ti(i,*) + basetime(i)
  maxt = 2.0*max(ti)
  n = n_elements(ti(*,0))-1

  while not done do begin

    ctime = min(ti(placement(0:n,0),currow))

    for i=0,numofvar do					$
      if ti(placement(i,0),currow(i)) eq ctime then begin
	dumdata(i) = single(placement(i,0),currow(i))
	ti(placement(i,0),currow(i)) = maxt
	if currow(i) lt numofpts(placement(i,0)) then 	$
	  currow(i) = currow(i) + 1
      endif else dumdata(i) = missing

    c_r_to_a, itime, ctime

    printf,1,format=format,itime,dumdata

    loc = where(currow ne numofpts(placement(*,0)),count)
    if (min(ti) eq maxt) or (count eq 0) then done = 1

  endwhile

  close,1

  RETURN

END

