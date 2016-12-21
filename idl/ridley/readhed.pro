
PRO readhed, numoffil, datafile, timeint, titles, units, 	$
	count, nrows, rowlen, type, data_struct

  k=0

  max = 0
  nrows = lonarr(numoffil)

  if n_elements(timeint) eq 0 then timeint = strarr(numoffil,2) 
  if n_elements(titles) eq 0 then titles = strarr(numoffil,500) 
  if n_elements(units) eq 0 then units = titles
  if n_elements(type) eq 0 then type = titles
  if n_elements(count) eq 0 then count = intarr(numoffil) 
  if n_elements(rowlen) eq 0 then rowlen = intarr(numoffil) 

  while k lt numoffil do begin

    openr, k+1, datafile(k)+'.hed'
    x$    = ''
    x     = strarr(200)

    i=1

    While (strpos(x(i-1),'NOTES') eq -1)  do begin

      readf,k+1,x$
      x(i) = x$
      i=i+1
      if (i gt n_elements(x)) then x = [x,strarr(200)]

    EndWhile

    nrows(k) = strmid(x(5),40,12)
    rowlen(k) = strmid(x(3),40,12)

    x=x(11:i-3)

    checker=0

    While (checker ne -1)  do begin

      readf,k+1,x$
      if (strpos(mklower(x$),'end') ne -1) then begin
        timeint(k,1)=x$
        checker=-1
      endif
                                                                        
      if (strpos(mklower(x$),'sta') ne -1) then timeint(k,0)=x$

    EndWhile

    close, k+1

    timeint(k,0)=strmid(timeint(k,0),14,25)
    timeint(k,1)=strmid(timeint(k,1),14,25)

    titles(k,0:i-14)=strmid(x,7,14)
    units(k,0:i-14)=strmid(x,21,14)	
    type(k,0:i-14)=strmid(x,57,13)
    count(k)=i-14

    if count(k) gt max then max = count(k)

    k=k+1

  EndWhile

  titles=titles(*,0:max)
  units=units(*,0:max)

  RETURN

END

