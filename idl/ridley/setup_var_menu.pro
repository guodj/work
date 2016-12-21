
PRO Setup_Var_Menu, numoffil, datafile, timeint, titles, units, 	$
	count, nrows, rowlen, colloc, coltype, colnele

  k=0

  max = 0
  nrows = lonarr(numoffil)

  timeint = strarr(numoffil,2) 
  titles = strarr(numoffil,500) 
  units = titles
  coltype = titles
  colloc = lonarr(numoffil,500)
  colnele = lonarr(numoffil,500,3)
  count = intarr(numoffil) 
  rowlen = intarr(numoffil) 

  while k lt numoffil do begin

    openr, k+61, datafile(k)+'.hed'
    x$    = ''
    x     = strarr(150)

    i=1

    While (strpos(x(i-1),'NOTES') eq -1)  do begin

      readf,k+61,x$
      x(i) = x$
      i=i+1

    EndWhile

    nrows(k) = strmid(x(5),40,12)
    rowlen(k) = strmid(x(3),40,12)

    x=x(11:i-3)

    checker=0

    While (checker ne -1)  do begin

      readf,k+61,x$
      if (strpos(mklower(x$),'end') ne -1) then begin
        timeint(k,1)=x$
        checker=-1
      endif
                                                                        
      if (strpos(mklower(x$),'sta') ne -1) then timeint(k,0)=x$

    EndWhile

    close, k+61

    timeint(k,0)=strmid(timeint(k,0),14,25)
    timeint(k,1)=strmid(timeint(k,1),14,25)

    titles(k,0:i-14)=strmid(x,7,14)
    units(k,0:i-14)=strmid(x,21,14)	
    coltype(k,0:i-14)=mklower(strmid(x,57,1))	
    colloc(k,0:i-14)=long(strmid(x,70,8))	

    for j=0,i-14 do begin
      temps = strmid(x(j),58,12)
      loc = strpos(temps,' ')
      if loc eq 0 then begin
        colnele(k,j,0) = 1
        colnele(k,j,1) = 0
        colnele(k,j,2) = 0
      endif else begin
        colnele(k,j,0) = long(temps)
        loc = strpos(temps,'x')
	if (loc eq -1) then begin
          colnele(k,j,1) = 0
          colnele(k,j,2) = 0
        endif else begin
          temps = strmid(temps,loc+1,strlen(temps)-loc-1)
          colnele(k,j,1) = long(temps)
          loc = strpos(temps,'x')
	  if (loc eq -1) then begin
            colnele(k,j,2) = 0
	  endif else begin
            temps = strmid(temps,loc+1,strlen(temps)-loc-1)
            colnele(k,j,1) = long(temps)
          endelse
        endelse
      endelse
    endfor

    count(k)=i-14

    if count(k) gt max then max = count(k)

    k=k+1

  EndWhile

  titles=titles(*,0:max)
  units=units(*,0:max)

  RETURN

END

