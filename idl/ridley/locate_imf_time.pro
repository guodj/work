
common ffinfo, header

day = 1
maxday = 31

year = 0
read, 'Enter year : ',year
if year gt 1900 then year = year - 1900
if year gt 100 then year = year - 100

read, 'Enter starting Day : ',day
read, 'Enter ending day   : ',maxday

outfile = ''
read, 'Enter output file name : ',outfile

openw,1, outfile

missing = -1.0e32

filename = 'imf_'+chopr('0'+tostr(year),2)

col_scal = [0,1,2,13]

while day le maxday do begin

  stime = [year,01,day,00,00,00]
  etime = [year,01,day+1,00,00,00]

  print, 'Reading data from file '+filename+'...'

  read_flat_scalor, stime, etime, col_scal, time, data_scal, nrows,	$
	filename = filename

  nrows = nrows(0)

  loc = where(data_scal(3,*) ne missing,count)

  ppp = 3
  space = 0.03

  if count gt nrows/2 then begin

    window,0,xsize=1000,ysize=800

    plotdumb

    pos_space,ppp,space,sizes, nx = 1

    pos = fltarr(ppp+1,4)

    for i=0,ppp-1 do begin
      get_position, ppp, space, sizes, i, tpos, 			$
                  xmargin = 0.05, ymargin=0.09, /rect
      pos(i,*) = tpos(*)
    endfor

    pos(ppp,*) = [pos(0,0),0.0,pos(0,2),0.05]

    data = fltarr(nrows)

    for col = 0, 2 do begin

      data(*) = data_scal(col,*)
      plot, data,							$
	  ystyle = 1, xstyle = 1,					$
	  /noerase, pos = pos(col,*), yrange = [-20,20],		$
	  max_value = 20.0, min_value=-20.0

      oplot, data*0.0, linestyle = 2

    endfor

    plot, [0,1], /nodata, /noerase, xstyle = 5, ystyle = 5, pos = pos(ppp,*)
    oplot, [0.01,0.01,0.99,0.99,0.01], [0.01,0.99,0.99,0.01,0.01]
    xyouts, 0.1, 0.5, 'Day Number : '+tostr(day)
    xyouts, 0.5, 0.5, 							$
        'Press Here to move to next time (second button to end)', 	$
        alignment=0.5

    command = 0

    basetime = min(time)
    picked   = 0
    utsec    = fltarr(100)

    while ((command+1) mod (ppp+1)) ne 0 do begin

      get_mouse_command, pos, command, ret_pos

      if command ne ppp*2+1 then command = command mod (ppp + 1)

      if command lt ppp then begin

        plot, [0,1], /nodata, /noerase, 				$
            xstyle = 5, ystyle = 5, pos = pos(command,*)

        oplot, [ret_pos(0), ret_pos(0)], [-100,100]

        c_r_to_a, itime, time(command,nrows*ret_pos(0))
        print, itime

; want to print time out to the file, but first we want to sort them 
;   so they are increasing.
; want to store the column which we picked also. This will be saved
;   as the "seconds" - since the data is 2 minute data, the seconds
;   are irrelavent, just have to remember to subtract them later.

        utsec(picked) = time(command,nrows*ret_pos(0)) - basetime + 	$
                      float(command)
        picked = picked + 1

      endif

    endwhile

  endif else begin

    command = ppp
    picked = 0

  endelse
  
; if we select next time, but not exit, display a message

  if command lt ppp*2 then begin
    day = day + 1 
    if (day le maxday) then begin
      plotdumb
      xyouts, 0.5, 0.5, 'Reading new data. Please wait...', /norm,	$
              alignment = 0.5
    endif
  endif else day = maxday + 1

; at the end of every day, sort the selected times and write them to the file

  if picked gt 0 then begin

    utsec = utsec(0:picked-1)
    sort, utsec
    for i=0,picked-1 do begin
      c_r_to_a, itime, basetime + double(utsec(i))
      printf,1, itime
    endfor

  endif

endwhile

close,1
wdelete

end
