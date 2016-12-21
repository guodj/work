
pro izmempot

  print, 'Enter group of files (ex *-north.pot) :'
  filein = ''
  read, filein

  filelist = findfile(filein)
  count = n_elements(filelist)

  data = fltarr(count,33,24) 
  time = strarr(count)
  date = strarr(count)
  imf = strarr(count,3)

  for i=0,count-1 do begin

    data1 = read_izmem(filelist(i))
    data(i,*,*) = data1.pot
    date(i) = data1.date
    time(i) = data1.time
    imf(i,*) = data1.imf

  endfor

  az  = fltarr(33,24)
  lat = fltarr(33,24)
  for i=0,32 do az(i,*) = 180-(findgen(24)+0.5)*360.0/24.0 
  for i=0,23 do lat(*,i) = findgen(33)+1.0

  fileout = ''

  que$ = ''
  read,'Subtract current from previous (y/n) ? ', que$
  if (strmid(que$,0,1) eq 'y') or (strmid(que$,0,1) eq 'Y') then 	$
    sub = 1 else sub = 0
  read, 'Enter ps file name (return for screen) : ', fileout

  lats = 33
  lons = 24

  condata = fltarr(3,lats*lons)

  maxran = 41.0

  if strlen(fileout) gt 0 then setdevice, fileout,l,4		$
  else window,1,xsize=1000,ysize=800

  plot, [0,1], /nodata, xstyle=5, ystyle=5, 			$
	Title = 'IZMEM Potentials', charsize = 1.5, pos=[0,0,1,1]

  space = 0.06
  ppp = 4
  pos_space, ppp, space, sizes

  aaa = !pi*findgen(361.0)/180.0
  xp = sin(aaa)
  yp = cos(aaa)
  maxran = 40.0

  for i=0,lats-1 do for j=0,lons-1 do begin
    condata(0,i*lons+j) = az(i,j)
    condata(1,i*lons+j) = lat(i,j)
  endfor

  if sub eq 0 then begin 
    mmdata = mm(data)
    mmdata(0) = fix(mmdata(0)/30.0 - 1.0) * 30.0
    mmdata(1) = fix(mmdata(1)/30.0 + 1.0) * 30.0
    if mmdata(1)-mmdata(0) lt 30.0 then mmdata(1) = mmdata(0) + 30.0
    levels = fix(findgen(30)*(mmdata(1)-mmdata(0))/29 + mmdata(0))
  endif else levels = 200.0*(findgen(24)+1.0)/25.0 - 100.0

  for n=sub,count-1 do begin
;  for n=0,count-1 do begin

    ln = n-1
;    ln = 0
    pn = n mod ppp

    get_position, count, space, sizes, pn, pos

    plot, maxran*xp, maxran*yp, xstyle=5, ystyle=5, pos=pos, /noerase
    oplot, 30.0*xp, 30.0*yp,linestyle=2
    oplot, 20.0*xp, 20.0*yp,linestyle=2
    oplot, 10.0*xp, 10.0*yp,linestyle=2
    oplot, [0.0,0.0],[-maxran,maxran]
    oplot, [-maxran,maxran], [0.0,0.0]

    xyouts, -1.0*maxran, 0.95*maxran, date(n), alignment = 0.0, charsize=0.9
    xyouts, -1.0*maxran, (-1.08)*maxran, imf(n,2), 		$
	alignment = 0.0, charsize=0.9
    xyouts, -1.0*maxran, (-1.0)*maxran, imf(n,1), 		$
	alignment = 0.0, charsize=0.9
    xyouts, -1.0*maxran, (-0.92)*maxran, imf(n,0), 		$
	alignment = 0.0, charsize=0.9

    if sub eq 0 then 						$
      xyouts, 1.0*maxran, 0.95*maxran, time(n), 		$
	alignment = 1.0, charsize=0.9				$
    else							$
      xyouts, 1.0*maxran, 0.95*maxran, strmid(time(n),0,5)+'-'+time(ln),  $
	alignment = 1.0, charsize=0.9

    xyouts, 0.0,maxran*(-1.08), '00', alignment=0.5, charsize=0.8
    xyouts, maxran*(1.01), -0.025*maxran, '06', charsize=0.8
    xyouts, 0.0,maxran*(1.02), '12', alignment=0.5, charsize=0.8
    xyouts, maxran*(-1.01), -0.025*maxran, '18', alignment=1.0, charsize=0.8

    for i=0,lats-1 do for j=0,lons-1 do begin
      if sub eq 1 then 							$
        condata(2,i*lons+j) = 						$
		100*(data(n,i,j)-data(ln,i,j))/				$
		(max(data(ln,*,*))-min(data(ln,*,*))) 			$
      else								$
        condata(2,i*lons+j) = data(n,i,j)
    endfor

    basic_contour, condata, maxran, pos, levels

    xyouts,maxran*(1.08),maxran*(-1.08), 				$
	'Max : '+tostr(max(condata(2,*))), alignment=1.0
    xyouts,maxran*(1.08),maxran*(-1.0),	 				$
	'Min : '+tostr(min(condata(2,*))), alignment=1.0

    if pn+1 eq ppp then begin
      if !d.name eq 'X' then prompt_for_next
      plot, [0,1], /nodata, xstyle=5, ystyle=5, 			$
	Title = 'IZMEM Potentials', charsize = 1.5, pos=[0,0,1,1]
    endif

  endfor

  if !d.name eq 'PS' then begin
    device, /close
    set_plot, 'X'
  endif

end

