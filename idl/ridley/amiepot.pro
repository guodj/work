pro amiepot

  filein = ''
  fileti = ''
  amiefile = ''
  fileout = ''
  line = ''

  read, 'Enter file name (suffix only) of the AMIE output : ', amiefile
  filein = amiefile+'.dat'
  fileti = amiefile+'.time'
  que$ = ''
  read,'Subtract current from previous (y/n) ? ', que$
  if (strmid(que$,0,1) eq 'y') or (strmid(que$,0,1) eq 'Y') then 	$
    sub = 1 else sub = 0
  read, 'Enter ps file name (return for screen) : ', fileout

  close,1
  openr,1,filein

  readf,1, line
  readf,1, line
  lons = fix(strmid(line,13,3))

  for i=0,7 do readf,1,line
  readf,1, line
  lats = fix(strmid(line,12,3))
  for i=0,5 do readf,1,line
  ltpos = fltarr(lons,lats)
  lnpos = fltarr(lons,lats)
  for i=1,lons do for j=1,lats do begin
    ltpos(i-1,j-1) = float(j-1)*(90.0-50.0)/float(lats-1)
    lnpos(i-1,j-1) = 180.0-float(i)*360.0/float(lons)
  endfor

  loc = where(lnpos lt 0.0, count)

  if count gt 0 then lnpos(loc) = lnpos(loc) + 360.0

  data = fltarr(1,lons,lats)
  fakedata = fltarr(1,6)
  count = 0
  while (not eof(1)) do begin

    readf,1,line
    print, line
    if count gt 0 then data = [data, fltarr(1,lons,lats)]

    for i=0,24 do for j=1,36,6 do begin
      readf,1, format = '(6E13.5)',fakedata
      data(count,(j-1):(j+4),i) = fakedata(0,*)
    endfor

    count = count + 1

  endwhile

  close,1,2
  openr,2,fileti

  readf,2, line
  readf,2, line
  year = strmid(line,1,4)
  mon = 'JanFebMarAprMayJunJulAugSepOctNovDec'
  date = strarr(count)
  time = strarr(count)
  for i=0,count-1 do begin
    readf,2,line
    month = fix(strmid(line,1,2))
    date(i) = strmid(mon,(month-1)*3,3)+' '+strmid(line,3,2)+', '+year
    time(i) = strmid(line,6,2)+':'+strmid(line,8,2)+' UT'
  endfor

  close, 2

  condata = fltarr(3,lats*lons)

  maxran = 41.0

  if strlen(fileout) gt 0 then setdevice, fileout,'p',4,0.83	$
  else window,1,xsize=1000,ysize=800

  plot, [0,1], /nodata, xstyle=5, ystyle=5

  space = 0.01
  ppp = 9
  pos_space, ppp, space, sizes

  aaa = !pi*findgen(361.0)/180.0
  xp = sin(aaa)
  yp = cos(aaa)
  maxran = 40.0

  for i=0,lons-1 do for j=0,lats-1 do begin
    condata(0,i*lats+j) = lnpos(i,j)
    condata(1,i*lats+j) = ltpos(i,j)
  endfor

  data(22,*,*) = (data(22,*,*))/1.0

   ntemp = 22

;  for n=sub,count-1 do begin
  for n=22,28 do begin
        ;23
    ln = n-1
    ln = ntemp
;    pn = (n-15) mod ppp

    if n gt ntemp then pn = (n-ntemp+2) mod ppp else pn = 1

    get_position, count, space, sizes, pn, pos

    if !d.name eq 'PS' then begin
      if (pn eq 0) or (n eq ntemp) then yadd = 1.0-pos(3)
      pos([1,3]) = pos([1,3])+yadd
    endif

    if n eq ntemp then begin

      plot, [0,1], /nodata, /noerase, pos = [0,0,1,1], 		$
	xstyle=5, ystyle=5
      xyouts, pos(0)-0.04, pos(3)-0.015, '(A)', charsize = 1.2
      xyouts, pos(0) - (pos(2)-pos(0)+space), pos(1)-space/2.0+0.005,'(B)', $
	charsize=1.2
      oplot, [pos(0)-(pos(2)-pos(0)+space)+0.05,pos(2)+(pos(2)-pos(0)+space)], $
	[pos(1)-space/2.0,pos(1)-space/2.0]

    endif

    plot, maxran*xp, maxran*yp, xstyle=5, ystyle=5, pos=pos, /noerase

    if (pn eq 0) or ((ntemp gt 0) and (n eq ntemp+1)) then begin
      t = 2.0*!pi*findgen(91)/360.0
      phi = (120.0+(n)*1.25)*!pi/180.0
      x = [0.0,maxran*cos(t+phi),0.0]
      y = [0.0,maxran*sin(t+phi),0.0]
      if !d.name eq 'PS' then c = 220 else c = 55
      polyfill,x,y, color=c
    endif

    oplot, maxran*xp, maxran*yp
    oplot, 30.0*xp, 30.0*yp,linestyle=2
    oplot, 20.0*xp, 20.0*yp,linestyle=2
    oplot, 10.0*xp, 10.0*yp,linestyle=2
    oplot, [0.0,0.0],[-maxran,maxran]
    oplot, [-maxran,maxran], [0.0,0.0]

    if pn eq 1 then 						$
    xyouts, 0.5, 1.03, date(n), alignment = 0.5, charsize=1.2, /norm

    if n gt ntemp then begin
      xyouts, -1.0*maxran, 0.90*maxran, time(n), 		$
        alignment = 0.0, charsize=0.9				
    endif else begin
      xyouts,-0.95*maxran, 0.90*maxran, 'Steady', 		$
	alignment = 0.0, charsize=0.9				
    endelse

    if (sub eq 1) and (n gt ntemp) then begin
      xyouts, 0.95*maxran, 0.90*maxran, 'Residual',  	$
	alignment = 1.0, charsize=0.7
;      xyouts, 1.0*maxran, 1.02*maxran, strmid(time(n),0,5)+'-'+time(ln),  $
;	alignment = 1.0, charsize=0.9
    endif

    if pos(1) lt 0.3 then 					$
      xyouts, 0.0,maxran*(-1.1), '00', alignment=0.5, charsize=0.8
    if pos(2) gt 0.7 then 					$
      xyouts, maxran*(1.01), -0.05*maxran, '06', charsize=0.8
    if pos(3) gt 0.9 then 					$
      xyouts, 0.0,maxran*(1.02), '12', alignment=0.5, charsize=0.8
    if pos(0) lt 0.3 then 					$
      xyouts, maxran*(-1.01), -0.05*maxran, '18', alignment=1.0, charsize=0.8

    if (sub eq 0) or (n eq ntemp) then begin 
      mmdata = mm(data)/1000.0
      mmdata(0) = fix(mmdata(0)/30.0 - 1.0) * 30.0
      mmdata(1) = fix(mmdata(1)/30.0 + 1.0) * 30.0
      if mmdata(1)-mmdata(0) lt 30.0 then mmdata(1) = mmdata(0) + 30.0
      levels = fix(findgen(30)*(mmdata(1)-mmdata(0))/29 + mmdata(0))
    endif else levels = 200.0*(findgen(24)+1.0)/25.0 - 100.0

    for i=0,lons-1 do for j=0,lats-1 do begin
      if (sub eq 1) and (n gt ntemp) then 				$
        condata(2,i*lats+j) = 						$
		100*(data(n,i,j)-data(ln,i,j))/				$
		(max(data(ln,*,*))-min(data(ln,*,*))) 			$
      else								$
        condata(2,i*lats+j) = data(n,i,j)/1000.0
    endfor

    basic_contour_new, condata, maxran, pos, levels

    if max(condata(2,*)) lt 100 then begin
      xyouts,maxran*(1.0),maxran*(-0.96), 				$
	  'Ma: '+tostr(max(condata(2,*))), alignment=1.0, charsize=0.85
    endif else begin
      xyouts,maxran*(1.0),maxran*(-0.96), 				$
	  'Ma:'+tostr(max(condata(2,*))), alignment=1.0, charsize=0.85
    endelse

    if min(condata(2,*)) gt -100 then begin
      xyouts,maxran*(-1.0),maxran*(-0.96),	 			$
	'Mi: '+tostr(min(condata(2,*))), alignment=0.0, charsize=0.85
    endif else begin
      xyouts,maxran*(-1.005),maxran*(-0.96),	 			$
	'Mi:'+tostr(min(condata(2,*))), alignment=0.0, charsize=0.85
    endelse

    if pn+1 eq ppp then begin
      if !d.name eq 'X' then prompt_for_next
      plot, [0,1], /nodata, xstyle=5, ystyle=5
    endif

  endfor

  if !d.name eq 'PS' then begin
    device, /close
    set_plot, 'X'
  endif

  return

end

