
filein = 'ionosphere_n000000.dat'
filein = ask('filename',filein)
line = ''

mr = float(ask('maximum range','40'))

if strpos(filein,'save') gt -1 then begin
  restore, filein
endif else begin

  openr,1,filein
  i_equals = -1
  n = 0
  while (i_equals lt 0) do begin
    readf,1,line
    ipos = strpos(line,"I=")
    if (ipos gt -1) then begin
      nlats = fix(strmid(line,ipos+2,4))
      nlons = fix(strmid(line,strpos(line,'J=')+2,5))
      i_equals = n
    endif
    n = n + 1
  endwhile
  close,1

  openr,1,filein

  for i=0,i_equals do readf,1,line

  point_lun, -1, position

  readf,1,line
  if (strpos(line,'DT=') eq -1) then point_lun, 1, position

  nvars = 18

  tmp = fltarr(nvars)
  data = fltarr(2,nvars,nlons,nlats)

  for n=0,1 do begin
    for j=0,nlons-1 do for i=0,nlats-1 do begin
      readf,1,tmp
      data(n,*,j,i) = tmp
    endfor
    if n eq 0 then begin
      readf,1,line
      point_lun, -1, position
      readf,1,line
      if (strpos(line,'I=') eq -1) then begin
        point_lun, 1, position
      endif else readf,1,line
    endif
  endfor
  close,1

endelse

hems = ["Northern Hemisphere", "Southern Hemisphere"]

for hem = 0, 1 do begin

  setdevice, filein+tostr(hem)+'.ps','p',4
  ppp = 4
  space = 0.05
  pos_space, ppp, space, sizes

  if hem eq 0 then begin
    loc = where(reform(data(0,3,0,*)) le mr)
    rang = reform(data(0,3,*,loc))
    lons = reform(data(0,4,*,loc))*!pi/180.0 + !pi/2
    xpos = rang*cos(lons)
    ypos = rang*sin(lons)
  endif else begin
    loc = where(reform(data(1,3,0,*)) gt 180.0-mr)
    rang = 180.0-reform(data(1,3,*,loc))
    lons = reform(data(0,4,*,loc))*!pi/180.0 + !pi/2
    xpos = rang*cos(lons)
    ypos = rang*sin(lons)
  endelse

  readct,ncolors, getenv("IDL_EXTRAS")+"white_red.ct"

  mini = min(data(hem,5:6,*,loc))/1.05
  maxi = max(data(hem,5:6,*,loc))*1.05

  mini = 0.0

  levels = (maxi-mini)*findgen(9)/8.0 + mini
  c_levels = (maxi-mini)*findgen(30)/29.0 + mini
  c_colors = (ncolors-1)*findgen(30)/29.0 + 1

  get_position, ppp, space, sizes, 0, pos
  pos([1,3]) =  pos([1,3]) + 0.075

  contour, data(hem,5,*,loc), xpos, ypos, /follow, nlevels=30, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,/cell_fill
  contour, data(hem,5,*,loc), xpos, ypos, /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase
  plotmlt, mr
  mini = min(data(hem,5,*,loc))
  maxi = max(data(hem,5,*,loc))
  maxs = "Max:"+string(maxi,format="(f5.2)")
  mins = "Min:"+string(mini,format="(f5.2)")
  xyouts, pos(0),pos(1)-0.02, mins, /norm
  xyouts, pos(2),pos(1)-0.02, maxs, /norm, align=1.0
  xyouts, (pos(0)+pos(2))/2.0,pos(3)+space/2.0,"Hall Conductance", 	$
	align=0.5, /norm, charsize = 1.25

  xyouts, 0.5, 1.01, hems(hem), alignment = 0.5, /norm

  get_position, ppp, space, sizes, 1, pos
  pos([1,3]) =  pos([1,3]) + 0.075

  contour, data(hem,6,*,loc), xpos, ypos, /follow, nlevels=30, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,	$
	/cell_fill, /noerase
  contour, data(hem,6,*,loc), xpos, ypos, /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase
  plotmlt, mr
  mini = min(data(hem,6,*,loc))
  maxi = max(data(hem,6,*,loc))
  maxs = "Max:"+string(maxi,format="(f5.2)")
  mins = "Min:"+string(mini,format="(f5.2)")
  xyouts, pos(0),pos(1)-0.02, mins, /norm
  xyouts, pos(2),pos(1)-0.02, maxs, /norm, align=1.0
  xyouts, (pos(0)+pos(2))/2.0,pos(3)+space/2.0,"Pedersen Conductance", 	$
	align=0.5, /norm, charsize = 1.25



readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"

maxi = max(abs(data(hem,7,*,loc)))*1.05
mini = -maxi
levels = (maxi-mini)*findgen(9)/8.0 + mini
c_levels = (maxi-mini)*findgen(30)/29.0 + mini
c_colors = (ncolors-1)*findgen(30)/29.0 + 1

get_position, ppp, space, sizes, 2, pos
contour, data(hem,7,*,loc), xpos, ypos, /follow, nlevels=30, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,	$
	/cell_fill, /noerase
contour, data(hem,7,*,loc), xpos, ypos, /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase
plotmlt, mr
mini = min(data(hem,7,*,loc))
maxi = max(data(hem,7,*,loc))
maxs = "Max:"+string(maxi,format="(e9.2)")
mins = "Min:"+string(mini,format="(e9.2)")
xyouts, pos(0),pos(1)-0.02, mins, /norm
xyouts, pos(2),pos(1)-0.02, maxs, /norm, align=1.0
  xyouts, (pos(0)+pos(2))/2.0,pos(3)+space/2.0,"Field Aligned Current", $
	align=0.5, /norm, charsize = 1.25


maxi = max(abs(data(hem,8,*,loc)))*1.05
mini = -maxi
levels = (maxi-mini)*findgen(9)/8.0 + mini
c_levels = (maxi-mini)*findgen(30)/29.0 + mini

  get_position, ppp, space, sizes, 3, pos
  contour, data(hem,8,*,loc), xpos, ypos, /follow, nlevels=30, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,	$
	/cell_fill, /noerase
  contour, data(hem,8,*,loc), xpos, ypos, /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase
  plotmlt, mr
  mini = min(data(hem,8,*,loc))
  maxi = max(data(hem,8,*,loc))
  maxs = "Max:"+string(maxi,format="(f7.2)")
  mins = "Min:"+string(mini,format="(f7.2)")
  xyouts, pos(0),pos(1)-0.02, mins, /norm
  xyouts, pos(2),pos(1)-0.02, maxs, /norm, align=1.0
  xyouts, (pos(0)+pos(2))/2.0,pos(3)+space/2.0,"Potential", 	$
	align=0.5, /norm, charsize = 1.25

  xyouts, -0.01, -0.01, filein, /norm, charsize = 0.5

closedevice

endfor

end    