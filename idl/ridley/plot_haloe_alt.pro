ctname = '/f/ridley/d.data/d.pem/d.climate/d.data/jim6.ct'

filein = 'haloe9_10.mag.save'
psfile = ask('ps file name','test.ps')

setdevice,psfile,'l',4,0.9

restore,filein

nlats = 45

dlat  = 180.0/float(nlats)

ppp = 2
space = 0.1
pos_space, ppp, space, sizes, ny=1
dy = float(!d.y_ch_size)/float(!d.y_size)
readct,ncolors,ctname
clevels = float(ncolors)*(findgen(30)+0.5)/30.0

height_loc = where(heights ge 80.0 and heights le 160.0,height_count)
lats = findgen(nlats+1)*dlat - 90.0

outbin = fltarr(nlats+1,height_count)

minlat = -90.0

for pic = 0,1 do begin

  k = 3

  if pic eq 0 then begin
    lat = latdata
  endif
  if pic eq 1 then begin
    lat = mlatdata
  endif

  outbin(*,*) = 0.0

  for iy=0,nlats do begin
    lat1 = minlat + float(iy)*dlat
    lat2 = minlat + float(iy+1)*dlat
    for ix=0,height_count-1 do begin
      loc = where( (lat ge lat1) and				$
	           (lat lt lat2)	and				$
	           (vdata(*,height_loc(ix)) ne 0.0),count)
        if count gt 0 then 					$
	  if count gt 1 then 					$
	    outbin(iy,ix) = mean(vdata(loc,height_loc(ix)))		$
	  else outbin(iy,ix) = vdata(loc(0),height_loc(ix))
    endfor

  endfor

mini = 1.0e-7
maxi = 3.0e-4

  dr = maxi-mini

  levels = (maxi-mini)*findgen(30)/29.0 + mini

;  for j=0,1 do begin
j = 1
    pn = pic

    get_position,ppp,space,sizes,pn,pos

    if pic eq 1 then begin
      pos([0,2]) = pos([0,2]) + 0.1
    endif

    pos([1,3]) = pos([1,3]) - 0.1

    temp = outbin
    locn = where(temp lt mini+0.01*dr and temp gt 0.0,countn)
    if countn gt 0 then temp(locn) = mini+0.01*dr
    locp = where(temp gt maxi,countp)
    if countp gt 0 then temp(locp) = maxi

    contour, temp, lats,heights(height_loc), /follow, nlevels = 30,    $
	xstyle = 1,		$
	ystyle = 1, /noerase, pos = pos, /cell_fill, $
	c_colors=clevels,levels=levels,min_value = mini,		$
	xtitle = 'Latitude', ytitle = 'Altitude'

      specheight = where(heights(height_loc) gt 130.0)
      alt = heights(height_loc(specheight(0)))
      specheight = height_loc(specheight(0))-height_loc(0)

      oplot, [-90,90],[alt,alt],linestyle=1

      dp = (pos(2) - pos(0))/8.0
      posct = [pos(0)+dp,pos(1)-dy*5.0,pos(2)-dp,pos(1)-dy*3.5]
      plot, [0,ncolors], [0,1], /noerase, pos = posct, 			$
	    xstyle=5,ystyle=5, /nodata
      xct = [0.0,0.0,1.0,1.0,0.0]
      yct = [0.0,1.0,1.0,0.0,0.0]
      index = indgen(ncolors)
      for i=0,ncolors-1 do                                      $
        polyfill, xct+float(index(i)), yct, color = i

      lmax = 10.0^(float(fix(alog10(maxi)))-1.0)
      lmaxs = tostr(alog10(lmax))
      lmaxs = '(x10!E'+lmaxs+'!N)'
      plot, [mini/lmax,maxi/lmax],[0,1],/noerase,/nodata,pos=posct,	$
	    xstyle = 1, ystyle = 5, xtitle = lmaxs, xticklen=-0.2
      plots, [mini/lmax,mini/lmax],[0,1]
      plots, [maxi/lmax,maxi/lmax],[0,1]

      pos(1) = pos(3)+0.01
      pos(3) = pos(1)+0.3

      plot, lats, temp(*,specheight), ytitle = 'NO MMR', xstyle = 1,	$
	yrange = [mini,maxi], pos = pos, /noerase, 			$
	xtickname = strarr(10)+' ', min_value = mini

        top = pos(3) + dy/2.0
	if pic eq 0 then 					$
          xyouts, mean(pos([0,2])), top, 'Geographic', /norm,		$
	    charsize=1.3, alignment = 0.5
	if pic eq 1 then 					$
          xyouts, mean(pos([0,2])), top, 'Magnetic', /norm,		$
	    charsize=1.3, alignment = 0.5

endfor

closedevice

end
