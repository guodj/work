
filein1 = 'ionosphere_n0000*.dat.save'
filelist1 = findfile(filein1)
nfiles1 = n_elements(filelist1)

filein2 = 'amie_only.*.save'
filelist2 = findfile(filein2)
nfiles2 = n_elements(filelist2)

mr = 40.0

vars = [7,8]
titles = ['Field Aligned Current','Potential']

nplots = n_elements(vars)

for n = 0,nfiles1-1 do begin

  filein = filelist1(n)
  print, "Reading file ", filein
  restore, filein

  if (n eq 0) then begin
    indata1  = fltarr(nplots,2,nfiles1,nlons,nlats)
    lats_n   = reform(data(0,3,*,*))
    lats_s   = 180.0 - reform(data(1,3,*,*))
    lons     = reform(data(0,4,*,*))*!pi/180.0 + !pi/2.0
  endif

  for i=0,nplots-1 do begin
    if (vars(i) gt -1) then begin
      for j=0,1 do indata1(i,j,n,*,*) = data(j,vars(i),*,*)
    endif else begin
      e1 = reform(data(*,9,*,*))/1000.0
      e2 = reform(data(*,10,*,*))/1000.0
      e3 = reform(data(*,11,*,*))/1000.0
      s  = reform(data(*,6,*,*))
      jh = s*(e1^2+e2^2+e3^3)
      for j=0,1 do begin
        indata1(i,j,n,*,*) = jh(j,*,*)
      endfor
    endelse
  endfor

endfor

for n = 0,nfiles2-1 do begin

  filein = filelist2(n)
  print, "Reading file ", filein
  restore, filein

  if (n eq 0) then begin
    indata2  = fltarr(nplots,2,nfiles2,nlons,nlats)
  endif

  for i=0,nplots-1 do begin
    if (vars(i) gt -1) then begin
      for j=0,1 do indata2(i,j,n,*,*) = data(j,vars(i),*,*)
    endif else begin
      e1 = reform(data(*,9,*,*))/1000.0
      e2 = reform(data(*,10,*,*))/1000.0
      e3 = reform(data(*,11,*,*))/1000.0
      s  = reform(data(*,6,*,*))
      jh = s*(e1^2+e2^2+e3^3)
      for j=0,1 do begin
        indata2(i,j,n,*,*) = jh(j,*,*)
      endfor
    endelse
  endfor

endfor

hems = ["Northern Hemisphere", "Southern Hemisphere"]

openw,1,'CPCP.Penetration'


for iplot = 0,nplots-1 do begin

  title = titles(iplot)

; determine color table

  mini1 = min(indata1(iplot,*,*,*,*))
  maxi1 = max(indata1(iplot,*,*,*,*))
  
  mini2 = min(indata2(iplot,*,*,*,*))
  maxi2 = max(indata2(iplot,*,*,*,*))
  
  mini = min([mini1,mini2])
  maxi = min([maxi1,maxi2])

  dmini = min(indata1(iplot,*,*,*,*)-indata2(iplot,*,*,*,*))
  dmaxi = max(indata1(iplot,*,*,*,*)-indata2(iplot,*,*,*,*))
  
  dmaxi = max([abs(dmini),dmaxi])
  dmini = -1.0*dmaxi

  for hem = 0, 0 do begin

    setdevice, 'variable'+tostr(iplot)+'_'+strmid(hems(hem),0,1)+'.ps','p',4

    if mini ge 0.0 then begin
      mini = 0.0
      ctname = getenv("IDL_EXTRAS")+"white_red.ct"
    endif else begin
      maxi = max([abs(mini),maxi])
      mini = -maxi
      ctname = getenv("IDL_EXTRAS")+"blue_white_red.ct"
    endelse
    readct,ncolors, ctname
    c_colors = (ncolors-1)*findgen(30)/29.0 + 1

    ppp = 12
    space = 0.05
    pos_space, ppp, space, sizes

    if hem eq 0 then begin
      loc = where(lats_n(0,*) le mr)
      xpos = lats_n(*,loc)*cos(lons(*,loc))
      ypos = lats_n(*,loc)*sin(lons(*,loc))
    endif else begin
      loc = where(lats_s(0,*) le mr)
      xpos = lats_s(*,loc)*cos(lons(*,loc))
      ypos = lats_s(*,loc)*sin(lons(*,loc))
    endelse

    if (strpos(titles(iplot), 'otential') gt -1) then begin
      levels = (maxi-mini)*findgen(21)/20.0 + mini 
      diff_levels = (dmaxi-dmini)*findgen(21)/20.0 + dmini
    endif else begin
      levels = (maxi-mini)*findgen(7)/6.0 + mini
      diff_levels = (dmaxi-dmini)*findgen(7)/6.0 + dmini
    endelse

    c_levels = (maxi-mini)*findgen(30)/29.0 + mini
    diff_c_levels = (dmaxi-dmini)*findgen(30)/29.0 + dmini

    for n = 0, nfiles1-1 do begin

; Plot the first variable in the first column:

      pn = (n*3) mod ppp

      if pn eq 0 then begin
        plotdumb
        xyouts, 0.0, 1.01, titles(iplot), /norm
        xyouts, 1.0, 1.01, hems(hem), alignment = 1.0, /norm
      endif

      get_position, ppp, space, sizes, pn, pos

      xyouts, pos(0)-0.05, (pos(1)+pos(3))/2.0, 'Time '+tostr(n),orient = 90, /norm, align=0.5
      if (pn lt 3) then $
	xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.05, filelist1(n), /norm, align = 0.5

      contour, indata1(iplot,hem,n,*,loc), xpos, ypos, /follow,          $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,      $
        /cell_fill, /noerase
      contour, indata1(iplot,hem,n,*,loc), xpos, ypos, /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase

      plotmlt, mr
      mini_tmp = min(indata1(iplot,hem,n,*,loc))
      maxi_tmp = max(indata1(iplot,hem,n,*,loc))
      if (abs(mini_tmp) gt 100.0) or (abs(maxi_tmp) gt 100.0) or        $
         (abs(maxi_tmp) lt 0.1) then begin
        maxs = "Max:"+string(maxi_tmp,format="(e7.1)")
        mins = "Min:"+string(mini_tmp,format="(e8.1)")
      endif else begin
        maxs = "Max:"+string(maxi_tmp,format="(f5.2)")
        mins = "Min:"+string(mini_tmp,format="(f6.2)")
      endelse
      xyouts, pos(0),pos(1)-0.02, mins, /norm, charsize = 0.8
      xyouts, pos(2),pos(1)-0.02, maxs, /norm, align=1.0, charsize = 0.8

      amie_cpcp = maxi_tmp - mini_tmp

; Plot the second variable in the second column:

      pn = (n*3 + 1) mod ppp

      get_position, ppp, space, sizes, pn, pos

      if (pn lt 3) then $
	xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.05, filelist2(n), /norm, align = 0.5

      contour, indata2(iplot,hem,n,*,loc), xpos, ypos, /follow,          $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,      $
        /cell_fill, /noerase
      contour, indata2(iplot,hem,n,*,loc), xpos, ypos, /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase

      plotmlt, mr
      mini_tmp = min(indata2(iplot,hem,n,*,loc))
      maxi_tmp = max(indata2(iplot,hem,n,*,loc))
      if (abs(mini_tmp) gt 100.0) or (abs(maxi_tmp) gt 100.0) or        $
         (abs(maxi_tmp) lt 0.1) then begin
        maxs = "Max:"+string(maxi_tmp,format="(e7.1)")
        mins = "Min:"+string(mini_tmp,format="(e8.1)")
      endif else begin
        maxs = "Max:"+string(maxi_tmp,format="(f5.2)")
        mins = "Min:"+string(mini_tmp,format="(f6.2)")
      endelse
      xyouts, pos(0),pos(1)-0.02, mins, /norm, charsize = 0.8
      xyouts, pos(2),pos(1)-0.02, maxs, /norm, align=1.0, charsize = 0.8

      amie_mike_cpcp = maxi_tmp - mini_tmp

; Plot the difference in the third column:

      pn = (n*3 + 2) mod ppp

      get_position, ppp, space, sizes, pn, pos

      if (pn lt 3) then $
	xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.05, 'Difference', /norm, align = 0.5

      difference = reform(indata1(iplot,hem,n,*,loc) - indata2(iplot,hem,n,*,loc))

      contour, difference, xpos, ypos, /follow,          $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = diff_c_levels, c_colors = c_colors,      $
        /cell_fill, /noerase
      contour, difference, xpos, ypos, /follow, levels=diff_levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase

      plotmlt, mr
      mini_tmp = min(difference)
      maxi_tmp = max(difference)
      if (abs(mini_tmp) gt 100.0) or (abs(maxi_tmp) gt 100.0) or        $
         (abs(maxi_tmp) lt 0.1) then begin
        maxs = "Max:"+string(maxi_tmp,format="(e7.1)")
        mins = "Min:"+string(mini_tmp,format="(e8.1)")
      endif else begin
        maxs = "Max:"+string(maxi_tmp,format="(f5.2)")
        mins = "Min:"+string(mini_tmp,format="(f6.2)")
      endelse
      xyouts, pos(0),pos(1)-0.02, mins, /norm, charsize = 0.8
      xyouts, pos(2),pos(1)-0.02, maxs, /norm, align=1.0, charsize = 0.8

      diff_cpcp = maxi_tmp - mini_tmp

      if (iplot eq 1) then printf,1, n, amie_cpcp, amie_mike_cpcp, diff_cpcp

    endfor

    closedevice

  endfor

endfor

save,indata1,indata2, file="patterns.save"

close,1

end    