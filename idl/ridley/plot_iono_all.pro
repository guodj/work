
filein = 'ionosphere_n00000*.dat.save'
filein = ask('filename',filein)
line = ''
filelist = findfile(filein)
nfiles = n_elements(filelist)

mr = float(ask('maximum range','40'))

vars = [5,7,8]
titles = ['Hall Conductance', 'Field Aligned Current','Potential']

nplots = n_elements(vars)

indices = fltarr(2,5,nfiles)

Re = 6371.2*1000.0

for n = 0,nfiles-1 do begin

  filein = filelist(n)
  print, "Reading file ", filein
  restore, filein

  if (n eq 0) then begin
    indata = fltarr(nplots,2,nfiles,nlons,nlats)
    lats_n   = reform(data(0,3,*,*))
    lats_s   = 180.0 - reform(data(1,3,*,*))
    lons     = reform(data(0,4,*,*))*!pi/180.0 + !pi/2.0
  endif

  for i=0,nplots-1 do begin
    if (vars(i) gt -1) then begin
      for j=0,1 do indata(i,j,n,*,*) = data(j,vars(i),*,*)
    endif else begin
      e1 = reform(data(*,9,*,*))/1000.0
      e2 = reform(data(*,10,*,*))/1000.0
      e3 = reform(data(*,11,*,*))/1000.0
      s  = reform(data(*,6,*,*))
      jh = s*(e1^2+e2^2+e3^3)
      for j=0,1 do begin
        indata(i,j,n,*,*) = jh(j,*,*)
        indices(j,0,n) = 0.0
        for k=0,nlons-2 do for l=1,nlats-1 do begin
          dr = (data(j,3,k,l) - data(j,3,k,l-1))*!pi/180.0
          dt = (data(j,4,k,l) - data(j,4,k+1,l))*!pi/180.0
          ave_lat = data(j,3,k,l)*!pi/180.0
          dt = dt*cos(ave_lat)
          area = dr*(Re) * dt*(Re)
          indices(j,0,n) = indices(j,0,n)+ area*jh(j,k,l)
        endfor
      endfor
    endelse
  endfor

  for j=0,1 do begin
    indices(j,1,n) = min(data(j,8,*,*))
    indices(j,2,n) = max(data(j,8,*,*))
    indices(j,3,n) = indices(j,2,n) - indices(j,1,n)
  endfor

endfor

indices(*,0,*) = indices(*,0,*)/1.0e9

hems = ["Northern Hemisphere", "Southern Hemisphere"]

if strpos(filelist(0),'3000') gt 0 then time = findgen(nfiles)  $
else time = findgen(nfiles)+1

save, indices, time, filename="indices.save"

for hem = 0, 1 do begin

  setdevice, 'indices'+'_'+strmid(hems(hem),0,1)+'.ps','p',4
  ppp = 2
  space = 0.05
  pos_space, ppp, space, sizes, ny = ppp

  plotdumb

  get_position, ppp, space, sizes, 0, pos, /rect
  plot, time, indices(hem,3,*), ytitle = 'Cross PC Potential (kV)', $
        /noerase, pos = pos, xrange = [0,60], yrange = [0,100]

  get_position, ppp, space, sizes, 1, pos, /rect
  plot, time, indices(hem,0,*), ytitle = 'Global Joule Heating (GW)', $
        /noerase, pos = pos, xrange = [0,60]

  closedevice

endfor

for iplot = 0,nplots-1 do begin

  title = titles(iplot)

  print, "Plotting variable ",title

; determine color table

  mini = min(indata(iplot,0,*,*,*))*1.05
  maxi = max(indata(iplot,0,*,*,*))*1.05
  
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

    ppp = min([nfiles,9])
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

    if iplot eq nplots-1 then begin
      levels = (maxi-mini)*findgen(21)/20.0 + mini
    endif else begin
      levels = (maxi-mini)*findgen(7)/6.0 + mini
    endelse

    print, n_elements(levels)

    c_levels = (maxi-mini)*findgen(30)/29.0 + mini

    for n = 0, nfiles-1 do begin

      pn = n mod ppp

      if pn eq 0 then begin
        plotdumb
        xyouts, 0.0, 1.01, titles(iplot), /norm
        xyouts, 1.0, 1.01, hems(hem), alignment = 1.0, /norm
      endif

      get_position, ppp, space, sizes, pn, pos

      contour, indata(iplot,hem,n,*,loc), xpos, ypos, /follow,          $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,      $
        /cell_fill, /noerase
      contour, indata(iplot,hem,n,*,loc), xpos, ypos, /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase

      plotmlt, mr
      mini_tmp = min(indata(iplot,hem,n,*,loc))
      maxi_tmp = max(indata(iplot,hem,n,*,loc))
      if (abs(mini_tmp) gt 100.0) or (abs(maxi_tmp) gt 100.0) or        $
         (abs(maxi_tmp) lt 0.1) then begin
        maxs = "Ma:"+string(maxi_tmp,format="(e8.2)")
        mins = "Mi:"+string(mini_tmp,format="(e9.2)")
      endif else begin
        maxs = "Ma:"+string(maxi_tmp,format="(f5.2)")
        mins = "Mi:"+string(mini_tmp,format="(f6.2)")
      endelse
      xyouts, pos(0),pos(1)-0.02, mins, /norm, charsize=0.8
      xyouts, pos(2),pos(1)-0.02, maxs, /norm, align=1.0,charsize=0.8

    endfor

    closedevice

  endfor

endfor

end    